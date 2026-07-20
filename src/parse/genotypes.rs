use std::collections::HashMap;
use rust_htslib::bcf::header::HeaderView;
use mongodb::bson::{doc, Bson, Document};
use rust_htslib::bcf::Record;
use rust_htslib::bcf::record::GenotypeAllele;
use crate::models::variant::VariantCategory;

use crate::models::sample::SampleInfo;

/// Validate that provided sample IDs match the samples present in the VCF header.
///
/// Checks that each requested sample exists at the expected VCF position.
pub fn validate_sample_mapping(
    header: &HeaderView,
    sample_mapping: &HashMap<String, SampleInfo>,
) -> Result<(), String> {
    let vcf_samples = header.samples();

    for (sample_id, sample_info) in sample_mapping {
        let index = sample_info.vcf_index;

        let vcf_sample = vcf_samples
            .get(index)
            .ok_or_else(|| {
                format!(
                    "Sample '{}' has VCF index {}, but the VCF contains only {} samples",
                    sample_id,
                    index,
                    vcf_samples.len()
                )
            })?;

        let vcf_sample = String::from_utf8_lossy(vcf_sample);

        if vcf_sample != *sample_id {
            return Err(format!(
                "Sample mismatch: CLI sample '{}' is assigned to VCF index {}, but VCF contains '{}'",
                sample_id,
                index,
                vcf_sample
            ));
        }
    }

    Ok(())
}


/// Parse genotype calls for selected samples.
pub fn parse_genotypes(
    record: &Record,
    sample_mapping: &HashMap<String, SampleInfo>,
    category: VariantCategory,
) -> Vec<Document> {
    let mut genotypes = Vec::new();

    for (sample_id, sample_info) in sample_mapping {
        let pos = sample_info.vcf_index;

        genotypes.push(parse_genotype(
            record,
            sample_id,
            &sample_info.display_name,
            pos,
            category
        ));
    }

    genotypes
}

/// Parse genotype information for a single sample.
/*
    SV specific format fields
    ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of paired-ends that support the event">
    ##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-ends that support the event">
    ##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts for the SV">
    ##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Raw high-quality read counts for the left control region">
    ##FORMAT=<ID=RCR,Number=1,Type=Integer,Description="Raw high-quality read counts for the right control region">
    ##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
    ##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
    ##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads that support the event">
    ##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number genotype for imprecise events">

    STR specific format fields
    ##FORMAT=<ID=LC,Number=1,Type=Float,Description="Locus coverage">
    ##FORMAT=<ID=REPCI,Number=1,Type=String,Description="Confidence interval for REPCN">
    ##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Number of repeat units spanned by the allele">
    ##FORMAT=<ID=SO,Number=1,Type=String,Description="Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat">
    ##FORMAT=<ID=ADFL,Number=1,Type=String,Description="Number of flanking reads consistent with the allele">
    ##FORMAT=<ID=ADIR,Number=1,Type=String,Description="Number of in-repeat reads consistent with the allele">
    ##FORMAT=<ID=ADSP,Number=1,Type=String,Description="Number of spanning reads consistent with the allele">

    TRGT
    ##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">
    ##FORMAT=<ID=ALLR,Number=.,Type=String,Description="Length range per allele">
    ##FORMAT=<ID=SD,Number=.,Type=Integer,Description="Number of spanning reads supporting per allele">
    ##FORMAT=<ID=MC,Number=.,Type=String,Description="Motif counts per allele">
    ##FORMAT=<ID=MS,Number=.,Type=String,Description="Motif spans per allele">
    ##FORMAT=<ID=AP,Number=.,Type=Float,Description="Allele purity per allele">
    ##FORMAT=<ID=AM,Number=.,Type=Float,Description="Mean methylation level per allele">
    ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">

    STRDROP
    ##FORMAT=<ID=SDP,Number=1,Type=Float,Description="Strdrop coverage sequencing depth level probability">
    ##FORMAT=<ID=EDR,Number=1,Type=Float,Description="Strdrop allele similarity Levenshtein edit distance ratio">
    ##FORMAT=<ID=SDR,Number=1,Type=Float,Description="Strdrop case average adjusted sequencing depth ratio">
    ##FORMAT=<ID=DROP,Number=1,Type=String,Description="Strdrop coverage drop detected, 1 for LowDepth">

    MEI specific format fields
    ##FORMAT=<ID=CLIP3,Number=1,Type=Float,Description="Number of soft clipped reads downstream of the breakpoint">
    ##FORMAT=<ID=CLIP5,Number=1,Type=Float,Description="Number of soft clipped reads upstream of the breakpoint">
    ##FORMAT=<ID=SP,Number=1,Type=Float,Description="Number of correctly mapped read pairs spanning breakpoint, useful for estimation of size of insertion">
    ##FORMAT=<ID=SP,Number=1,Type=Float,Description="Number of correctly mapped read pairs spanning breakpoint, useful for estimation of size of insertion">
*/
/// Extracts the GT field for the sample identified by its VCF index
/// and returns a MongoDB document containing the genotype call.
fn parse_genotype(
    record: &Record,
    sample_id: &str,
    display_name: &str,
    pos: usize,
    category: VariantCategory,
) -> Document {
    let mut gt_call = doc! {
        "sample_id": sample_id,
        "display_name": display_name,
    };

    if record.format(b"GT").integer().is_ok() {
        let genotypes = record.genotypes().expect("Could not read genotypes");
        let genotype = genotypes.get(pos);

        let allele_1 = genotype_allele_to_string(genotype.get(0));
        let allele_2 = genotype_allele_to_string(genotype.get(1));

        let phase_sep = match genotype.get(1) {
            Some(GenotypeAllele::Phased(_)) 
            | Some(GenotypeAllele::PhasedMissing) => "|",
            _ => "/",
        };

        gt_call.insert(
            "genotype_call",
            format!("{}{}{}", allele_1, phase_sep, allele_2),
        );
    }

    // STR-specific fields
    let mut spanning_ref = None;
    let mut spanning_alt = None;
    let mut flanking_ref = None;
    let mut flanking_alt = None;
    let mut inrepeat_ref = None;
    let mut inrepeat_alt = None;
    let mut sd_ref = None;
    let mut sd_alt = None;

    if category == VariantCategory::Str {
        if let Some(so) = get_str_so(record, pos) {
            gt_call.insert("so", so);
        }

        (spanning_ref, spanning_alt) =
            parse_format_entry(record, pos, b"ADSP");

        (flanking_ref, flanking_alt) =
            parse_format_entry(record, pos, b"ADFL");

        (inrepeat_ref, inrepeat_alt) =
            parse_format_entry(record, pos, b"ADIR");

        // TRGT long read STR specific
        let (_, mc_alt) = parse_format_entry_trgt_mc(record, pos);
        gt_call.insert("alt_mc", mc_alt);

        (sd_ref, sd_alt) =
            parse_format_entry(record, pos, b"SD");

        // STRdrop long read STR specific
        if let Some(sdp) = parse_format_entry_single_float(record, pos, b"SDP") {
            gt_call.insert("sdp", sdp);
        }

        if let Some(edr) = parse_format_entry_single_float(record, pos, b"EDR") {
            gt_call.insert("edr", edr);
        }

        if let Some(sdr) = parse_format_entry_single_float(record, pos, b"SDR") {
            gt_call.insert("sdr", sdr);
        }

        if let Some(drop) = parse_format_entry_single_string(record, pos, b"DROP") {
            gt_call.insert("drop", drop);
        }
    }

    // MEI-specific fields
    let (spanning_mei_ref, clip5_alt, clip3_alt) = get_mei_reads(record, pos);

    // SV-specific fields
    let (paired_end_ref, paired_end_alt) = get_paired_ends(record, pos);
    let (split_read_ref, split_read_alt) = get_split_reads(record, pos);

    let alt_depth = get_alt_depth(
        record,
        pos,
        paired_end_alt,
        split_read_alt,
        spanning_alt,
        flanking_alt,
        inrepeat_alt,
        sd_alt,
        clip5_alt,
        clip3_alt,
    );

    gt_call.insert("alt_depth", alt_depth);

    let ref_depth = get_ref_depth(
        record,
        pos,
        paired_end_ref,
        split_read_ref,
        spanning_ref,
        flanking_ref,
        inrepeat_ref,
        sd_ref,
        spanning_mei_ref,
    );

    gt_call.insert("ref_depth", ref_depth);

    let read_depth = get_read_depth(record, pos, alt_depth, ref_depth);
    gt_call.insert("read_depth", read_depth);

    let alt_frequency = get_alt_frequency(record, pos);
    gt_call.insert("alt_frequency", alt_frequency);

    let genotype_quality = get_genotype_quality(record, pos);
    gt_call.insert("genotype_quality", genotype_quality);

    let ffpm = get_ffpm_info(record, pos);
    gt_call.insert("ffpm", ffpm);

    gt_call.insert("split_read", split_read_alt);

    let copy_number = get_copy_number(record, pos);
    gt_call.insert("copy_number", copy_number);


    let mut gt_obj = doc! {
        "sample_id": gt_call.get("sample_id"),
        "display_name": gt_call.get("display_name"),
        "genotype_call": gt_call.get("genotype_call"),
        "allele_depths": [
            gt_call.get("ref_depth"),
            gt_call.get("alt_depth"),
        ],
        "read_depth": gt_call.get("read_depth"),
        "alt_frequency": gt_call
            .get("alt_frequency")
            .cloned()
            .unwrap_or(Bson::Int32(-1)),
        "genotype_quality": gt_call.get("genotype_quality"),
    };


    for format_tag in [
        "alt_mc",
        "copy_number",
        "edr",
        "ffpm",
        "sdp",
        "sdr",
        "so",
        "split_read",
    ] {
        if let Some(value) = gt_call.get(format_tag) {
            gt_obj.insert(format_tag, value.clone());
        }
    }

    gt_obj
}

/// Converts a genotype allele to its string representation.
///
/// Maps allele indexes to their VCF genotype representation.
/// Missing or unsupported alleles are represented as ".".
fn genotype_allele_to_string(allele: Option<&GenotypeAllele>) -> &'static str {
    match allele {
        Some(GenotypeAllele::Unphased(0)) | Some(GenotypeAllele::Phased(0)) => "0",
        Some(GenotypeAllele::Unphased(1)) | Some(GenotypeAllele::Phased(1)) => "1",
        Some(GenotypeAllele::UnphasedMissing)
        | Some(GenotypeAllele::PhasedMissing)
        | None => ".",
        _ => ".",
    }
}

/// Parse a FORMAT entry containing reference and alternative values.
///
/// The FORMAT field can contain either one value (ALT only) or two values
/// (REF and ALT). Missing or negative values are returned as None.
fn parse_format_entry(
    record: &Record,
    sample_pos: usize,
    format_entry_name: &[u8],
) -> (Option<i32>, Option<i32>) {
    let mut ref_value = None;
    let mut alt_value = None;

    let values = match record.format(format_entry_name).integer() {
        Ok(values) => values,
        Err(_) => return (None, None),
    };

    let sample_values = values.get(sample_pos);

    if let Some(values) = sample_values {
        if values.len() > 1 {
            if values[0] >= 0 {
                ref_value = Some(values[0]);
            }
            if values[1] >= 0 {
                alt_value = Some(values[1]);
            }
        } else if values.len() == 1 && values[0] >= 0 {
            alt_value = Some(values[0]);
        }
    }

    (ref_value, alt_value)
}

/// Get the STR Sequence Ontology (SO) annotation for a sample.
///
/// Returns None if the SO FORMAT field is missing or unavailable.
fn get_str_so(record: &Record, pos: usize) -> Option<String> {
    let so_values = match record.format(b"SO").string() {
        Ok(values) => values,
        Err(_) => return None,
    };

    let so = so_values.get(pos)?;

    if so.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(so).to_string())
    }
}

/// Get the PathologicStruc INFO annotation from a VCF record.
///
/// Returns the motif indexes contributing to pathological size.
/// Returns None when the INFO field is missing or cannot be parsed.
fn get_pathologic_struc(record: &Record) -> Option<Vec<usize>> {
    let values = record.info(b"PathologicStruc").string().ok()??;

    let value = values.first()?;

    let value = String::from_utf8_lossy(value);

    Some(
        value
            .trim_matches(|c| c == '[' || c == ']')
            .split(',')
            .filter_map(|x| x.trim().parse::<usize>().ok())
            .collect(),
    )
}

/// Parse TRGT MC FORMAT entry.
///
/// The MC FORMAT field contains motif counts for each allele. Motif counts
/// are separated by "_" and alleles are separated by ",". If a
/// PathologicStruc INFO field is present, only the specified motif indexes
/// contribute to the pathological size.
///
/// # Arguments
///
/// * `record` - The VCF record containing MC FORMAT and PathologicStruc INFO.
/// * `pos` - The sample index in the VCF record.
///
/// # Returns
///
/// A tuple containing:
/// * reference motif count
/// * alternative motif count
fn parse_format_entry_trgt_mc(
    record: &Record,
    pos: usize,
) -> (Option<i32>, Option<i32>) {
    let mut mc_ref = None;
    let mut mc_alt = None;

    let mc_values = match record.format(b"MC").string() {
        Ok(values) => values,
        Err(_) => return (mc_ref, mc_alt),
    };

    let mc = match mc_values.get(pos) {
        Some(value) if !value.is_empty() => String::from_utf8_lossy(value),
        _ => return (mc_ref, mc_alt),
    };

    let genotypes = match record.genotypes() {
        Ok(genotypes) => genotypes,
        Err(_) => return (mc_ref, mc_alt),
    };

    let genotype = genotypes.get(pos);

    // Find which genotype position corresponds to the reference allele (0)
    let ref_idx = genotype
        .iter()
        .position(|allele| allele.index() == Some(0));

    let pathologic_struc = get_pathologic_struc(record);

    for (idx, allele) in mc.split(',').enumerate() {
        let motifs: Vec<&str> = allele.split('_').collect();

        let count = if motifs.len() > 1 {
            motifs
                .iter()
                .enumerate()
                .filter(|(motif_idx, _)| {
                    pathologic_struc
                        .as_ref()
                        .map(|indexes| indexes.contains(motif_idx))
                        .unwrap_or(true)
                })
                .filter_map(|(_, value)| value.parse::<i32>().ok())
                .sum()
        } else if allele == "." {
            0
        } else {
            allele.parse::<i32>().unwrap_or(0)
        };

        if Some(idx) == ref_idx {
            mc_ref = Some(count);
        } else {
            mc_alt = Some(count);
        }
    }

    (mc_ref, mc_alt)
}

/// Parse a single string FORMAT entry.
///
/// Returns the value for the selected sample, or None if the FORMAT field
/// is missing or empty.
fn parse_format_entry_single_string(
    record: &Record,
    pos: usize,
    format_entry_name: &[u8],
) -> Option<String> {
    let values = record.format(format_entry_name).string().ok()?;

    let value = values.get(pos)?;

    if value.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(value).to_string())
    }
}

/// Parse a single floating-point FORMAT entry.
///
/// Returns the value for the selected sample, or None if the FORMAT field
/// is missing or cannot be parsed.
fn parse_format_entry_single_float(
    record: &Record,
    pos: usize,
    format_entry_name: &[u8],
) -> Option<f64> {
    let values = record.format(format_entry_name).float().ok()?;

    let value = values.get(pos)?.first()?;

    if *value >= 0.0 {
        Some(*value as f64)
    } else {
        None
    }
}

/// Get MEI caller read details from FORMAT fields.
///
/// Returns:
/// * number of reference spanning reads (`SP`)
/// * number of alternative 5' clipped reads (`CLIP5`)
/// * number of alternative 3' clipped reads (`CLIP3`)
///
/// Missing fields or invalid/negative values are returned as `None`.
fn get_mei_reads(
    record: &Record,
    pos: usize,
) -> (Option<i32>, Option<i32>, Option<i32>) {
    let spanning_ref = parse_format_entry_single_integer(record, pos, b"SP");
    let clip5_alt = parse_format_entry_single_integer(record, pos, b"CLIP5");
    let clip3_alt = parse_format_entry_single_integer(record, pos, b"CLIP3");

    (spanning_ref, clip5_alt, clip3_alt)
}

/// Parse a single integer FORMAT entry.
///
/// Returns the value for the selected sample, or None if the field is
/// missing, invalid, or contains a negative value.
fn parse_format_entry_single_integer(
    record: &Record,
    pos: usize,
    format_entry_name: &[u8],
) -> Option<i32> {
    let values = record.format(format_entry_name).integer().ok()?;

    let value = values.get(pos)?.first()?;

    if *value >= 0 {
        Some(*value)
    } else {
        None
    }
}


/// Get paired-end read support from SV FORMAT fields.
///
/// Returns:
/// * reference paired-end support
/// * alternative paired-end support
///
/// Values are extracted from PE, PR, DV, and DR FORMAT fields.
fn get_paired_ends(
    record: &Record,
    pos: usize,
) -> (Option<i32>, Option<i32>) {
    let mut paired_end_ref = None;
    let mut paired_end_alt = None;

    // PE: Number of paired-end reads supporting the variant
    if let Ok(values) = record.format(b"PE").integer() {
        if let Some(value) = values.get(pos).and_then(|v| v.first()) {
            if *value >= 0 {
                paired_end_alt = Some(*value);
            }
        }
    }

    // PR: Number of paired-end reads supporting ref and alt alleles
    if let Ok(values) = record.format(b"PR").integer() {
        if let Some(sample_values) = values.get(pos) {
            if let Some(ref_value) = sample_values.first() {
                if *ref_value >= 0 {
                    paired_end_ref = Some(*ref_value);
                }
            }

            if let Some(alt_value) = sample_values.get(1) {
                if *alt_value >= 0 {
                    paired_end_alt = Some(*alt_value);
                }
            }
        }
    }

    // DV: Number of paired-end reads supporting the event
    if let Ok(values) = record.format(b"DV").integer() {
        if let Some(value) = values.get(pos).and_then(|v| v.first()) {
            if *value >= 0 {
                paired_end_alt = Some(*value);
            }
        }
    }

    // DR: Number of paired-end reads supporting the reference
    if let Ok(values) = record.format(b"DR").integer() {
        if let Some(value) = values.get(pos).and_then(|v| v.first()) {
            if *value >= 0 {
                paired_end_ref = Some(*value);
            }
        }
    }

    (paired_end_ref, paired_end_alt)
}

/// Get split-read support from SV FORMAT fields.
///
/// Returns:
/// * reference split-read support
/// * alternative split-read support
///
/// Values are extracted from SR, RV, and RR FORMAT fields.
fn get_split_reads(
    record: &Record,
    pos: usize,
) -> (Option<i32>, Option<i32>) {
    let mut split_read_ref = None;
    let mut split_read_alt = None;

    // SR: Number of split reads supporting ref and alt alleles
    if let Ok(values) = record.format(b"SR").integer() {
        if let Some(sample_values) = values.get(pos) {
            let (mut alt_value, mut ref_value) = (None, None);

            if sample_values.len() == 1 {
                alt_value = sample_values.first().copied();
            }

            if sample_values.len() == 2 {
                ref_value = sample_values.first().copied();
                alt_value = sample_values.get(1).copied();
            }

            if let Some(value) = alt_value {
                if value >= 0 {
                    split_read_alt = Some(value);
                }
            }

            if let Some(value) = ref_value {
                if value >= 0 {
                    split_read_ref = Some(value);
                }
            }
        }
    }

    // RV: Number of split reads supporting the event
    if let Ok(values) = record.format(b"RV").integer() {
        if let Some(value) = values.get(pos).and_then(|v| v.first()) {
            if *value >= 0 {
                split_read_alt = Some(*value);
            }
        }
    }

    // RR: Number of split reads supporting the reference
    if let Ok(values) = record.format(b"RR").integer() {
        if let Some(value) = values.get(pos).and_then(|v| v.first()) {
            if *value >= 0 {
                split_read_ref = Some(*value);
            }
        }
    }

    (split_read_ref, split_read_alt)
}

/// Extract alternative allele depth from the AD FORMAT field.
fn get_gt_allele_depth(
    record: &Record,
    pos: usize,
    allele_index: usize,
) -> Option<i32> {
    let depths = record.format(b"AD").integer().ok()?;

    depths
        .get(pos)
        .and_then(|sample| sample.get(allele_index).copied())
        .filter(|value| *value >= 0)
}

/// Get alternative read depth.
///
/// First tries to use the genotype-derived alternative depth.
/// If unavailable, falls back to caller-specific FORMAT fields.
fn get_alt_depth(
    record: &Record,
    pos: usize,
    paired_end_alt: Option<i32>,
    split_read_alt: Option<i32>,
    spanning_alt: Option<i32>,
    flanking_alt: Option<i32>,
    inrepeat_alt: Option<i32>,
    sd_alt: Option<i32>,
    clip5_alt: Option<i32>,
    clip3_alt: Option<i32>,
) -> i32 {
    let mut alt_depth = get_gt_allele_depth(record, pos, 1).unwrap_or(-1);

    if alt_depth != -1 {
        return alt_depth;
    }

    // VD: Number of variant supporting reads
    if let Ok(values) = record.format(b"VD").integer() {
        if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
            alt_depth = *value;
        }
    }

    let alt_items: &[&[Option<i32>]] = &[
        &[sd_alt],
        &[paired_end_alt, split_read_alt],
        &[clip5_alt, clip3_alt],
        &[spanning_alt, flanking_alt, inrepeat_alt],
    ];

    for items in alt_items {
        if items.iter().any(|item| item.is_some()) {
            alt_depth = items
                .iter()
                .filter_map(|item| *item)
                .filter(|value| *value != 0)
                .sum();
        }
    }

    alt_depth
}

/// Get reference read depth.
///
/// First tries to use the genotype-derived reference depth.
/// If unavailable, falls back to caller-specific FORMAT fields.
fn get_ref_depth(
    record: &Record,
    pos: usize,
    paired_end_ref: Option<i32>,
    split_read_ref: Option<i32>,
    spanning_ref: Option<i32>,
    flanking_ref: Option<i32>,
    inrepeat_ref: Option<i32>,
    sd_ref: Option<i32>,
    spanning_mei_ref: Option<i32>,
) -> i32 {
    let mut ref_depth = get_gt_allele_depth(record, pos, 0).unwrap_or(-1);

    if ref_depth != -1 {
        return ref_depth;
    }

    let ref_items: &[&[Option<i32>]] = &[
        &[sd_ref],
        &[paired_end_ref, split_read_ref],
        &[spanning_ref, flanking_ref, inrepeat_ref],
    ];

    for items in ref_items {
        if items.iter().any(|item| item.is_some()) {
            ref_depth = items
                .iter()
                .filter_map(|item| *item)
                .filter(|value| *value != 0)
                .sum();
        }
    }

    // MEI SP can add additional reference spanning support
    if let Some(value) = spanning_mei_ref {
        ref_depth += value;
    }

    ref_depth
}


/// Get total read depth.
///
/// First tries to use the genotype-derived read depth.
/// If unavailable, falls back to DP, LC, or the sum of the
/// reference and alternative depths.
fn get_read_depth(
    record: &Record,
    pos: usize,
    alt_depth: i32,
    ref_depth: i32,
) -> i32 {
    let mut read_depth = -1;

    if let Ok(values) = record.format(b"DP").integer() {
        if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
            read_depth = *value;
        }
    }

    if read_depth == -1 {
        // DP: Total read depth
        if let Ok(values) = record.format(b"DP").integer() {
            if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
                read_depth = *value;
            }
        }
        // LC: Locus coverage
        else if let Ok(values) = record.format(b"LC").float() {
            if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
                read_depth = value.round() as i32;
            }
        }
        // Fall back to the sum of alt and ref depths
        else if alt_depth != -1 || ref_depth != -1 {
            read_depth = 0;

            if alt_depth != -1 {
                read_depth += alt_depth;
            }

            if ref_depth != -1 {
                read_depth += ref_depth;
            }
        }
    }

    read_depth
}

/// Get alternative allele frequency.
///
/// Prioritises the caller-provided AF FORMAT field if available.
/// Otherwise calculates it from genotype allele depths.
fn get_alt_frequency(
    record: &Record,
    pos: usize,
) -> f32 {
    // AF FORMAT field (caller-provided)
    if let Ok(values) = record.format(b"AF").float() {
        if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
            return *value;
        }
    }

    // Fallback: calculate from genotype allele depths
    let alt_depth = get_gt_allele_depth(record, pos, 1).unwrap_or(-1);
    let ref_depth = get_gt_allele_depth(record, pos, 0).unwrap_or(-1);

    if alt_depth == -1 || ref_depth == -1 {
        return -1.0;
    }

    let total_depth = alt_depth + ref_depth;

    if total_depth == 0 {
        return 0.0;
    }

    alt_depth as f32 / total_depth as f32
}

/// Get FFPM information from FORMAT tags.
///
/// Returns the fusion fragments per million value if available.
fn get_ffpm_info(
    record: &Record,
    pos: usize,
) -> Option<i32> {
    if let Ok(values) = record.format(b"FFPM").integer() {
        if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
            return Some(*value);
        }
    }

    None
}

/// Get genotype quality (GQ) for a sample.
///
/// Returns -1 if genotype quality is missing or unavailable.
fn get_genotype_quality(
    record: &Record,
    pos: usize,
) -> i32 {
    if let Ok(values) = record.format(b"GQ").integer() {
        if let Some(value) = values.get(pos).and_then(|sample| sample.first()) {
            return *value;
        }
    }

    -1
}

/// Get copy number from the CN FORMAT field.
///
/// Returns None if CN is missing, invalid, or represents a missing value.
fn get_copy_number(
    record: &Record,
    pos: usize,
) -> Option<i32> {
    let values = record.format(b"CN").float().ok()?;

    let cn_value = values
        .get(pos)
        .and_then(|sample| sample.first())?;

    if cn_value.is_nan() || *cn_value < 0.0 {
        return None;
    }

    Some(*cn_value as i32)
}