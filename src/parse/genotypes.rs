use std::collections::HashMap;
use rust_htslib::bcf::header::HeaderView;
use mongodb::bson::{doc, Document};
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
///
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

    // SV-specific fields

    // MEI-specific fields
    let (spanning_mei_ref, clip5_alt, clip3_alt) = get_mei_reads(record, pos);

    // Derived fields

    gt_call
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