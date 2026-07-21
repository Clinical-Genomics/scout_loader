use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::Record;
use std::collections::HashMap;
use crate::models::cytoband::Cytoband;
use crate::models::variant::{Coordinates, VariantCategory};
use crate::parse::info::parse_info_int;
use lazy_static::lazy_static;
use regex::Regex;

lazy_static! {
    static ref BND_ALT_PATTERN: Regex =
        Regex::new(r".*[\],\[](.*?):(.*?)[\],\[]").unwrap();
}

const SV_TYPES: &[&str] = &["ins", "del", "dup", "cnv", "inv", "bnd"];


/// Finds the cytoband overlapping a genomic coordinate.
///
/// # Arguments
///
/// * `cytobands` - Cytoband annotations indexed by chromosome.
/// * `chrom` - Normalized chromosome name.
/// * `pos` - Genomic position (1-based).
///
/// # Returns
///
/// The cytoband name if the position overlaps an interval, otherwise an empty
/// string.
pub fn get_cytoband_coordinates(
    cytobands: &HashMap<String, Vec<Cytoband>>,
    chrom: &str,
    pos: u64,
) -> Option<String> {
    cytobands
        .get(chrom)
        .and_then(|bands| {
            bands.iter().find(|band| {
                pos >= band.start && pos < band.end
            })
        })
        .map(|band| band.name.clone())
}


/// Normalizes a chromosome name by removing an optional "chr" prefix.
///
/// The prefix is removed in a case-insensitive manner (e.g. "chr1",
/// "CHR1", and "Chr1" are all normalized to "1").
fn normalize_chromosome(chromosome: &str) -> String {
    if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") {
        chromosome[3..].to_string()
    } else {
        chromosome.to_string()
    }
}

/// Return the end chromosome for a translocation.
///
/// BND variants can represent translocations between different chromosomes.
/// The chromosome of the mate breakpoint is extracted from the ALT allele.
/// If no chromosome can be extracted, the original chromosome is returned.
fn get_end_chrom(alt: &str, chrom: &str) -> String {
    if !alt.contains(':') {
        return chrom.to_string();
    }

    if let Some(bnd_match) = BND_ALT_PATTERN.captures(alt) {
        if let Some(other_chrom) = bnd_match.get(1) {
            return normalize_chromosome(other_chrom.as_str());
        }
    }

    chrom.to_string()
}

/// Find SV type for structural variants.
///
/// The SVTYPE INFO tag is deprecated in the VCF standard but may still
/// be present. If available, it is preferred. Otherwise, the type is
/// inferred from the ALT allele.
fn get_svtype(record: &Record, alt: &str, alt_len: usize) -> String {
    if let Ok(Some(values)) = record.info(b"SVTYPE").string() {
        if let Some(value) = values.first() {
            let mut svtype = String::from_utf8_lossy(value).to_lowercase();

            if svtype == "sgl" {
                svtype = "bnd".to_string();
            }

            return svtype;
        }
    }

    let alt_type = alt
        .trim_start_matches('<')
        .trim_end_matches('>')
        .to_lowercase();

    if SV_TYPES.contains(&alt_type.as_str()) {
        return alt_type;
    }

    if alt.contains('.') && alt_len > 1 {
        return "bnd".to_string();
    }

    unreachable!("Unable to determine SV type")
}

/// Return the end coordinate for a structural variant.
///
/// The END INFO field is usually sufficient, but some callers set END
/// equal to POS for variants such as insertions. In those cases SVLEN
/// can be used as a fallback.
///
/// Breakends (BNDs) require special handling because the end coordinate
/// can be encoded in the ALT allele pattern.
fn sv_end(
    pos: u64,
    alt: &str,
    svend: Option<i64>,
    svlen: Option<i64>,
) -> u64 {
    let mut end = svend.map(|value| value as u64);

    if alt.contains(':') {
        if let Some(captures) = BND_ALT_PATTERN.captures(alt) {
            if let Some(position) = captures.get(2) {
                end = position.as_str().parse::<u64>().ok();
            }
        }
    } else if alt.contains('.') && alt.len() > 1 {
        end = Some(pos);
    }

    if end.is_none() {
        if let Some(svlen) = svlen {
            end = Some((pos as i64 + svlen) as u64);
        }
    }

    end.unwrap_or(pos)
}

/// Return the length of a structural variant.
///
/// Returns a very large value for variants spanning different molecules.
/// Uses SVLEN when available. If no length can be determined, returns -1.
fn sv_length(
    pos: u64,
    end: u64,
    chrom: &str,
    end_chrom: &str,
    svlen: Option<i64>,
) -> i64 {
    if chrom != end_chrom {
        return 100_000_000_000;
    }

    if let Some(svlen) = svlen {
        return svlen.abs();
    }

    if end == 0 || end == pos {
        return -1;
    }

    end as i64 - pos as i64
}


/// Parse genomic coordinates and variant classification information.
///
/// Extracts chromosome, position, reference and alternative alleles from a VCF
/// record and computes variant-specific coordinates:
/// - SNVs and indels use the reference/alternative allele lengths.
/// - SVs, fusions, and cancer SVs use SV-specific END, SVLEN, and BND handling.
/// - MEIs use insertion length information.
///
/// Also extracts optional mate IDs and cytoband annotations.
///
/// # Arguments
///
/// * `record` - VCF record containing variant information.
/// * `header` - VCF header used to resolve chromosome names.
/// * `cytobands` - Cytoband definitions used to annotate start and end positions.
/// * `variant_type` - Variant category used to determine coordinate calculation.
///
/// # Returns
///
/// A `Coordinates` object containing chromosome, position, end, length,
/// sub-category, mate information, and cytoband annotations.
pub fn parse_coordinates(
    record: &Record,
    header: &HeaderView,
    cytobands: &HashMap<String, Vec<Cytoband>>,
    category: &VariantCategory,
) -> Coordinates {
    let rid = record.rid().expect("missing chromosome");

    let chrom = String::from_utf8_lossy(
        header.rid2name(rid).expect("unknown chromosome"),
    )
    .to_string();

    let chrom = normalize_chromosome(&chrom);

    let position = (record.pos() + 1) as u64;

    let reference = String::from_utf8_lossy(record.alleles()[0]).to_string();

    let alternative = record
        .alleles()
        .get(1)
        .map(|allele| String::from_utf8_lossy(allele).to_string())
        .unwrap_or_else(|| ".".to_string());

    let ref_len = reference.len();
    let alt_len = alternative.len();

    let mut sub_category = "snv".to_string();
    let mut end = record.end() as u64;
    let mut end_chrom = chrom.clone();
    let mut length = alt_len as i64;

    match category {
        VariantCategory::Sv | VariantCategory::CancerSv | VariantCategory::Fusion => {
            let svtype = get_svtype(record, &alternative, alt_len);

            sub_category = svtype.clone();

            if sub_category == "bnd" {
                end_chrom = get_end_chrom(&alternative, &chrom);
            }

            end = sv_end(
                position,
                &alternative,
                parse_info_int(record, b"END").map(|v| v as i64),
                parse_info_int(record, b"SVLEN").map(|v| v as i64),
            );

            length = sv_length(
                position,
                end,
                &chrom,
                &end_chrom,
                record
                    .info(b"SVLEN")
                    .integer()
                    .ok()
                    .flatten()
                    .and_then(|values| values.first().copied())
                    .map(|value| value as i64),
                    );
        }

        VariantCategory::Mei => {
            sub_category = "mei".to_string();

            length = alt_len as i64;

            if ref_len != alt_len {
                length = (ref_len as i64 - alt_len as i64).abs();
            }
        }

        _ => {
            if ref_len != alt_len {
                sub_category = "indel".to_string();
                length = (ref_len as i64 - alt_len as i64).abs();
            }
        }
    }

    let mate_id = record
        .info(b"MATEID")
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values
                .first()
                .map(|value| String::from_utf8_lossy(value).to_string())
        });

    let cytoband_start = get_cytoband_coordinates(
        cytobands,
        &chrom,
        position,
    );

    let cytoband_end = get_cytoband_coordinates(
        cytobands,
        &end_chrom,
        end,
    );

    Coordinates {
        chromosome: chrom,
        position,
        end,
        end_chrom,
        length,
        sub_category,
        mate_id,
        cytoband_start,
        cytoband_end,
    }
}