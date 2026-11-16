use crate::models::cytoband::Cytoband;
use crate::models::variant::{Coordinates, VariantCategory};
use crate::parse::info::parse_info_int;
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::bcf::Record;
use rust_htslib::bcf::header::HeaderView;
use std::collections::HashMap;

lazy_static! {
    static ref BND_ALT_PATTERN: Regex = Regex::new(r".*[\],\[](.*?):(.*?)[\],\[]").unwrap();
}

const SV_TYPES: &[&str] = &["ins", "del", "dup", "cnv", "inv", "bnd"];

/// Returns the cytoband containing the given genomic position.
///
/// The position is expected to be a 1-based genomic coordinate.
pub fn get_cytoband_coordinates(
    cytobands: &HashMap<String, Vec<Cytoband>>,
    chrom: &str,
    pos: i32,
) -> Option<String> {
    let chrom = chrom.trim_start_matches("chr");

    cytobands.get(chrom).and_then(|bands| {
        bands
            .iter()
            .find(|band| pos >= 0 && pos >= band.start && pos < band.end)
            .map(|band| band.name.clone())
    })
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

    if let Some(bnd_match) = BND_ALT_PATTERN.captures(alt)
        && let Some(other_chrom) = bnd_match.get(1)
    {
        return normalize_chromosome(other_chrom.as_str());
    }

    chrom.to_string()
}

/// Find SV type for structural variants.
///
/// The SVTYPE INFO tag is deprecated in the VCF standard but may still
/// be present. If available, it is preferred. Otherwise, the type is
/// inferred from the ALT allele.
fn get_svtype(record: &Record, alt: &str, alt_len: usize) -> String {
    if let Ok(Some(values)) = record.info(b"SVTYPE").string()
        && let Some(value) = values.first()
    {
        let mut svtype = String::from_utf8_lossy(value).to_lowercase();

        if svtype == "sgl" {
            svtype = "bnd".to_string();
        }

        return svtype;
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

/// Determines the end position of a structural variant.
///
/// The end position is resolved in the following order:
/// - For breakends, use the position encoded in the ALT allele.
/// - For symbolic alleles containing `.`, use the variant position.
/// - Otherwise, use `SVLEN` relative to the variant position when available.
/// - Fall back to the variant position if no end position can be determined.
fn sv_end(pos: i32, alt: &str, svend: Option<i32>, svlen: Option<i32>) -> i32 {
    let mut end = svend;

    if alt.contains(':') {
        if let Some(captures) = BND_ALT_PATTERN.captures(alt)
            && let Some(position) = captures.get(2)
        {
            end = position.as_str().parse::<i32>().ok();
        }
    } else if alt.contains('.') && alt.len() > 1 {
        end = Some(pos);
    }

    if end.is_none()
        && let Some(svlen) = svlen
    {
        end = Some(pos + svlen);
    }

    end.unwrap_or(pos)
}

/// Return the length of a structural variant.
///
/// Returns a very large value for variants spanning different molecules.
/// Uses SVLEN when available. If no length can be determined, returns -1.
fn sv_length(pos: i32, end: i32, chrom: &str, end_chrom: &str, svlen: Option<i32>) -> i32 {
    if chrom != end_chrom {
        return i32::MAX;
    }

    if let Some(svlen) = svlen {
        return svlen.abs();
    }

    if end == 0 || end == pos {
        return -1;
    }

    end - pos
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
/// * `category` - Variant category used to determine coordinate calculation.
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
    let chrom = match record.rid() {
        Some(rid) if header.contig_count() > 0 => {
            String::from_utf8_lossy(header.rid2name(rid).expect("unknown chromosome")).to_string()
        }
        _ => {
            let chrom = record
                .info(b"CHRA")
                .string()
                .expect("invalid CHRA")
                .expect("missing CHRA");

            String::from_utf8_lossy(chrom[0]).to_string()
        }
    };

    let chrom = normalize_chromosome(&chrom);

    let position = (record.pos() + 1) as i32;

    let reference = String::from_utf8_lossy(record.alleles()[0]).to_string();

    let alternative = record
        .alleles()
        .get(1)
        .map(|allele| String::from_utf8_lossy(allele).to_string())
        .unwrap_or_else(|| ".".to_string());

    let ref_len = reference.len();
    let alt_len = alternative.len();

    let mut sub_category = "snv".to_string();
    let mut end = record.end() as i32;
    let mut end_chrom = chrom.clone();
    let mut length = alt_len as i32;

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
                parse_info_int(record, b"END"),
                parse_info_int(record, b"SVLEN"),
            );

            length = if chrom != end_chrom {
                1_000_000_000
            } else {
                sv_length(
                    position,
                    end,
                    &chrom,
                    &end_chrom,
                    record
                        .info(b"SVLEN")
                        .integer()
                        .ok()
                        .flatten()
                        .and_then(|values| values.first().copied()),
                )
            };
        }

        VariantCategory::Mei => {
            sub_category = "mei".to_string();

            length = if ref_len != alt_len {
                (ref_len as i32 - alt_len as i32).abs()
            } else {
                alt_len as i32
            };
        }

        _ => {
            if ref_len != alt_len {
                sub_category = "indel".to_string();
                length = (ref_len as i32 - alt_len as i32).abs();
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

    let cytoband_start = get_cytoband_coordinates(cytobands, &chrom, position);

    let cytoband_end = get_cytoband_coordinates(cytobands, &end_chrom, end);

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

#[cfg(test)]
mod tests {
    use super::*;

    // Ported from scout's `tests/parse/test_parse_coordinates.py`.

    #[test]
    fn test_get_sv_length_small_ins() {
        // GIVEN an insertion with the whole sequence in ALT (POS == END)
        assert_eq!(sv_length(144343218, 144343218, "1", "1", Some(296)), 296);
    }

    #[test]
    fn test_get_sv_length_large_ins_no_length() {
        // GIVEN an imprecise insertion without SVLEN
        assert_eq!(sv_length(133920667, 133920667, "1", "1", None), -1);
    }

    #[test]
    fn test_get_sv_length_translocation() {
        // GIVEN a translocation between different chromosomes
        assert_eq!(sv_length(726044, 0, "1", "10", None), i32::MAX);
    }

    #[test]
    fn test_get_sv_length_translocation_same_chrom() {
        assert_eq!(
            sv_length(726044, 16457990, "1", "1", None),
            16457990 - 726044
        );
    }

    #[test]
    fn test_get_sv_length_cnvnator_del() {
        // GIVEN a CNVnator deletion with negative SVLEN
        assert_eq!(sv_length(1, 10000, "1", "1", Some(-10000)), 10000);
    }

    #[test]
    fn test_get_sv_length_del_no_length() {
        assert_eq!(sv_length(869314, 870246, "1", "1", None), 870246 - 869314);
    }

    #[test]
    fn test_get_sv_length_end_zero() {
        assert_eq!(sv_length(869314, 0, "1", "1", None), -1);
    }

    #[test]
    fn test_get_end_tiddit_translocation() {
        // GIVEN a TIDDIT translocation, END comes from the ALT field
        assert_eq!(sv_end(724779, "N[hs37d5:12060532[", None, None), 12060532);
    }

    #[test]
    fn test_get_end_bnd_alt_overrides_svend() {
        assert_eq!(sv_end(724779, "]1:5000]N", Some(1), None), 5000);
    }

    #[test]
    fn test_get_end_deletion() {
        assert_eq!(sv_end(869314, "<DEL>", Some(870246), None), 870246);
    }

    #[test]
    fn test_get_end_sgl_bnd() {
        // GIVEN a single-end BND, END is set to POS
        assert_eq!(sv_end(123455, "GTCA.", None, None), 123455);
    }

    #[test]
    fn test_get_end_from_svlen() {
        assert_eq!(sv_end(1000, "<DUP>", None, Some(250)), 1250);
    }

    #[test]
    fn test_get_end_fallback_to_pos() {
        assert_eq!(sv_end(1000, "<INS>", None, None), 1000);
    }

    #[test]
    fn test_get_end_chrom() {
        assert_eq!(get_end_chrom("N[hs37d5:12060532[", "1"), "hs37d5");
        assert_eq!(get_end_chrom("]CHR2:100]N", "1"), "2");
        assert_eq!(get_end_chrom("<DEL>", "1"), "1");
        assert_eq!(get_end_chrom("N:", "1"), "1");
    }

    #[test]
    fn test_normalize_chromosome() {
        assert_eq!(normalize_chromosome("chr1"), "1");
        assert_eq!(normalize_chromosome("Chr1"), "1");
        assert_eq!(normalize_chromosome("CHRX"), "X");
        assert_eq!(normalize_chromosome("MT"), "MT");
        assert_eq!(normalize_chromosome("ch"), "ch");
    }
}
