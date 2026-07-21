use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::Record;
use std::collections::HashMap;
use crate::models::variant::Coordinates;
use crate::models::cytoband::Cytoband;


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


/// Extracts variant coordinates and basic variant information from a VCF record.
///
/// The chromosome name is normalized by removing a possible `chr` prefix.
/// Positions are converted from the VCF 0-based representation to 1-based
/// coordinates.
///
/// # Arguments
///
/// * `record` - A VCF record containing variant information.
/// * `header` - The VCF header used to resolve chromosome identifiers.
/// * `cytobands` - The VCF header used to resolve chromosome identifiers.
///
/// # Returns
///
/// A [`Coordinates`] object containing:
///
/// * chromosome name
/// * position
/// * end coordinate
/// * end chromosome
/// * variant length
/// * variant sub-category
/// * mate ID when present
/// * cytoband start
/// * cytoband end
pub fn parse_coordinates(
    record: &Record,
    header: &HeaderView,
    cytobands: &HashMap<String, Vec<Cytoband>>
) -> Coordinates {
    let rid = record.rid().expect("missing chromosome");

    let chrom = String::from_utf8_lossy(
        header.rid2name(rid).expect("unknown chromosome")
    )
    .to_string();

    let chrom = normalize_chromosome(&chrom);

    let position = (record.pos() + 1) as u64;

    let reference = String::from_utf8_lossy(record.alleles()[0])
        .to_string();

    let alternative = record
        .alleles()
        .get(1)
        .map(|allele| String::from_utf8_lossy(allele).to_string())
        .unwrap_or_else(|| ".".to_string());

    let ref_len = reference.len();
    let alt_len = alternative.len();

    let mut sub_category = "snv".to_string();

    let end = record.end() as u64;

    let length = if ref_len != alt_len {
        sub_category = "indel".to_string();
        (ref_len as i64 - alt_len as i64).abs()
    } else {
        alt_len as i64
    };

    let mut end_chrom = chrom.clone();

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
        chromosome: chrom.clone(),
        position,
        end,
        end_chrom,
        length,
        sub_category,
        mate_id: mate_id,
        cytoband_start,
        cytoband_end
    }
}