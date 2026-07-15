use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::Record;
use crate::models::variant::Coordinates;


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
///
/// # Returns
///
/// A [`Coordinates`] object containing:
///
/// * chromosome name
/// * position
/// * end coordinate
/// * reference and alternative alleles
/// * variant length
/// * variant sub-category
/// * mate ID when present
pub fn parse_coordinates(
    record: &Record,
    header: &HeaderView,
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

    Coordinates {
        chromosome: chrom.clone(),
        position,
        end,
        end_chrom: chrom,
        length,
        sub_category,
        mate_id: None,
        cytoband_start: None,
        cytoband_end: None,
    }
}