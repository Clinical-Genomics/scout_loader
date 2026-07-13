use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::Record;


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


/// Extracts chromosome, position, and end coordinates from a VCF record.
///
/// # Arguments
///
/// * `record` - A VCF record containing variant information.
/// * `header` - The VCF header used to resolve chromosome identifiers.
///
/// # Returns
///
/// A tuple containing:
/// * chromosome name
/// * 1-based variant position
/// * variant end coordinate
pub fn parse_coordinates(
    record: &Record,
    header: &HeaderView,
) -> (String, u64, u64) {
    let rid = record.rid().expect("missing chromosome");

    let chromosome = String::from_utf8_lossy(
        header.rid2name(rid).expect("unknown chromosome")
    )
    .to_string();
    let chromosome = normalize_chromosome(&chromosome);

    let position = (record.pos() + 1) as u64;

    let end = position;

    (chromosome, position, end)
}