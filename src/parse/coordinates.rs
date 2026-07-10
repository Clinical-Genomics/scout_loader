use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::Record;

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

    let position = (record.pos() + 1) as u64;

    let end = position;

    (chromosome, position, end)
}