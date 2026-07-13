use rust_htslib::bcf::{Read, Record};
use rust_htslib::bcf::header::HeaderView;

/// Parses the FILTER field from a VCF record.
///
/// Converts filter identifiers stored in the VCF record into their corresponding
/// filter names using the VCF header. If the record has no filters, returns
/// "PASS" as the default filter value.
///
/// # Arguments
///
/// * `record` - The VCF record containing the filter information.
/// * `header` - The VCF header used to resolve filter identifiers into names.
///
/// # Returns
///
/// A vector of filter names associated with the variant.
pub fn parse_filters(
    record: &Record,
    header: &HeaderView,
) -> Vec<String> {
    let filters: Vec<String> = record
        .filters()
        .map(|id| {
            let name = header.id_to_name(id);
            String::from_utf8_lossy(&name).to_string()
        })
        .collect();

    if filters.is_empty() {
        vec!["PASS".to_string()]
    } else {
        filters
    }
}