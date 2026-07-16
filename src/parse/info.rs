use rust_htslib::bcf::Record;
use std::str::FromStr;

/// Parses a VCF INFO field into the requested type.
///
/// Returns `None` if the INFO field is missing or if the value cannot be
/// converted to the requested type.
///
/// # Arguments
///
/// * `record` - VCF record containing the INFO annotation.
/// * `tag` - INFO field name.
///
/// # Type Parameters
///
/// * `T` - Target type implementing `FromStr`.
///
/// # Returns
///
/// The parsed value wrapped in `Some`, or `None` if parsing fails.
pub fn parse_info<T>(
    record: &Record,
    tag: &[u8],
) -> Option<T>
where
    T: FromStr,
{
    record
        .info(tag)
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values.first().map(|value| {
                String::from_utf8_lossy(value).to_string()
            })
        })
        .and_then(|value| value.parse::<T>().ok())
}