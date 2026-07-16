use rust_htslib::bcf::Record;

/// Parse a Float INFO field from a VCF record.
///
/// Returns the first value if the field exists, otherwise None.
/// INFO fields can be arrays in VCF, so only the first value is used.
pub fn parse_info_float(record: &Record, key: &[u8]) -> Option<f64> {
    record
        .info(key)
        .float()
        .ok()
        .flatten()
        .and_then(|values| values.iter().next().copied())
        .map(|value| value as f64)
}

/// Parse an Integer INFO field from a VCF record.
///
/// Returns the first value if the field exists, otherwise None.
pub fn parse_info_int(record: &Record, key: &[u8]) -> Option<i32> {
    record
        .info(key)
        .integer()
        .ok()
        .flatten()
        .and_then(|values| values.iter().next().copied())
}


/// Parse a String INFO field from a VCF record.
///
/// Returns the first value if the field exists, otherwise None.
pub fn parse_info_string(record: &Record, key: &[u8]) -> Option<String> {
    match record.info(key).string() {
        Ok(Some(values)) => values
            .first()
            .map(|value| String::from_utf8_lossy(value).to_string()),
        _ => None,
    }
}