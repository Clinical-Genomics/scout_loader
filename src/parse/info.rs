use mongodb::bson::{self, Document};
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

/// Insert a float INFO field into a MongoDB document if the field exists.
///
/// The value is read from the VCF record using `vcf_key` and inserted as a
/// BSON double under `mongo_key`. Missing INFO fields are ignored.
pub fn insert_info_float(
    record: &Record,
    doc: &mut Document,
    vcf_key: &[u8],
    mongo_key: &str,
) {
    if let Some(value) = parse_info_float(record, vcf_key) {
        doc.insert(
            mongo_key,
            bson::Bson::Double(value),
        );
    }
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

/// Insert an integer INFO field into a MongoDB document if the field exists.
///
/// The value is read from the VCF record using `vcf_key` and inserted as a
/// BSON integer under `mongo_key`. Missing INFO fields are ignored.
pub fn insert_info_int(
    record: &Record,
    doc: &mut Document,
    vcf_key: &[u8],
    mongo_key: &str,
) {
    if let Some(value) = parse_info_int(record, vcf_key) {
        doc.insert(
            mongo_key,
            bson::Bson::Int32(value),
        );
    }
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


/// Insert a string INFO field into a MongoDB document if the field exists.
///
/// The value is read from the VCF record using `vcf_key` and inserted into
/// the document under `mongo_key`. Missing INFO fields are ignored.
pub fn insert_info_string(
    record: &Record,
    doc: &mut Document,
    vcf_key: &[u8],
    mongo_key: &str,
) {
    if let Some(value) = parse_info_string(record, vcf_key) {
        doc.insert(
            mongo_key,
            bson::Bson::String(value),
        );
    }
}