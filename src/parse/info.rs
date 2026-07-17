use mongodb::bson::{self, Document};
use rust_htslib::bcf::Record;
use mongodb::bson::Bson;

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

/// Parse a multi-valued string INFO field.
///
/// Returns all string values associated with the INFO tag,
/// or `None` if the field is missing or cannot be parsed.
///
/// This is intended for INFO fields with `Number > 1`, such as
/// `MEINFO` (`NAME,START,END,POLARITY`).
pub fn parse_info_string_array(record: &Record, key: &[u8]) -> Option<Vec<String>> {
    let values = record.info(key).string().ok().flatten()?;

    Some(
        values
            .iter()
            .map(|v| String::from_utf8_lossy(v).into_owned())
            .collect::<Vec<String>>(),
    )
}


/// Parse the SCOUT_CUSTOM INFO field.
///
/// Input format:
/// "key1|val1,key2|val2"
///
/// Returns:
/// [["key1", "val1"], ["key2", "val2"]]
///
/// Missing or malformed values return None.
pub fn parse_custom_data(custom_str: Option<String>) -> Option<Bson> {
    let custom_str = custom_str?;

    let pairs: Vec<Bson> = custom_str
        .split(',')
        .map(|pair| {
            let values: Vec<Bson> = pair
                .split('|')
                .map(|value| Bson::String(value.to_string()))
                .collect();

            Bson::Array(values)
        })
        .collect();

    Some(Bson::Array(pairs))
}