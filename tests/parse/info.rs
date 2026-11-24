use crate::vcf_record::record;
use mongodb::bson::{Bson, Document, doc};
use rust_htslib::bcf::Record;
use scout_loader::parse::info::{
    insert_info_float, insert_info_int, insert_info_string, parse_custom_data, parse_info_float,
    parse_info_int, parse_info_string, parse_info_string_array, percent_decode,
};

fn info_record(info: &str) -> Record {
    record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=FLOAT,Number=.,Type=Float,Description=\"Float\">",
            "##INFO=<ID=INT,Number=.,Type=Integer,Description=\"Integer\">",
            "##INFO=<ID=STR,Number=1,Type=String,Description=\"String\">",
            "##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"NAME,START,END,POLARITY\">",
        ],
        &format!("1\t100\t.\tA\tG\t.\tPASS\t{info}"),
    )
}

#[test]
fn test_parse_info_values() {
    // GIVEN a record with multi-valued Float and Integer fields and a String
    let record = info_record("FLOAT=0.5,0.25;INT=3,4;STR=abc");

    // THEN the first value of each field is returned
    assert_eq!(parse_info_float(&record, b"FLOAT"), Some(0.5));
    assert_eq!(parse_info_int(&record, b"INT"), Some(3));
    assert_eq!(parse_info_string(&record, b"STR").as_deref(), Some("abc"));
}

#[test]
fn test_parse_info_missing_or_wrong_type() {
    // GIVEN a record without the fields
    let record = info_record(".");

    // THEN missing, undeclared and wrongly typed fields give None
    assert_eq!(parse_info_float(&record, b"FLOAT"), None);
    assert_eq!(parse_info_int(&record, b"INT"), None);
    assert_eq!(parse_info_string(&record, b"STR"), None);
    assert_eq!(parse_info_string_array(&record, b"MEINFO"), None);
    assert_eq!(parse_info_string(&record, b"UNDECLARED"), None);

    let record = info_record("STR=abc");
    assert_eq!(parse_info_float(&record, b"STR"), None);
    assert_eq!(parse_info_int(&record, b"STR"), None);
}

#[test]
fn test_parse_info_string_array() {
    // GIVEN an MEI record with a 4-valued MEINFO
    let record = info_record("MEINFO=ALU,-23,262,+");

    // THEN all values are returned
    assert_eq!(
        parse_info_string_array(&record, b"MEINFO"),
        Some(vec![
            "ALU".to_string(),
            "-23".to_string(),
            "262".to_string(),
            "+".to_string()
        ])
    );
}

#[test]
fn test_insert_info_values() {
    // GIVEN a record with all field types
    let record = info_record("FLOAT=0.5;INT=3;STR=abc");
    let mut document = Document::new();

    // WHEN inserting present and missing fields
    insert_info_float(&record, &mut document, b"FLOAT", "float");
    insert_info_int(&record, &mut document, b"INT", "int");
    insert_info_string(&record, &mut document, b"STR", "string");
    insert_info_string(&record, &mut document, b"MEINFO", "missing");

    // THEN only present fields are added, with the right types
    assert_eq!(
        document,
        doc! { "float": 0.5, "int": 3_i32, "string": "abc" }
    );
}

#[test]
fn test_parse_custom_data() {
    // GIVEN a SCOUT_CUSTOM string with two key/value pairs
    let parsed = parse_custom_data(Some("key1|val1,key2|val2".to_string()));

    // THEN it is split into pairs
    assert_eq!(
        parsed,
        Some(Bson::Array(vec![
            Bson::Array(vec!["key1".into(), "val1".into()]),
            Bson::Array(vec!["key2".into(), "val2".into()]),
        ]))
    );

    // AND a missing value gives None
    assert_eq!(parse_custom_data(None), None);
}

#[test]
fn test_percent_decode() {
    assert_eq!(
        percent_decode("GeneReviews%20Internet%202019-11-07"),
        "GeneReviews Internet 2019-11-07"
    );
    assert_eq!(
        percent_decode("Paucar%2C%20Nilsson%2C%20Engvall%20et%20al"),
        "Paucar, Nilsson, Engvall et al"
    );
    assert_eq!(percent_decode("a%3Bb%3Dc%25d"), "a;b=c%d");
    assert_eq!(percent_decode("%2c"), ",");
}

#[test]
fn test_percent_decode_keeps_plain_and_invalid_values() {
    assert_eq!(percent_decode("AD"), "AD");
    assert_eq!(
        percent_decode("Spinocerebellar ataxia"),
        "Spinocerebellar ataxia"
    );
    assert_eq!(percent_decode("100%"), "100%");
    assert_eq!(percent_decode("%2"), "%2");
    assert_eq!(percent_decode("%zz"), "%zz");
    assert_eq!(percent_decode("%C3%A9"), "é");
}
