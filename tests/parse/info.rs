use crate::vcf_record::record;
use mongodb::bson::{Bson, Document, doc};
use rust_htslib::bcf::Record;
use scout_loader::parse::info::{
    insert_info_float, insert_info_int, insert_info_string, parse_custom_data, parse_info_float,
    parse_info_int, parse_info_string, parse_info_string_array, percent_decode,
};
use scout_loader::parse::{
    fusions::set_fusion_info, loqusdb_frequencies::add_loqus_archive_frequencies,
    rank_scores::parse_rank_score_other, severity::set_severity_predictions, strs::set_str_info,
};

const FLOAT_ANNOTATION_META: &[&str] = &[
    "##contig=<ID=1,length=249250621>",
    "##INFO=<ID=CADD,Number=1,Type=Float,Description=\"CADD\">",
    "##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description=\"CADD PHRED\">",
    "##INFO=<ID=SPIDEX,Number=1,Type=Float,Description=\"SPIDEX\">",
    "##INFO=<ID=clinical_genomics_loqusFrq,Number=1,Type=Float,Description=\"Loqus frequency\">",
    "##INFO=<ID=Frq,Number=1,Type=Float,Description=\"Archive frequency\">",
    "##INFO=<ID=Cancer_Germline_Frq,Number=1,Type=Float,Description=\"Cancer germline frequency\">",
    "##INFO=<ID=Cancer_Somatic_Frq,Number=1,Type=Float,Description=\"Cancer somatic frequency\">",
    "##INFO=<ID=Cancer_Somatic_Panel_Frq,Number=1,Type=Float,Description=\"Cancer panel frequency\">",
    "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"Fusion score\">",
    "##INFO=<ID=HGNC_ID_A,Number=1,Type=Float,Description=\"HGNC ID\">",
    "##INFO=<ID=TRANSCRIPT_ID_A,Number=1,Type=String,Description=\"Transcript ID\">",
    "##INFO=<ID=EXON_NUMBER_A,Number=1,Type=Float,Description=\"Exon\">",
    "##INFO=<ID=MivmirScore,Number=1,Type=Float,Description=\"Mivmir score\">",
    "##INFO=<ID=GicamScore,Number=1,Type=Float,Description=\"Gicam score\">",
    "##INFO=<ID=SweGenMean,Number=1,Type=Float,Description=\"SweGen mean\">",
    "##INFO=<ID=SweGenStd,Number=1,Type=Float,Description=\"SweGen standard deviation\">",
];

fn info_record(info: &str) -> Record {
    record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=FLOAT,Number=.,Type=Float,Description=\"Float\">",
            "##INFO=<ID=INT,Number=.,Type=Integer,Description=\"Integer\">",
            "##INFO=<ID=STR,Number=1,Type=String,Description=\"String\">",
            "##INFO=<ID=STR_DOT,Number=.,Type=String,Description=\"String\">",
            "##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description=\"Custom annotations\">",
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
fn test_parse_info_float_missing_values() {
    // GIVEN missing floats, including the first element of a multi-valued field
    for info in ["FLOAT=.", "FLOAT=.,0.5", "FLOAT=NaN"] {
        let record = info_record(info);
        let mut document = Document::new();

        // WHEN parsing and inserting the missing value
        let value = parse_info_float(&record, b"FLOAT");
        insert_info_float(&record, &mut document, b"FLOAT", "float");

        // THEN it is absent, rather than NaN or a later array element
        assert_eq!(value, None, "{info}");
        assert!(document.is_empty(), "{info}");
    }
}

#[test]
fn test_missing_float_annotations_are_omitted() {
    // GIVEN explicitly missing Float annotations and a fusion transcript
    let record = record(
        FLOAT_ANNOTATION_META,
        "1\t100\t.\tA\tG\t.\tPASS\tCADD=.;CADD_PHRED=.;SPIDEX=.;clinical_genomics_loqusFrq=.;Frq=.;Cancer_Germline_Frq=.;Cancer_Somatic_Frq=.;Cancer_Somatic_Panel_Frq=.;SCORE=.;HGNC_ID_A=1316;TRANSCRIPT_ID_A=ENST000001;EXON_NUMBER_A=.;MivmirScore=.;GicamScore=.;SweGenMean=.;SweGenStd=.",
    );
    let mut severity = Document::new();
    let mut archive = Document::new();
    let mut fusion = Document::new();
    let mut str_info = Document::new();

    // WHEN building the annotation documents
    set_severity_predictions(&mut severity, &record, &[]);
    add_loqus_archive_frequencies(&record, &mut archive, None);
    set_fusion_info(&record, &mut fusion);
    set_str_info(&record, &mut str_info);

    // THEN missing annotations are not stored as NaN or the string "NaN"
    assert_eq!(severity, doc! { "cadd_score": 0.0 });
    assert!(archive.is_empty());
    assert_eq!(parse_rank_score_other(&record), None);
    assert!(!str_info.contains_key("str_swegen_mean"));
    assert!(!str_info.contains_key("str_swegen_std"));
    assert!(!fusion.contains_key("rank_score"));
    assert!(!fusion.contains_key("fusion_score"));
    let gene = fusion.get_array("genes").unwrap()[0].as_document().unwrap();
    let transcript = gene.get_array("transcripts").unwrap()[0]
        .as_document()
        .unwrap();
    assert_eq!(
        transcript,
        &doc! {
            "transcript_id": "ENST000001",
            "hgnc_id": 1316_i32,
            "is_canonical": false,
        }
    );
}

#[test]
fn test_missing_float_annotations_preserve_fallbacks() {
    // GIVEN missing primary scores and frequencies with valid alternatives
    for (cadd_info, expected_cadd) in [
        ("CADD=.;CADD_PHRED=.", 24.0),
        ("CADD=.;CADD_PHRED=30", 30.0),
    ] {
        let record = record(
            FLOAT_ANNOTATION_META,
            &format!(
                "1\t100\t.\tA\tG\t.\tPASS\t{cadd_info};SPIDEX=-2.5;clinical_genomics_loqusFrq=.;Frq=0.25"
            ),
        );
        let mut severity = Document::new();
        let mut archive = Document::new();

        // WHEN parsing severity and archive frequencies
        set_severity_predictions(&mut severity, &record, &[doc! { "cadd": 24.0 }]);
        add_loqus_archive_frequencies(&record, &mut archive, None);

        // THEN the next INFO key or transcript is used and negative scores survive
        assert_eq!(
            severity,
            doc! { "cadd_score": expected_cadd, "spidex": -2.5 }
        );
        assert_eq!(archive, doc! { "local_obs_old_freq": 0.25 });
    }
}

#[test]
fn test_parse_info_float_preserves_zero_and_negative_values() {
    // GIVEN valid scores, which need not be positive
    for value in [0.0, -1.0, -2.5] {
        let record = info_record(&format!("FLOAT={value}"));

        // WHEN parsing the float THEN the original value is preserved
        assert_eq!(parse_info_float(&record, b"FLOAT"), Some(value));
    }
}

#[test]
fn test_parse_info_string_commas() {
    // GIVEN String INFO fields with literal commas, declared with Number=1 and Number=.
    // (rust-htslib splits String values on commas whatever Number is declared)
    let record = info_record("STR=a,b,c;STR_DOT=x,,y");

    // THEN the full value is returned, like cyvcf2 does in scout, including empty parts
    assert_eq!(parse_info_string(&record, b"STR").as_deref(), Some("a,b,c"));
    assert_eq!(
        parse_info_string(&record, b"STR_DOT").as_deref(),
        Some("x,,y")
    );
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
fn test_parse_custom_data_from_record() {
    // GIVEN a SCOUT_CUSTOM annotation with two key/value pairs
    // (from scout's customannotations_one.vcf.gz)
    let record = info_record("SCOUT_CUSTOM=key1|val1,key2|val2");

    // WHEN parsing it as the loader does
    let custom = parse_custom_data(parse_info_string(&record, b"SCOUT_CUSTOM"));

    // THEN both pairs are kept
    assert_eq!(
        custom,
        Some(Bson::Array(vec![
            Bson::Array(vec!["key1".into(), "val1".into()]),
            Bson::Array(vec!["key2".into(), "val2".into()]),
        ]))
    );
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
