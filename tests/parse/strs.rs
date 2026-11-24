//! STR annotation tests for old and new stranger output (stranger#124).
//!
//! Older stranger versions write catalog values with Python's `str()`, e.g.
//! `PathologicStruc=[0, 1, 2, 6]` and `SourceDisplay=GeneReviews Internet 2019-11-07`.
//! Newer versions write `PathologicStruc` as an Integer list and percent-encode
//! free text. Both must give the same variant document.

use crate::vcf_record::record;
use mongodb::bson::{Document, doc};
use rust_htslib::bcf::Record;
use scout_loader::parse::strs::set_str_info;

const OLD_STRUC_HEADER: &str = "##INFO=<ID=PathologicStruc,Number=1,Type=String,Description=\"Struc components counting towards pathogenicity as list of motif numbers\">";
const NEW_STRUC_HEADER: &str = "##INFO=<ID=PathologicStruc,Number=.,Type=Integer,Description=\"Struc components counting towards pathogenicity as list of motif numbers\">";

/// Build a stranger-annotated STR record with the given PathologicStruc header and INFO.
fn str_record(pathologic_struc_header: &str, info: &str) -> Record {
    record(
        &[
            "##contig=<ID=X,length=156040895>",
            pathologic_struc_header,
            "##INFO=<ID=SourceDisplay,Number=1,Type=String,Description=\"Source for variant definition, display\">",
            "##INFO=<ID=Source,Number=1,Type=String,Description=\"Source collection for variant definition\">",
            "##INFO=<ID=SourceId,Number=1,Type=String,Description=\"Source id for variant definition\">",
            "##INFO=<ID=Disease,Number=1,Type=String,Description=\"Associated disorder\">",
            "##INFO=<ID=InheritanceMode,Number=1,Type=String,Description=\"Main mode of inheritance for disorder\">",
            "##INFO=<ID=DisplayRU,Number=1,Type=String,Description=\"Display repeat unit familiar to clinician\">",
        ],
        &format!("X\t147912051\t.\tC\t<STR>\t.\tPASS\t{info}"),
    )
}

fn str_info(record: &Record) -> Document {
    let mut variant = Document::new();
    set_str_info(record, &mut variant);
    variant
}

fn expected() -> Document {
    doc! {
        "str_display_ru": "CGG",
        "str_pathologic_struc": "[0, 1, 2, 6]",
        "str_disease": "Fragile X syndrome",
        "str_inheritance_mode": "XLD",
        "str_source": {
            "display": "GeneReviews Internet 2019-11-07",
            "type": "GeneReviews",
            "id": "NBK1384",
        },
    }
}

#[test]
fn test_set_str_info_old_stranger_format() {
    // GIVEN STR annotations as written by older stranger versions, with spaces
    let record = str_record(
        OLD_STRUC_HEADER,
        "PathologicStruc=[0, 1, 2, 6];SourceDisplay=GeneReviews Internet 2019-11-07;Source=GeneReviews;SourceId=NBK1384;Disease=Fragile X syndrome;InheritanceMode=XLD;DisplayRU=CGG",
    );

    // THEN the values are stored as before
    assert_eq!(str_info(&record), expected());
}

#[test]
fn test_set_str_info_new_stranger_format() {
    // GIVEN the same annotations as written after stranger#124:
    // an Integer list and percent-encoded free text
    let record = str_record(
        NEW_STRUC_HEADER,
        "PathologicStruc=0,1,2,6;SourceDisplay=GeneReviews%20Internet%202019-11-07;Source=GeneReviews;SourceId=NBK1384;Disease=Fragile%20X%20syndrome;InheritanceMode=XLD;DisplayRU=CGG",
    );

    // THEN the stored values are the same as for the old format
    assert_eq!(str_info(&record), expected());
}

#[test]
fn test_set_str_info_single_pathologic_struc() {
    // GIVEN a single pathologic motif index in both formats
    let old = str_record(OLD_STRUC_HEADER, "PathologicStruc=[0]");
    let new = str_record(NEW_STRUC_HEADER, "PathologicStruc=0");

    // THEN both are stored as a one-element list
    assert_eq!(str_info(&old), doc! { "str_pathologic_struc": "[0]" });
    assert_eq!(str_info(&new), doc! { "str_pathologic_struc": "[0]" });
}

#[test]
fn test_set_str_info_percent_encoded_comma() {
    // GIVEN a source display with encoded commas and spaces
    let record = str_record(
        NEW_STRUC_HEADER,
        "SourceDisplay=Paucar%2C%20Nilsson%2C%20Engvall%20et%20al%20medRxiv%202023-10-03",
    );

    // THEN it is decoded for display
    assert_eq!(
        str_info(&record)
            .get_document("str_source")
            .unwrap()
            .get_str("display"),
        Ok("Paucar, Nilsson, Engvall et al medRxiv 2023-10-03")
    );
}
