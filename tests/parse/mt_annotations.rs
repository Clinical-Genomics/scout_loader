use mongodb::bson::Document;

use scout_loader::parse::mt_annotations::{set_hmtvar, set_mitomap_associated_diseases};

use crate::vcf_record::record;

const MT_HEADERS: &[&str] = &[
    "##contig=<ID=MT,length=16569>",
    "##INFO=<ID=MitomapAssociatedDiseases,Number=1,Type=String,Description=\"Mitomap associated diseases\">",
    "##INFO=<ID=HmtVar,Number=1,Type=String,Description=\"HmtVar variant id\">",
];

/// Set the MT annotations of a record with the given INFO column on a new variant document.
fn mt_annotations(info: &str) -> Document {
    let record = record(MT_HEADERS, &format!("MT\t3243\t.\tA\tG\t.\tPASS\t{info}"));
    let mut variant = Document::new();
    set_mitomap_associated_diseases(&record, &mut variant);
    set_hmtvar(&record, &mut variant);
    variant
}

#[test]
fn test_set_mitomap_associated_diseases() {
    // GIVEN a Mitomap disease annotation with underscores
    let variant = mt_annotations("MitomapAssociatedDiseases=MELAS_/_Leigh_Syndrome");

    // THEN the underscores are replaced with spaces
    assert_eq!(
        variant.get_str("mitomap_associated_diseases"),
        Ok("MELAS / Leigh Syndrome")
    );
}

#[test]
fn test_set_hmtvar() {
    // GIVEN a HmtVar variant id
    let variant = mt_annotations("HmtVar=2963");

    // THEN it is stored as an integer
    assert_eq!(variant.get_i32("hmtvar_variant_id"), Ok(2963));
}

#[test]
fn test_set_hmtvar_not_numeric() {
    // GIVEN a non-numeric HmtVar value
    let variant = mt_annotations("HmtVar=unknown");

    // THEN no HmtVar id is set
    assert!(!variant.contains_key("hmtvar_variant_id"));
}

#[test]
fn test_mt_annotations_missing() {
    // GIVEN a record without MT annotations
    let variant = mt_annotations(".");

    // THEN no MT annotations are set
    assert!(variant.is_empty());
}
