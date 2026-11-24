use scout_loader::parse::clinvar_pathogenic::{is_pathogenic_info, is_pathogenic_transcripts};

use crate::vcf_record::{record, vep_record};

const CLNSIG_HEADER: &str =
    "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">";

#[test]
fn test_is_pathogenic_transcripts() {
    // GIVEN a transcript with a conflicting ClinVar classification in CLINVAR_CLNSIG
    let (record, field_indices) = vep_record(
        "Allele|CLNSIG|CLINVAR_CLNSIG",
        "T|Benign|Conflicting_classifications_of_pathogenicity",
        &[],
        "",
    );

    // THEN the transcripts count as pathogenic
    assert!(is_pathogenic_transcripts(&record, &field_indices));
}

#[test]
fn test_is_pathogenic_transcripts_benign() {
    // GIVEN transcripts with only benign ClinVar classifications
    let (record, field_indices) = vep_record("Allele|CLNSIG", "T|Benign,T|Likely_benign", &[], "");

    // THEN the transcripts do not count as pathogenic
    assert!(!is_pathogenic_transcripts(&record, &field_indices));
}

#[test]
fn test_is_pathogenic_transcripts_second_transcript() {
    // GIVEN a benign first transcript and a likely pathogenic second transcript
    let (record, field_indices) =
        vep_record("Allele|CLNSIG", "T|Benign,T|likely_pathogenic", &[], "");

    // THEN the transcripts count as pathogenic
    assert!(is_pathogenic_transcripts(&record, &field_indices));
}

#[test]
fn test_is_pathogenic_transcripts_without_clnsig_fields() {
    // GIVEN a CSQ without ClinVar fields
    let (record, field_indices) = vep_record("Allele|SYMBOL", "T|Pathogenic", &[], "");

    // THEN the transcripts do not count as pathogenic
    assert!(!is_pathogenic_transcripts(&record, &field_indices));
}

#[test]
fn test_is_pathogenic_info() {
    // GIVEN CLNSIG INFO values, including multi-valued and numeric ones
    let is_pathogenic = |clnsig: &str| {
        let record = record(
            &[CLNSIG_HEADER],
            &format!("1\t100\t.\tA\tG\t.\tPASS\tCLNSIG={clnsig}"),
        );
        is_pathogenic_info(&record)
    };

    // THEN pathogenic, likely pathogenic, conflicting and 4/5 values count as pathogenic
    assert!(is_pathogenic("Pathogenic"));
    assert!(is_pathogenic("Benign/Likely_pathogenic"));
    assert!(is_pathogenic(
        "Conflicting_interpretations_of_pathogenicity"
    ));
    assert!(is_pathogenic("2|5"));

    // THEN benign and uncertain values do not
    assert!(!is_pathogenic("Benign"));
    assert!(!is_pathogenic("Uncertain_significance"));
}

#[test]
fn test_is_pathogenic_info_missing() {
    // GIVEN a record without CLNSIG
    let record = record(&[CLNSIG_HEADER], "1\t100\t.\tA\tG\t.\tPASS\t.");

    // THEN the variant does not count as pathogenic
    assert!(!is_pathogenic_info(&record));
}
