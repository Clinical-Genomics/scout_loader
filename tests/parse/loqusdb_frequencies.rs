use mongodb::bson::{Document, doc};

use scout_loader::parse::loqusdb_frequencies::add_loqus_archive_frequencies;

use crate::vcf_record::record;

const LOQUSDB_HEADERS: &[&str] = &[
    "##INFO=<ID=Obs,Number=1,Type=Integer,Description=\"The number of observations for the variant\">",
    "##INFO=<ID=Hom,Number=1,Type=Integer,Description=\"The number of observed homozygotes\">",
    "##INFO=<ID=Frq,Number=1,Type=Float,Description=\"Observed frequency of the variant\">",
    "##INFO=<ID=clinical_genomics_loqusObs,Number=1,Type=Integer,Description=\"Observations\">",
    "##INFO=<ID=clinical_genomics_loqusFrq,Number=1,Type=Float,Description=\"Frequency\">",
    "##INFO=<ID=clin_obs,Number=1,Type=Integer,Description=\"Observations\">",
    "##INFO=<ID=Cancer_Germline_Obs,Number=1,Type=Integer,Description=\"Cancer germline observations\">",
    "##INFO=<ID=Cancer_Germline_Hom,Number=1,Type=Integer,Description=\"Cancer germline homozygotes\">",
    "##INFO=<ID=Cancer_Germline_Frq,Number=1,Type=Float,Description=\"Cancer germline frequency\">",
    "##INFO=<ID=Cancer_Somatic_Obs,Number=1,Type=Integer,Description=\"Cancer somatic observations\">",
    "##INFO=<ID=Cancer_Somatic_Panel_Obs,Number=1,Type=Integer,Description=\"Cancer somatic panel observations\">",
];

/// Add the LoqusDB archive frequencies of a record with the given INFO column to a new variant.
fn loqus_frequencies(info: &str, local_archive_info: Option<&Document>) -> Document {
    let record = record(
        LOQUSDB_HEADERS,
        &format!("1\t100\t.\tA\tG\t.\tPASS\t{info}"),
    );
    let mut variant = Document::new();
    add_loqus_archive_frequencies(&record, &mut variant, local_archive_info);
    variant
}

#[test]
fn test_add_loqus_archive_frequencies() {
    // GIVEN germline observations, homozygotes and frequency, and archive metadata
    let archive = doc! {
        "Description": "Local germline archive",
        "Date": "2026-11-01",
        "NrCases": 1200,
    };

    // WHEN adding the archive frequencies
    let variant = loqus_frequencies("Obs=12;Hom=2;Frq=0.01", Some(&archive));

    // THEN the observations and archive metadata are set
    assert_eq!(variant.get_i32("local_obs_old"), Ok(12));
    assert_eq!(variant.get_i32("local_obs_hom_old"), Ok(2));
    assert!((variant.get_f64("local_obs_old_freq").unwrap() - 0.01).abs() < 1e-6);
    assert_eq!(
        variant.get_str("local_obs_old_desc"),
        Ok("Local germline archive")
    );
    assert_eq!(variant.get_str("local_obs_old_date"), Ok("2026-11-01"));
    assert_eq!(variant.get_i32("local_obs_old_nr_cases"), Ok(1200));
}

#[test]
fn test_add_loqus_archive_frequencies_clinical_genomics_keys() {
    // GIVEN observations under the clinical_genomics_loqus keys
    // WHEN adding the archive frequencies
    let variant = loqus_frequencies(
        "clinical_genomics_loqusObs=7;clinical_genomics_loqusFrq=0.002",
        None,
    );

    // THEN they are used as the old local observations
    assert_eq!(variant.get_i32("local_obs_old"), Ok(7));
    assert!((variant.get_f64("local_obs_old_freq").unwrap() - 0.002).abs() < 1e-6);
}

#[test]
fn test_add_loqus_archive_frequencies_sentinel_falls_back() {
    // GIVEN Obs set to the -1 sentinel and a clin_obs count
    // WHEN adding the archive frequencies
    let variant = loqus_frequencies("Obs=-1;Hom=-1;clin_obs=4", None);

    // THEN the sentinel is skipped and clin_obs is used
    assert_eq!(variant.get_i32("local_obs_old"), Ok(4));
    assert!(!variant.contains_key("local_obs_hom_old"));
}

#[test]
fn test_add_loqus_archive_frequencies_cancer() {
    // GIVEN cancer germline, somatic and somatic panel observations
    // WHEN adding the archive frequencies
    let variant = loqus_frequencies(
        "Cancer_Germline_Obs=3;Cancer_Germline_Hom=1;Cancer_Germline_Frq=0.05;\
         Cancer_Somatic_Obs=9;Cancer_Somatic_Panel_Obs=-1",
        None,
    );

    // THEN each cancer archive gets its own keys
    assert_eq!(variant.get_i32("local_obs_cancer_germline_old"), Ok(3));
    assert_eq!(variant.get_i32("local_obs_cancer_germline_hom_old"), Ok(1));
    assert!(
        (variant
            .get_f64("local_obs_cancer_germline_old_freq")
            .unwrap()
            - 0.05)
            .abs()
            < 1e-6
    );
    assert_eq!(variant.get_i32("local_obs_cancer_somatic_old"), Ok(9));

    // THEN sentinel values are not stored
    assert!(!variant.contains_key("local_obs_cancer_somatic_panel_old"));

    // THEN no germline keys are set
    assert!(!variant.contains_key("local_obs_old"));
}

#[test]
fn test_add_loqus_archive_frequencies_missing() {
    // GIVEN a record without LoqusDB annotations and no archive metadata
    // WHEN adding the archive frequencies
    let variant = loqus_frequencies(".", None);

    // THEN nothing is added
    assert!(variant.is_empty());
}
