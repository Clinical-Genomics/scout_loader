use scout_loader::parse::genetic_models::parse_genetic_models;

use crate::vcf_record::record;

const GENETIC_MODELS_HEADER: &str =
    "##INFO=<ID=GeneticModels,Number=1,Type=String,Description=\"Genetic models\">";

#[test]
fn test_parse_genetic_models() {
    // GIVEN a record with two genetic models for one family
    let record = record(
        &[GENETIC_MODELS_HEADER],
        "1\t100\t.\tA\tT\t.\t.\tGeneticModels=internal_id:XD_dn|XR_dn",
    );

    // WHEN parsing the genetic models
    let genetic_models = parse_genetic_models(&record);

    // THEN both models are returned
    assert_eq!(
        genetic_models,
        vec!["XD_dn".to_string(), "XR_dn".to_string()]
    );
}

#[test]
fn test_parse_genetic_models_missing() {
    // GIVEN a record without genetic models
    let record = record(&[GENETIC_MODELS_HEADER], "1\t100\t.\tA\tT\t.\t.\t.");

    // THEN no models are returned
    assert!(parse_genetic_models(&record).is_empty());
}
