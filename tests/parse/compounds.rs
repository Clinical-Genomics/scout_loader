use scout_loader::parse::compounds::parse_compounds;
use scout_loader::parse::ids::parse_document_id;

const CASE_ID: &str = "case_id";

#[test]
fn test_parse_compounds() {
    // GIVEN a compound annotation with two compounds for one family
    let compound_info = Some("internal_id:Y_14923735_C_T>-7|Y_14898429_A_T>-3".to_string());

    // WHEN parsing the compounds
    let compounds = parse_compounds(compound_info, "clinical", CASE_ID);

    // THEN both compounds are returned with their display names and scores
    assert_eq!(compounds.len(), 2);

    assert_eq!(compounds[0].display_name, "Y_14923735_C_T");
    assert_eq!(compounds[0].combined_score, -7.0);

    assert_eq!(compounds[1].display_name, "Y_14898429_A_T");
    assert_eq!(compounds[1].combined_score, -3.0);
}

#[test]
fn test_parse_compound_no_score() {
    // GIVEN a compound without a score
    let compound_info = Some("internal_id:7_117175579_AT_A".to_string());

    // WHEN parsing the compounds
    let compounds = parse_compounds(compound_info, "clinical", CASE_ID);

    // THEN the score defaults to zero
    assert_eq!(compounds.len(), 1);
    assert_eq!(compounds[0].display_name, "7_117175579_AT_A");
    assert_eq!(compounds[0].combined_score, 0.0);
}

#[test]
fn test_parse_compound_no_compound() {
    // GIVEN no compound annotation
    // WHEN parsing the compounds
    let compounds = parse_compounds(None, "clinical", CASE_ID);

    // THEN no compounds are returned
    assert!(compounds.is_empty());
}

#[test]
fn test_parse_compound_without_family_separator() {
    // GIVEN a compound annotation without the family:compounds separator
    let compound_info = Some("7_117175579_AT_A>32".to_string());

    // WHEN parsing the compounds
    let compounds = parse_compounds(compound_info, "clinical", CASE_ID);

    // THEN the entry is skipped
    assert!(compounds.is_empty());
}

#[test]
fn test_parse_compound_id_matches_variant_id() {
    // GIVEN a compound with a chr-prefixed name
    let compound_info = Some("internal_id:chr7_117175579_AT_A>32".to_string());

    // WHEN parsing the compounds for a research variant
    let compounds = parse_compounds(compound_info, "research", CASE_ID);

    // THEN the chr prefix is removed from the display name
    assert_eq!(compounds[0].display_name, "7_117175579_AT_A");

    // THEN the compound points at the document id of the referenced variant
    assert_eq!(
        compounds[0].variant,
        parse_document_id("7", &117175579, "AT", "A", "research", CASE_ID)
    );

    // THEN a clinical compound points at a different variant
    let clinical = parse_compounds(
        Some("internal_id:7_117175579_AT_A>32".to_string()),
        "clinical",
        CASE_ID,
    );
    assert_ne!(clinical[0].variant, compounds[0].variant);
}

#[test]
fn test_parse_compounds_multiple_families() {
    // GIVEN compound annotations for two families
    let compound_info =
        Some("fam1:7_117175579_AT_A>32,fam2:7_117175580_T_G>28|7_117175590_C_G>20".to_string());

    // WHEN parsing the compounds
    let compounds = parse_compounds(compound_info, "clinical", CASE_ID);

    // THEN the compounds of all families are returned
    let names: Vec<&str> = compounds.iter().map(|c| c.display_name.as_str()).collect();
    assert_eq!(
        names,
        vec!["7_117175579_AT_A", "7_117175580_T_G", "7_117175590_C_G"]
    );
}
