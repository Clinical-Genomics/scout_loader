use scout_loader::parse::compounds::parse_compounds;

#[test]
fn test_parse_compounds() {
    let compound_info = Some("internal_id:Y_14923735_C_T>-7|Y_14898429_A_T>-3".to_string());

    let compounds = parse_compounds(compound_info, "clinical");

    assert_eq!(compounds.len(), 2);

    assert_eq!(compounds[0].display_name, "Y_14923735_C_T");
    assert_eq!(compounds[0].combined_score, -7.0);

    assert_eq!(compounds[1].display_name, "Y_14898429_A_T");
    assert_eq!(compounds[1].combined_score, -3.0);
}
