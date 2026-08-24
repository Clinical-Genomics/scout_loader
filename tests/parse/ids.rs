use scout_loader::parse::ids::parse_document_id;

#[test]
fn parse_document_id_snv_and_sv_are_different() {
    let snv_id = parse_document_id("1", &123456, "A", "T", "clinical", "case_123");

    let sv_id = parse_document_id("1", &123500, "N", "<DEL>", "clinical", "case_123");

    assert_ne!(snv_id, sv_id);
}
