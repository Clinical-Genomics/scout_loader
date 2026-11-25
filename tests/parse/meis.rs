use mongodb::bson::{Document, doc};

use scout_loader::parse::meis::set_mei_info;

use crate::vcf_record::record;

const MEINFO_HEADER: &str = "##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">";
const MEINFO_ANY_NUMBER_HEADER: &str = "##INFO=<ID=MEINFO,Number=.,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">";

/// Parse the MEI info of a record with the given MEINFO header and INFO column.
fn mei_info(header: &str, info: &str) -> Document {
    let record = record(
        &[header],
        &format!("1\t21898411\t.\tA\t<INS:ME>\t18\t.\t{info}"),
    );
    let mut variant = Document::new();
    set_mei_info(&record, &mut variant);
    variant
}

#[test]
fn test_set_mei_info() {
    // GIVEN an ALU insertion with name, start, end and polarity
    let variant = mei_info(MEINFO_HEADER, "SVTYPE=INS;MEINFO=ALU,21898411,21898412,+");

    // THEN the name and polarity are set
    assert_eq!(variant, doc! { "mei_name": "ALU", "mei_polarity": "+" });
}

#[test]
fn test_set_mei_info_unknown_polarity() {
    // GIVEN an insertion with an unknown (NA) polarity, as written by retroseq
    let variant = mei_info(MEINFO_HEADER, "MEINFO=ALU,21898411,21898412,NA");

    // THEN the polarity is kept as given
    assert_eq!(variant.get_str("mei_polarity"), Ok("NA"));
}

#[test]
fn test_set_mei_info_malformed() {
    // GIVEN a MEINFO with three instead of four fields
    let variant = mei_info(MEINFO_ANY_NUMBER_HEADER, "MEINFO=ALU,21898411,21898412");

    // THEN no MEI info is set
    assert!(variant.is_empty());
}

#[test]
fn test_set_mei_info_missing() {
    // GIVEN a record without MEINFO
    let variant = mei_info(MEINFO_HEADER, ".");

    // THEN no MEI info is set
    assert!(variant.is_empty());
}
