use crate::vcf_record::record;
use mongodb::bson::Document;
use scout_loader::parse::strs::set_str_info;

#[test]
fn test_set_str_info_with_commas() {
    // GIVEN TRGT/stranger STR annotations containing commas
    // (from scout's NIST.trgt.stranger.vcf.gz)
    let record = record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=MOTIFS,Number=.,Type=String,Description=\"Motifs that the tandem repeat is composed of\">",
            "##INFO=<ID=PathologicStruc,Number=1,Type=String,Description=\"Struc components counting towards pathogenicity\">",
            "##INFO=<ID=SourceDisplay,Number=1,Type=String,Description=\"Source for variant definition, display\">",
        ],
        "1\t100\t.\tA\tG\t.\tPASS\tMOTIFS=AAAAT,GAAAT;PathologicStruc=[0, 1, 2, 6];SourceDisplay=Paucar, Nilsson, Engvall et al medRxiv 2023-10-03",
    );
    let mut variant = Document::new();

    // WHEN setting the STR info
    set_str_info(&record, &mut variant);

    // THEN the full values are stored, as scout stores them
    assert_eq!(variant.get_str("str_motifs"), Ok("AAAAT,GAAAT"));
    assert_eq!(variant.get_str("str_pathologic_struc"), Ok("[0, 1, 2, 6]"));
    assert_eq!(
        variant
            .get_document("str_source")
            .unwrap()
            .get_str("display"),
        Ok("Paucar, Nilsson, Engvall et al medRxiv 2023-10-03")
    );
}
