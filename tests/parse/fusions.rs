use mongodb::bson::{Bson, Document, doc};
use rust_htslib::bcf::Record;

use scout_loader::parse::fusions::set_fusion_info;

use crate::vcf_record::record;

/// INFO headers as written by nf-core/rnafusion (see scout's `demo/fusion_data.vcf.gz`).
const FUSION_HEADERS: &[&str] = &[
    "##contig=<ID=1,length=249250621>",
    "##contig=<ID=2,length=243199373>",
    "##INFO=<ID=GENEA,Number=.,Type=String,Description=\"Gene A\">",
    "##INFO=<ID=GENEB,Number=.,Type=String,Description=\"Gene B\">",
    "##INFO=<ID=ORIENTATION,Number=.,Type=String,Description=\"Strand1 and strand2 directions\">",
    "##INFO=<ID=FOUND_DB,Number=.,Type=String,Description=\"Databases in which the fusion has been found\">",
    "##INFO=<ID=TOOL_HITS,Number=.,Type=Integer,Description=\"Number of tools that found the fusion\">",
    "##INFO=<ID=SCORE,Number=.,Type=Float,Description=\"Score from fusionreport\">",
    "##INFO=<ID=FRAME_STATUS,Number=.,Type=String,Description=\"Frame status of the fusion\">",
    "##INFO=<ID=TRANSCRIPT_ID_A,Number=.,Type=String,Description=\"Transcript id A \">",
    "##INFO=<ID=TRANSCRIPT_ID_B,Number=.,Type=String,Description=\"Transcript id B\">",
    "##INFO=<ID=HGNC_ID_A,Number=.,Type=Float,Description=\"HGNC id A\">",
    "##INFO=<ID=HGNC_ID_B,Number=.,Type=Float,Description=\"HGNC id A\">",
    "##INFO=<ID=EXON_NUMBER_A,Number=.,Type=Float,Description=\"Exon number A\">",
    "##INFO=<ID=EXON_NUMBER_B,Number=.,Type=Float,Description=\"Exon number B\">",
];

/// Build a fusion record with the given INFO column.
fn fusion_record(info: &str) -> Record {
    record(
        FUSION_HEADERS,
        &format!("2\t42295516\t.\tN\tN]2:29223528]\t.\tPASS\t{info}"),
    )
}

/// Parse the fusion info of a record into a new variant document.
fn fusion_info(info: &str) -> Document {
    let mut variant = Document::new();
    set_fusion_info(&fusion_record(info), &mut variant);
    variant
}

#[test]
fn test_set_fusion_info() {
    // GIVEN an in-frame EML4-ALK fusion found by three tools and in four databases
    let variant = fusion_info(
        "GENEA=EML4;GENEB=ALK;ORIENTATION=+,-;FOUND_DB=COSMIC, Mitelman, FusionGDB, FusionGDB2;\
         TOOL_HITS=3;SCORE=0.688;FRAME_STATUS=INFRAME;TRANSCRIPT_ID_A=ENST00000401738;\
         TRANSCRIPT_ID_B=ENST00000389048;HGNC_ID_A=1316;HGNC_ID_B=427;EXON_NUMBER_A=14;\
         EXON_NUMBER_B=20",
    );

    // THEN the fusion annotations are set
    assert_eq!(variant.get_i32("tool_hits"), Ok(3));
    assert_eq!(
        variant.get_array("found_db").unwrap(),
        &vec![
            Bson::from("COSMIC"),
            Bson::from("Mitelman"),
            Bson::from("FusionGDB"),
            Bson::from("FusionGDB2"),
        ]
    );
    assert!((variant.get_f64("fusion_score").unwrap() - 0.688).abs() < 1e-6);
    assert_eq!(
        variant.get_f64("rank_score"),
        variant.get_f64("fusion_score")
    );
    assert_eq!(variant.get_str("orientation"), Ok("+,-"));
    assert_eq!(variant.get_str("frame_status"), Ok("INFRAME"));

    // THEN both partners are added as genes with their transcript and exon
    assert_eq!(
        variant.get_array("genes").unwrap(),
        &vec![
            Bson::Document(doc! {
                "hgnc_symbol": "EML4",
                "hgnc_id": 1316,
                "transcripts": [{
                    "transcript_id": "ENST00000401738",
                    "hgnc_id": 1316,
                    "exon": "14",
                    "is_canonical": false,
                }],
            }),
            Bson::Document(doc! {
                "hgnc_symbol": "ALK",
                "hgnc_id": 427,
                "transcripts": [{
                    "transcript_id": "ENST00000389048",
                    "hgnc_id": 427,
                    "exon": "20",
                    "is_canonical": false,
                }],
            }),
        ]
    );
    assert_eq!(
        variant.get_array("hgnc_ids").unwrap(),
        &vec![Bson::Int32(1316), Bson::Int32(427)]
    );
    assert_eq!(
        variant.get_array("hgnc_symbols").unwrap(),
        &vec![Bson::from("EML4"), Bson::from("ALK")]
    );
    assert_eq!(
        variant.get_array("fusion_genes").unwrap(),
        &vec![Bson::from("EML4"), Bson::from("ALK")]
    );
}

#[test]
fn test_set_fusion_info_without_transcripts() {
    // GIVEN a fusion without database hits, frame status or transcript information
    let variant = fusion_info(
        "GENEA=VPS72;GENEB=TMOD4;HGNC_ID_A=11644;HGNC_ID_B=11874;ORIENTATION=-,-;\
         FOUND_DB=[];TOOL_HITS=1;SCORE=0.062",
    );

    // THEN an empty database list is not stored
    assert!(!variant.contains_key("found_db"));

    // THEN a missing frame status is stored as an empty string
    assert_eq!(variant.get_str("frame_status"), Ok(""));

    // THEN the partners are added as genes without transcripts
    let genes = variant.get_array("genes").unwrap();
    assert_eq!(genes.len(), 2);
    for gene in genes {
        assert!(
            gene.as_document()
                .unwrap()
                .get_array("transcripts")
                .unwrap()
                .is_empty()
        );
    }
}

#[test]
fn test_set_fusion_info_nan_values() {
    // GIVEN a fusion with nan placeholders and an exon number of 0 for partner A
    let variant = fusion_info(
        "GENEA=PAX3;GENEB=FOXO1;HGNC_ID_A=8617;HGNC_ID_B=3819;ORIENTATION=nan,nan;\
         FRAME_STATUS=nan;TRANSCRIPT_ID_A=nan;TRANSCRIPT_ID_B=ENST00000379561;\
         EXON_NUMBER_A=0;EXON_NUMBER_B=2;TOOL_HITS=2;SCORE=0.5",
    );

    // THEN the nan placeholders are stored as empty strings
    assert_eq!(variant.get_str("orientation"), Ok(""));
    assert_eq!(variant.get_str("frame_status"), Ok(""));

    // THEN partner A gets no transcript, partner B keeps its transcript and exon
    let genes = variant.get_array("genes").unwrap();
    let transcripts = |index: usize| {
        genes[index]
            .as_document()
            .unwrap()
            .get_array("transcripts")
            .unwrap()
            .clone()
    };
    assert!(transcripts(0).is_empty());
    assert_eq!(
        transcripts(1),
        vec![Bson::Document(doc! {
            "transcript_id": "ENST00000379561",
            "hgnc_id": 3819,
            "exon": "2",
            "is_canonical": false,
        })]
    );
}

#[test]
fn test_set_fusion_info_partner_without_hgnc_id() {
    // GIVEN a fusion where partner B has no HGNC id
    let variant = fusion_info("GENEA=EML4;GENEB=LINC01234;HGNC_ID_A=1316;TOOL_HITS=1;SCORE=0.1");

    // THEN only partner A is added as a gene
    assert_eq!(
        variant.get_array("hgnc_ids").unwrap(),
        &vec![Bson::Int32(1316)]
    );
    assert_eq!(
        variant.get_array("hgnc_symbols").unwrap(),
        &vec![Bson::from("EML4")]
    );
    assert_eq!(variant.get_array("genes").unwrap().len(), 1);

    // THEN both partners are listed as fusion genes
    assert_eq!(
        variant.get_array("fusion_genes").unwrap(),
        &vec![Bson::from("EML4"), Bson::from("LINC01234")]
    );
}
