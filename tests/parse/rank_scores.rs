use mongodb::bson::{Bson, doc};
use scout_loader::parse::rank_scores::{
    parse_rank_result, parse_rank_score_other, parse_rank_scores, parse_score_entry,
};

#[test]
fn test_parse_rank_scores() {
    // GIVEN a variant with a rank score and a normalized rank score
    let record = rank_record("RankScore=internal_id:-20;RankScoreNormalized=internal_id:-0.5");

    // WHEN parsing the rank scores
    let (rank_score, norm_rank_score) = parse_rank_scores(&record);

    // THEN both scores are parsed
    assert_eq!(rank_score, Some(-20.0));
    assert_eq!(norm_rank_score, Some(-0.5));
}

#[test]
fn test_parse_score_entry() {
    assert_eq!(parse_score_entry(Some("internal_id:-20")), Some("-20"));
    assert_eq!(parse_score_entry(Some("internal_id:12.5")), Some("12.5"));
    assert_eq!(parse_score_entry(None), None);
    assert_eq!(parse_score_entry(Some("invalid")), None);
}

fn rank_record(info: &str) -> rust_htslib::bcf::Record {
    crate::vcf_record::record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=RankScore,Number=.,Type=String,Description=\"Rank score\">",
            "##INFO=<ID=RankScoreNormalized,Number=.,Type=String,Description=\"Normalized rank score\">",
            "##INFO=<ID=RankResult,Number=.,Type=String,Description=\"Rank result\">",
            "##INFO=<ID=MivmirScore,Number=1,Type=Float,Description=\"Mivmir score\">",
            "##INFO=<ID=MivmirExplanation,Number=1,Type=String,Description=\"Mivmir explanation\">",
            "##INFO=<ID=GicamScore,Number=1,Type=Float,Description=\"Gicam score\">",
        ],
        &format!("1\t100\t.\tA\tT\t.\tPASS\t{info}"),
    )
}

#[test]
fn test_parse_rank_scores_missing_or_invalid() {
    // GIVEN no rank scores, and rank scores without a score part or a number
    let (missing, missing_norm) = parse_rank_scores(&rank_record("."));
    let (invalid, invalid_norm) = parse_rank_scores(&rank_record(
        "RankScore=case_1;RankScoreNormalized=case_1:abc",
    ));

    // THEN None is returned
    assert_eq!((missing, missing_norm), (None, None));
    assert_eq!((invalid, invalid_norm), (None, None));
}

#[test]
fn test_parse_rank_result() {
    // GIVEN a RankResult with a header of categories
    let header: Vec<String> = ["Consequence", "Frequency", "Inheritance_Models"]
        .iter()
        .map(|category| category.to_string())
        .collect();
    let record = rank_record("RankResult=5|-1|x");

    // WHEN parsing the rank result
    let parsed = parse_rank_result(&record, &header).unwrap();

    // THEN numeric scores are paired with their category, invalid scores skipped
    assert_eq!(
        parsed,
        vec![
            Bson::Document(doc! { "category": "Consequence", "score": 5 }),
            Bson::Document(doc! { "category": "Frequency", "score": -1 }),
        ]
    );

    // AND a missing RankResult gives None
    assert_eq!(parse_rank_result(&rank_record("."), &header), None);
}

#[test]
fn test_parse_rank_score_mivmir_gicam() {
    // GIVEN a variant with Mivmir (with explanation) and Gicam scores
    let record = rank_record(
        "MivmirScore=0.0045;GicamScore=0.076;MivmirExplanation=[CSQ_CLINVAR_CLNSIG=0.48]",
    );

    // WHEN parsing other rank scores
    let parsed = parse_rank_score_other(&record).unwrap();

    // THEN both scores are parsed (VCF Floats are 32-bit)
    assert_eq!(
        parsed.get_document("Gicam").unwrap(),
        &doc! { "value": f64::from(0.076_f32) }
    );
    let mivmir = parsed.get_document("Mivmir").unwrap();
    assert_eq!(mivmir.get_f64("value"), Ok(f64::from(0.0045_f32)));
    assert_eq!(
        mivmir.get_document("desc"),
        Ok(&doc! { "CSQ_CLINVAR_CLNSIG": 0.48 })
    );
}

#[test]
fn test_parse_mivmir_explanation_with_several_entries() {
    // GIVEN a Mivmir explanation with several comma-separated entries
    let record = rank_record(
        "MivmirScore=0.0045;MivmirExplanation=[CSQ_CLINVAR_CLNSIG=0.48,CSQ_CLINVAR_CLNREVSTAT=0.24,]",
    );

    // WHEN parsing other rank scores
    let parsed = parse_rank_score_other(&record).unwrap();

    // THEN all entries are part of the description
    assert_eq!(
        parsed.get_document("Mivmir").unwrap().get_document("desc"),
        Ok(&doc! { "CSQ_CLINVAR_CLNSIG": 0.48, "CSQ_CLINVAR_CLNREVSTAT": 0.24 })
    );
}

#[test]
fn test_parse_rank_score_other_missing() {
    // GIVEN a variant without other rank scores
    // THEN None is returned
    assert_eq!(parse_rank_score_other(&rank_record(".")), None);
}
