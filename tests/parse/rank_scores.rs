use rust_htslib::bcf::{Read, Reader};
use scout_loader::parse::rank_scores::{parse_rank_scores, parse_score_entry};
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_parse_rank_scores() {
    let vcf = b"##fileformat=VCFv4.2
##INFO=<ID=RankScore,Number=1,Type=String,Description=\"Rank score\">
##INFO=<ID=RankScoreNormalized,Number=1,Type=String,Description=\"Normalized rank score\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tT\t.\t.\tRankScore=internal_id:-20;RankScoreNormalized=internal_id:-0.5
";

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(vcf).unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();
    let mut record = reader.empty_record();

    let _ = reader.read(&mut record).unwrap();

    let (rank_score, norm_rank_score) = parse_rank_scores(&record);

    assert_eq!(rank_score, -20.0);
    assert_eq!(norm_rank_score, -0.5);
}

#[test]
fn test_parse_score_entry() {
    assert_eq!(parse_score_entry(Some("internal_id:-20")), Some("-20"));
    assert_eq!(parse_score_entry(Some("internal_id:12.5")), Some("12.5"));
    assert_eq!(parse_score_entry(None), None);
    assert_eq!(parse_score_entry(Some("invalid")), None);
}
