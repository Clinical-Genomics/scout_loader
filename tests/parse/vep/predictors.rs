use mongodb::bson::Document;
use rust_htslib::bcf::{Read, Reader};
use scout_loader::parse::severity::set_severity_predictions;
use scout_loader::parse::vep::transcripts::parse_vep_transcripts;
use std::collections::HashMap;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_parse_alphamissense() {
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: FEATURE|AM_CLASS|AM_PATHOGENICITY\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t6448812\tchr1_6448812_C_T\tC\tT\t57\tPASS\tCSQ=ENST00000632593|likely_benign|0.1248,ENST00000636330|likely_benign|0.1248,ENST00000645284|likely_benign|0.1248
";

    let mut temp_file = NamedTempFile::new().unwrap();
    temp_file.write_all(vcf.as_bytes()).unwrap();

    let mut reader = Reader::from_path(temp_file.path()).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    let field_indices = HashMap::from([
        ("FEATURE".to_string(), 0),
        ("AM_CLASS".to_string(), 1),
        ("AM_PATHOGENICITY".to_string(), 2),
    ]);

    let mut variant = Document::new();

    let parsed_transcripts = parse_vep_transcripts(&record, &field_indices, &mut variant);

    assert_eq!(parsed_transcripts.len(), 3);

    for transcript in &parsed_transcripts {
        assert_eq!(transcript.get_str("am_class").unwrap(), "likely_benign");
        assert_eq!(transcript.get_f64("am_pathogenicity").unwrap(), 0.1248);
    }

    set_severity_predictions(&mut variant, &record, &parsed_transcripts);

    assert_eq!(variant.get_str("am_class").unwrap(), "likely_benign");
    assert_eq!(variant.get_f64("am_pathogenicity").unwrap(), 0.1248);
}
