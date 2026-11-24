use mongodb::bson::Document;
use scout_loader::parse::severity::set_severity_predictions;
use scout_loader::parse::vep::transcripts::parse_vep_transcripts;

use crate::vcf_record::vep_record;

#[test]
fn test_parse_alphamissense() {
    // GIVEN AlphaMissense predictions on three transcripts
    let (record, field_indices) = vep_record(
        "FEATURE|AM_CLASS|AM_PATHOGENICITY",
        "ENST00000632593|likely_benign|0.1248,ENST00000636330|likely_benign|0.1248,ENST00000645284|likely_benign|0.1248",
        &[],
        "",
    );

    let mut variant = Document::new();

    // WHEN parsing the transcripts
    let parsed_transcripts = parse_vep_transcripts(&record, &field_indices, &mut variant);

    // THEN each transcript gets the AlphaMissense class and score
    assert_eq!(parsed_transcripts.len(), 3);

    for transcript in &parsed_transcripts {
        assert_eq!(transcript.get_str("am_class").unwrap(), "likely_benign");
        assert_eq!(transcript.get_f64("am_pathogenicity").unwrap(), 0.1248);
    }

    // WHEN setting severity predictions on the variant
    set_severity_predictions(&mut variant, &record, &parsed_transcripts);

    // THEN the variant gets the AlphaMissense prediction
    assert_eq!(variant.get_str("am_class").unwrap(), "likely_benign");
    assert_eq!(variant.get_f64("am_pathogenicity").unwrap(), 0.1248);
}

#[test]
fn test_parse_avi_scores() {
    // GIVEN a VCF with AVI_PHRED and AVI_RAW in the CSQ field, including a
    // transcript with multiple values, a negative raw score and a transcript
    // without scores
    let (record, field_indices) = vep_record(
        "FEATURE|AVI_PHRED|AVI_RAW",
        "ENST00000632593|12.5|-0.75,ENST00000636330|21.3&18.1|1.92&2.05,ENST00000645284||",
        &[],
        "",
    );

    let mut variant = Document::new();

    // WHEN parsing the transcripts
    let parsed_transcripts = parse_vep_transcripts(&record, &field_indices, &mut variant);

    // THEN each transcript gets the highest of its AVI scores
    assert_eq!(parsed_transcripts.len(), 3);
    assert_eq!(parsed_transcripts[0].get_f64("avi_phred").unwrap(), 12.5);
    assert_eq!(parsed_transcripts[0].get_f64("avi_raw").unwrap(), -0.75);
    assert_eq!(parsed_transcripts[1].get_f64("avi_phred").unwrap(), 21.3);
    assert_eq!(parsed_transcripts[1].get_f64("avi_raw").unwrap(), 2.05);

    // THEN transcripts without AVI scores do not get the keys
    assert!(!parsed_transcripts[2].contains_key("avi_phred"));
    assert!(!parsed_transcripts[2].contains_key("avi_raw"));

    // WHEN setting severity predictions on the variant
    set_severity_predictions(&mut variant, &record, &parsed_transcripts);

    // THEN the variant gets the highest AVI scores across transcripts
    assert_eq!(variant.get_f64("avi_phred").unwrap(), 21.3);
    assert_eq!(variant.get_f64("avi_raw").unwrap(), 2.05);
}

#[test]
fn test_parse_avi_scores_missing() {
    // GIVEN a VCF without AVI scores in the CSQ field
    let (record, field_indices) =
        vep_record("FEATURE|AVI_PHRED|AVI_RAW", "ENST00000632593||", &[], "");

    let mut variant = Document::new();

    // WHEN parsing transcripts and setting severity predictions
    let parsed_transcripts = parse_vep_transcripts(&record, &field_indices, &mut variant);
    set_severity_predictions(&mut variant, &record, &parsed_transcripts);

    // THEN no AVI keys are written on the variant
    assert!(!variant.contains_key("avi_phred"));
    assert!(!variant.contains_key("avi_raw"));
}
