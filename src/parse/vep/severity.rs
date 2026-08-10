use mongodb::bson::Document;
use rust_htslib::bcf::Record;

use crate::parse::info::parse_info_float;

/// Get the highest score for a field across parsed transcripts.
///
/// Returns `None` if no transcript contains a valid numeric value.
fn get_highest_transcript_score(transcripts: &[Document], key: &str) -> Option<f64> {
    transcripts
        .iter()
        .filter_map(|tx| tx.get_f64(key).ok())
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
}

/// Parse the CADD score from a VCF record and parsed transcripts.
///
/// Checks the `CADD` and `CADD_PHRED` INFO fields first. If neither is
/// available, returns the highest CADD score found in the transcripts.
/// Returns `0.0` if no CADD score is available.
fn parse_cadd(record: &Record, transcripts: &[Document]) -> f64 {
    for key in [b"CADD".as_slice(), b"CADD_PHRED".as_slice()] {
        if let Some(value) = parse_info_float(record, key) {
            return value;
        }
    }

    transcripts
        .iter()
        .filter_map(|tx| tx.get_f64("cadd").ok())
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(0.0)
}

/// Set severity prediction scores on a parsed variant.
///
/// Adds the CADD and SPIDEX scores from the VCF INFO fields. When parsed
/// transcripts are available, also adds the highest REVEL rank score and
/// highest REVEL raw score found across all transcripts.
///
/// The following fields are added when values are available:
/// - `cadd_score`
/// - `spidex`
/// - `revel_score`
/// - `revel`
pub fn set_severity_predictions(
    variant: &mut Document,
    record: &Record,
    parsed_transcripts: &[Document],
) {
    variant.insert("cadd_score", parse_cadd(record, parsed_transcripts));

    if let Some(value) = parse_info_float(record, b"SPIDEX") {
        variant.insert("spidex", value);
    }

    if !parsed_transcripts.is_empty() {
        if let Some(value) = get_highest_transcript_score(parsed_transcripts, "revel_rankscore") {
            variant.insert("revel_score", value);
        }

        if let Some(value) = get_highest_transcript_score(parsed_transcripts, "revel_raw_score") {
            variant.insert("revel", value);
        }
    }
}
