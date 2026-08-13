use crate::parse::info::{parse_info_float, parse_info_string};
use mongodb::bson::{Bson, Document, doc};
use rust_htslib::bcf::Record;

pub const MIMVIR_SCORE_KEY: &str = "MivmirScore";
pub const MIMVIR_SCORE_DESC: &str = "MivmirExplanation";
pub const GICAM_SCORE_KEY: &str = "GicamScore";

/// Parses the rank score annotations for a variant from a VCF record.
///
/// This function extracts the `RankScore` and `RankScoreNormalized` INFO
/// annotations from the VCF record and retrieves the values corresponding
/// to the provided case ID.
///
/// Missing or invalid scores are replaced with default values:
/// `0` for the integer rank score and `0.0` for the normalized rank score.
///
/// # Arguments
///
/// * `record` - VCF record containing the rank score annotations.
/// * `case_id` - Case identifier used to select the corresponding scores.
///
/// # Returns
///
/// A tuple containing:
///
/// * `rank_score` - The integer rank score for the case (`i32`).
/// * `norm_rank_score` - The normalized rank score for the case (`f64`).
pub fn parse_rank_scores(record: &Record, case_id: &str) -> (i32, f64) {
    let rank_score_entry = record
        .info(b"RankScore")
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values
                .first()
                .map(|value| String::from_utf8_lossy(value).to_string())
        });

    let norm_rank_score_entry = record
        .info(b"RankScoreNormalized")
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values
                .first()
                .map(|value| String::from_utf8_lossy(value).to_string())
        });

    let rank_score = parse_score_entry(rank_score_entry.as_deref(), case_id)
        .and_then(|score| score.parse::<i32>().ok())
        .unwrap_or(0);

    let norm_rank_score = parse_score_entry(norm_rank_score_entry.as_deref(), case_id)
        .and_then(|score| score.parse::<f64>().ok())
        .unwrap_or(0.0);

    (rank_score, norm_rank_score)
}

/// Extracts a score value for a specific case from a raw VCF annotation.
///
/// The annotation contains scores for multiple cases, separated by commas.
/// Each entry follows the format `<case_id>:<score>`.
///
/// # Arguments
///
/// * `score_entry` - Optional raw score annotation string.
/// * `case_id` - Case identifier used to select the corresponding score.
///
/// # Returns
///
/// The score value as a string slice if the case ID is found, otherwise `None`.
pub fn parse_score_entry<'a>(score_entry: Option<&'a str>, case_id: &str) -> Option<&'a str> {
    let score_entry = score_entry?;

    for family_info in score_entry.split(',') {
        let mut split_info = family_info.split(':');

        let Some(entry_case_id) = split_info.next() else {
            continue;
        };

        let Some(score) = split_info.next() else {
            continue;
        };

        if entry_case_id == case_id {
            return Some(score);
        }
    }

    None
}

/// Parse the `RankResult` INFO field into category/score entries.
///
/// The categories are taken from the `RankResult` header and the scores
/// from the corresponding pipe-separated values in the INFO field.
pub fn parse_rank_result(record: &Record, rank_results_header: &[String]) -> Option<Vec<Bson>> {
    let rank_result = parse_info_string(record, b"RankResult")?;

    let results = rank_results_header
        .iter()
        .zip(rank_result.split('|'))
        .filter_map(|(category, score)| {
            score.parse::<i32>().ok().map(|score| {
                Bson::Document(doc! {
                    "category": category,
                    "score": score,
                })
            })
        })
        .collect();

    Some(results)
}

/// Parse additional rank scores from a VCF record.
///
/// Parses the Mivmir and Gicam scores. The Mivmir score can optionally
/// include an explanation from the `MivmirExplanation` INFO field.
///
/// Returns `None` if no additional rank scores are present.
pub fn parse_rank_score_other(record: &Record) -> Option<Document> {
    let mut rank_scores = Document::new();

    if let Some(value) = parse_info_float(record, MIMVIR_SCORE_KEY.as_bytes()) {
        let mut mivmir = Document::new();
        mivmir.insert("value", Bson::Double(value));

        if let Some(desc) = parse_rank_score_description(record, MIMVIR_SCORE_DESC.as_bytes()) {
            mivmir.insert("desc", desc);
        }

        rank_scores.insert("Mivmir", Bson::Document(mivmir));
    }

    if let Some(value) = parse_info_float(record, GICAM_SCORE_KEY.as_bytes()) {
        let mut gicam = Document::new();
        gicam.insert("value", Bson::Double(value));

        rank_scores.insert("Gicam", Bson::Document(gicam));
    }

    if rank_scores.is_empty() {
        None
    } else {
        Some(rank_scores)
    }
}

/// Parse the explanation associated with a rank score.
///
/// The VCF INFO value is expected to contain comma-separated `key=value`
/// pairs, optionally enclosed in square brackets. Values are converted
/// to floating-point numbers and stored in a BSON document.
///
/// Returns `None` if the explanation is missing or cannot be parsed.
fn parse_rank_score_description(record: &Record, key: &[u8]) -> Option<Bson> {
    let raw = parse_info_string(record, key)?;

    let raw = raw
        .trim_start_matches('[')
        .trim_end_matches(']')
        .trim_end_matches(',');

    if raw.is_empty() {
        return None;
    }

    let mut description = Document::new();

    for item in raw.split(',') {
        let (key, value) = item.split_once('=')?;
        let value = value.parse::<f64>().ok()?;

        description.insert(key, Bson::Double(value));
    }

    if description.is_empty() {
        None
    } else {
        Some(Bson::Document(description))
    }
}
