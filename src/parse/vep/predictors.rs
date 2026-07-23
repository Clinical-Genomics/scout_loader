use mongodb::bson::{Bson, Document};
use crate::HashMap;

/// Extract a prediction from VEP transcript annotations.
///
/// The prediction fields are typically formatted as
/// `prediction(score)`, for example `deleterious(0.01)`.
/// Returns only the prediction label. If none of the provided fields
/// are available or contain a value, `"unknown"` is returned.
pub fn get_prediction(entry: &HashMap<String, String>, fields: &[&str]) -> String {
    for field in fields {
        if let Some(value) = entry.get(*field).filter(|value| !value.is_empty()) {
            return value
                .split('(')
                .next()
                .unwrap_or("unknown")
                .to_string();
        }
    }

    "unknown".to_string()
}

/// Parse SpliceAI annotations from a VEP transcript entry.
///
/// Extracts SpliceAI delta scores and positions from VEP CSQ fields.
/// The maximum delta score is stored together with its corresponding
/// position. Also stores a summary of all splice predictions.
pub fn parse_transcripts_spliceai(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let spliceai_positions = [
        ("SPLICEAI_PRED_DP_AG", "spliceai_dp_ag"),
        ("SPLICEAI_PRED_DP_AL", "spliceai_dp_al"),
        ("SPLICEAI_PRED_DP_DG", "spliceai_dp_dg"),
        ("SPLICEAI_PRED_DP_DL", "spliceai_dp_dl"),
    ];

    let spliceai_delta_scores = [
        ("SPLICEAI_PRED_DS_AG", "spliceai_ds_ag"),
        ("SPLICEAI_PRED_DS_AL", "spliceai_ds_al"),
        ("SPLICEAI_PRED_DS_DG", "spliceai_ds_dg"),
        ("SPLICEAI_PRED_DS_DL", "spliceai_ds_dl"),
    ];

    for (source, target) in spliceai_positions {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            if let Ok(position) = value.parse::<i32>() {
                transcript.insert(
                    target,
                    Bson::Int32(position),
                );
            }
        }
    }

    for (source, target) in spliceai_delta_scores {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            if let Ok(score) = value.parse::<f64>() {
                transcript.insert(
                    target,
                    Bson::Double(score),
                );
            }
        }
    }

    let spliceai_pairs = [
        ("spliceai_ds_ag", "spliceai_dp_ag"),
        ("spliceai_ds_al", "spliceai_dp_al"),
        ("spliceai_ds_dg", "spliceai_dp_dg"),
        ("spliceai_ds_dl", "spliceai_dp_dl"),
    ];

    let mut max_score: Option<f64> = None;
    let mut max_position: Option<i32> = None;
    let mut predictions = Vec::new();

    for (score_key, position_key) in spliceai_pairs {
        let score = transcript.get_f64(score_key).ok();
        let position = transcript.get_i32(position_key).ok();

        if let Some(score) = score {
            if max_score.map_or(true, |current| score > current) {
                max_score = Some(score);
                max_position = position;
            }
        }

        predictions.push(format!(
            "{} {} {} {}",
            score_key,
            score.map(|x| x.to_string()).unwrap_or("-".to_string()),
            position_key,
            position.map(|x| x.to_string()).unwrap_or("-".to_string()),
        ));
    }

    if let Some(score) = max_score {
        transcript.insert(
            "spliceai_score",
            Bson::Double(score),
        );

        if let Some(position) = max_position {
            transcript.insert(
                "spliceai_position",
                Bson::Int32(position),
            );
        }

        transcript.insert(
            "spliceai_prediction",
            Bson::Array(
                predictions
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }
}