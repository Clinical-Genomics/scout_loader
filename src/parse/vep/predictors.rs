use crate::parse::vep::transcripts::VepFields;
use mongodb::bson::{Bson, Document};

/// Extract a prediction from VEP transcript annotations.
///
/// The prediction fields are typically formatted as
/// `prediction(score)`, for example `deleterious(0.01)`.
/// Returns only the prediction label. If none of the provided fields
/// are available or contain a value, `"unknown"` is returned.
pub fn get_prediction(entry: &VepFields<'_>, fields: &[&str]) -> String {
    for field in fields {
        if let Some(value) = entry.get(field).filter(|value| !value.is_empty()) {
            return value.split('(').next().unwrap_or("unknown").to_string();
        }
    }

    "unknown".to_string()
}

/// Extracts SpliceAI delta scores and positions from VEP CSQ fields.
/// Stores the maximum delta score with its corresponding position,
/// along with a summary of all splice predictions.
pub fn parse_transcripts_spliceai(transcript: &mut Document, entry: &VepFields<'_>) {
    let spliceai_fields = [
        ("AG", "SPLICEAI_PRED_DS_AG", "SPLICEAI_PRED_DP_AG"),
        ("AL", "SPLICEAI_PRED_DS_AL", "SPLICEAI_PRED_DP_AL"),
        ("DG", "SPLICEAI_PRED_DS_DG", "SPLICEAI_PRED_DP_DG"),
        ("DL", "SPLICEAI_PRED_DS_DL", "SPLICEAI_PRED_DP_DL"),
    ];

    let mut max_score: Option<f64> = None;
    let mut max_position: Option<i32> = None;
    let mut predictions = Vec::new();
    let mut has_prediction = false;

    for (suffix, score_field, position_field) in spliceai_fields {
        let score = entry
            .get(score_field)
            .filter(|value| !value.is_empty())
            .and_then(|value| value.parse::<f64>().ok());

        let position = entry
            .get(position_field)
            .filter(|value| !value.is_empty())
            .and_then(|value| value.parse::<i32>().ok());

        if let Some(score) = score {
            transcript.insert(format!("spliceai_ds_{suffix}"), Bson::Double(score));

            if max_score.is_none_or(|current| score > current) {
                max_score = Some(score);
                max_position = position;
            }
        }

        if let Some(position) = position {
            transcript.insert(format!("spliceai_dp_{suffix}"), Bson::Int32(position));
        }

        if score.is_some_and(|value| value != 0.0) || position.is_some_and(|value| value != 0) {
            has_prediction = true;
        }

        let suffix = suffix.to_lowercase();

        predictions.push(format!(
            "spliceai_ds_{suffix} {} spliceai_dp_{suffix} {}",
            score
                .filter(|value| *value != 0.0)
                .map_or_else(|| "-".to_string(), |value| value.to_string()),
            position
                .filter(|value| *value != 0)
                .map_or_else(|| "-".to_string(), |value| value.to_string()),
        ));
    }

    if let Some(score) = max_score {
        transcript.insert("spliceai_score", Bson::Double(score));

        if let Some(position) = max_position {
            transcript.insert("spliceai_position", Bson::Int32(position));
        }
    }

    transcript.insert(
        "spliceai_prediction",
        if has_prediction {
            Bson::Array(predictions.into_iter().map(Bson::String).collect())
        } else {
            Bson::Null
        },
    );
}
