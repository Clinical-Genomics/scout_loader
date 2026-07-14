use rust_htslib::bcf::Record;

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
pub fn parse_rank_scores(
    record: &Record,
    case_id: &str,
) -> (i32, f64) {
    let rank_score_entry = record
        .info(b"RankScore")
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values.first().map(|value| {
                String::from_utf8_lossy(value).to_string()
            })
        });

    let norm_rank_score_entry = record
        .info(b"RankScoreNormalized")
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values.first().map(|value| {
                String::from_utf8_lossy(value).to_string()
            })
        });

    let rank_score = parse_score_entry(
        rank_score_entry.as_deref(),
        case_id,
    )
    .and_then(|score| score.parse::<i32>().ok())
    .unwrap_or(0);

    let norm_rank_score = parse_score_entry(
        norm_rank_score_entry.as_deref(),
        case_id,
    )
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
pub fn parse_score_entry<'a>(
    score_entry: Option<&'a str>,
    case_id: &str,
) -> Option<&'a str> {
    let Some(score_entry) = score_entry else {
        return None;
    };

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




