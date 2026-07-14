use rust_htslib::bcf::Record;

/// Parses the rank score and normalized rank score for a variant.
///
/// This function extracts the `RankScore` and `RankScoreNormalized` INFO
/// annotations from the VCF record and uses `parse_rank_score` to retrieve
/// the values corresponding to the provided case ID.
///
/// Missing or invalid scores are returned as `None`.
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
/// * `rank_score` - The primary rank score for the case, if available.
/// * `norm_rank_score` - The normalized rank score for the case, if available.
pub fn parse_rank_scores(
    record: &Record,
    case_id: &str,
) -> (Option<i32>, Option<i32>) {
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

    let rank_score = parse_rank_score(
        rank_score_entry.as_deref(),
        case_id,
    );

    let norm_rank_score = parse_rank_score(
        norm_rank_score_entry.as_deref(),
        case_id,
    );

    (rank_score, norm_rank_score)
}

/// Parses a rank score for a specific case from a raw rank score entry.
///
/// The rank score entry contains scores for multiple cases, separated by
/// commas. Each entry is expected to follow the format `<case_id>:<score>`.
/// The function extracts the score corresponding to the provided case ID.
///
/// # Arguments
///
/// * `rank_score_entry` - Optional raw rank score annotation string.
/// * `case_id` - Case identifier used to select the corresponding score.
///
/// # Returns
///
/// The parsed score as `Some(i32)` if the case ID is found and the score can
/// be converted successfully, otherwise `None`.
pub fn parse_rank_score(
    rank_score_entry: Option<&str>,
    case_id: &str,
) -> Option<i32> {
    let Some(rank_score_entry) = rank_score_entry else {
        return None;
    };

    for family_info in rank_score_entry.split(',') {
        let split_info: Vec<&str> = family_info.split(':').collect();

        if split_info.len() < 2 {
            continue;
        }

        if split_info[0] == case_id {
            return split_info[1].parse::<i32>().ok();
        }
    }

    None
}




