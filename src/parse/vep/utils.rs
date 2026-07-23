use crate::HashMap;

/// Return the highest float value from a string with numbers possibly
/// separated by `&`.
///
/// Invalid values are ignored. Returns `None` if the input is empty or
/// no valid float values are found.
pub fn get_highest_float_score_in_string(value: &str) -> Option<f64> {
    value
        .split('&')
        .filter_map(|part| part.trim().parse::<f64>().ok())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
}

/// Extract a sequence annotation from a VEP transcript entry.
///
/// The VEP fields `HGVSC` and `HGVSP` are formatted as
/// `transcript:sequence`. This function returns only the sequence part.
/// Returns `None` if the field is missing or does not contain `:`.
pub fn get_sequence_aux(entry: &HashMap<String, String>, name: &str) -> Option<String> {
    let sequence_entry = entry
        .get(name)?
        .split(':')
        .collect::<Vec<&str>>();

    if sequence_entry.len() > 1 {
        Some(sequence_entry.last()?.to_string())
    } else {
        None
    }
}
