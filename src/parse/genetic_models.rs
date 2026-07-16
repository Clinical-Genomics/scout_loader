use rust_htslib::bcf::Record;

/// Parses the genetic models for a specific case from a VCF record.
///
/// The `GeneticModels` INFO field contains information for one or more cases,
/// separated by commas. Each case entry is expected to follow the format:
///
/// ```text
/// <case_id>:<model1>|<model2>|...
/// ```
///
/// The function extracts the genetic models associated with the provided
/// case ID.
///
/// # Arguments
///
/// * `record` - VCF record containing the `GeneticModels` INFO annotation.
/// * `case_id` - Case identifier used to select the corresponding models.
///
/// # Returns
///
/// A vector containing the genetic models for the specified case. Returns an
/// empty vector if no models are found or the INFO field is missing.
pub fn parse_genetic_models(
    record: &Record,
    case_id: &str,
) -> Vec<String> {
    let models_info = record
        .info(b"GeneticModels")
        .string()
        .ok()
        .flatten()
        .and_then(|values| {
            values
                .first()
                .map(|value| String::from_utf8_lossy(value).to_string())
        });

    let Some(models_info) = models_info else {
        return Vec::new();
    };

    for family_info in models_info.split(',') {
        let split_info: Vec<&str> = family_info.split(':').collect();

        if split_info.len() < 2 {
            continue;
        }

        if split_info[0] == case_id {
            return split_info[1]
                .split('|')
                .map(|model| model.to_string())
                .collect();
        }
    }

    Vec::new()
}