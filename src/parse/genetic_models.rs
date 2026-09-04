use rust_htslib::bcf::Record;

/// Parses genetic models from the `GeneticModels` INFO field of a VCF record.
///
/// The `GeneticModels` field may contain model annotations for one or more
/// cases, separated by commas. Each case entry is expected to follow the
/// format:
///
/// ```text
/// <case_id>:<model1>|<model2>|...
/// ```
///
/// The case identifier is ignored when parsing. This keeps the parsing logic
/// independent of the case being loaded, which is important when loading
/// cloned cases with a different case ID.
///
/// # Arguments
///
/// * `record` - VCF record containing the `GeneticModels` INFO annotation.
///
/// # Returns
///
/// A vector containing all genetic models found in the `GeneticModels` field.
/// Returns an empty vector if the INFO field is missing or contains no valid
/// model entries.
pub fn parse_genetic_models(record: &Record) -> Vec<String> {
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

    models_info
        .split(',')
        .filter_map(|family_info| family_info.split_once(':'))
        .flat_map(|(_, models)| models.split('|'))
        .map(str::to_string)
        .collect()
}
