use crate::parse::info::parse_info_string;
use rust_htslib::bcf::Record;
use std::collections::HashMap;

/// Checks whether the VCF INFO field contains a pathogenic ClinVar classification.
///
/// Returns `true` if `CLNSIG` contains a pathogenic classification,
/// otherwise `false`.
pub fn is_pathogenic_info(record: &Record) -> bool {
    let Some(clnsig) = parse_info_string(record, b"CLNSIG") else {
        return false;
    };

    clnsig.split(['|', '&', ',']).any(is_pathogenic_value)
}

/// Checks whether any VEP transcript contains a pathogenic ClinVar classification.
///
/// Extracts only the ClinVar significance fields from the VCF `CSQ` INFO
/// annotation, without constructing full transcript documents.
///
/// # Arguments
///
/// * `record` - VCF record containing the `CSQ` INFO annotation.
/// * `field_indices` - Mapping of VEP field names to their positions in `CSQ`.
///
/// # Returns
///
/// `true` if any transcript has a pathogenic, likely pathogenic, or conflicting
/// pathogenic classification; otherwise, `false`.
pub fn is_pathogenic_transcripts(record: &Record, field_indices: &HashMap<String, usize>) -> bool {
    let (Some(&clnsig_idx), Some(&clinvar_clnsig_idx)) = (
        field_indices.get("clnsig"),
        field_indices.get("clinvar_clnsig"),
    ) else {
        return false;
    };

    let Ok(Some(csq)) = record.info(b"CSQ").string() else {
        return false;
    };

    csq.iter().any(|value| {
        String::from_utf8_lossy(value).split(',').any(|transcript| {
            let fields: Vec<&str> = transcript.split('|').collect();

            [clnsig_idx, clinvar_clnsig_idx]
                .iter()
                .filter_map(|&idx| fields.get(idx))
                .any(|value| is_pathogenic_value(value))
        })
    })
}

/// Checks whether a ClinVar significance value indicates pathogenicity.
///
/// Recognizes pathogenic, likely pathogenic, conflicting pathogenic
/// classifications, and their numeric codes.
fn is_pathogenic_value(value: &str) -> bool {
    value
        .to_lowercase()
        .split(['/', '|', '&', ','])
        .any(|term| {
            matches!(
                term.trim().trim_start_matches('_'),
                "pathogenic"
                    | "likely_pathogenic"
                    | "conflicting_classifications_of_pathogenicity"
                    | "conflicting_interpretations_of_pathogenicity"
                    | "conflicting_interpretations"
                    | "4"
                    | "5"
                    | "8"
            )
        })
}
