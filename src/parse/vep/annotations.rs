use mongodb::bson::{Bson, Document};
use crate::HashMap;
use std::collections::HashSet;
use crate::models::consequence::SO_TERMS;


/// Get strand information from a VEP transcript entry.
pub fn get_strand(entry: &HashMap<String, String>) -> Option<String> {
    match entry.get("STRAND").map(String::as_str) {
        Some("1") => Some("+".to_string()),
        Some("-1") => Some("-".to_string()),
        _ => None,
    }
}

/// Get regional annotations from Sequence Ontology consequence terms.
///
/// Maps functional annotations (SO consequence terms) to their broader
/// functional regions using the predefined Sequence Ontology term mapping.
///
/// Unknown annotations are ignored.
///
/// # Arguments
///
/// * `functional_annotations` - A list of Sequence Ontology consequence terms
///   extracted from a VEP transcript annotation.
///
/// # Returns
///
/// A vector containing the regional annotations corresponding to the provided
/// functional annotations.
pub fn get_regional_annotation(functional_annotations: &[String]) -> Vec<String> {
    functional_annotations
        .iter()
        .filter_map(|annotation| {
            SO_TERMS
                .get(annotation.as_str())
                .map(|term| term.region.to_string())
        })
        .collect()
}


/// Parse MANE transcript annotations from a VEP transcript entry.
///
/// Extracts MANE Select and MANE Plus Clinical transcript identifiers from
/// VEP v103/MANE v0.92 annotations. Falls back to the older `MANE` field
/// used by previous VEP versions.
pub fn parse_mane_annotations(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if entry.contains_key("MANE_SELECT") {
        if let Some(mane_select) = entry.get("MANE_SELECT").filter(|v| !v.is_empty()) {
            transcript.insert(
                "mane_select_transcript",
                Bson::String(mane_select.clone()),
            );
        }

        if let Some(mane_plus_clinical) = entry
            .get("MANE_PLUS_CLINICAL")
            .filter(|v| !v.is_empty())
        {
            transcript.insert(
                "mane_plus_clinical_transcript",
                Bson::String(mane_plus_clinical.clone()),
            );
        }
    } else if entry.contains_key("MANE") {
        if let Some(mane) = entry.get("MANE").filter(|v| !v.is_empty()) {
            transcript.insert(
                "mane_select_transcript",
                Bson::String(mane.clone()),
            );
        }
    }
}


/// Parse genomic superdups fractional match values from a VEP transcript entry.
pub fn parse_superdups_fracmatch(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if let Some(superdups_fractmatch) = entry
        .get("GENOMIC_SUPERDUPS_FRAC_MATCH")
        .filter(|value| !value.is_empty())
    {
        let values: Vec<Bson> = superdups_fractmatch
            .split('&')
            .filter_map(|fractmatch| {
                fractmatch.parse::<f64>().ok().map(Bson::Double)
            })
            .collect();

        if !values.is_empty() {
            transcript.insert(
                "superdups_fracmatch",
                Bson::Array(values),
            );
        }
    }
}

/// Parse ClinVar annotations from a VEP transcript entry.
///
/// Extracts ClinVar variation identifiers, clinical significance, review
/// status, and clinical significance terms when available.
pub fn parse_clinvar_annotations(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let clinvar_id = entry
        .get("CLINVAR_CLNVID")
        .or_else(|| entry.get("CLINVAR"));

    if let Some(clinvar_id) = clinvar_id {
        transcript.insert(
            "clinvar_clnvid",
            Bson::String(clinvar_id.clone()),
        );

        if let Some(clnsig) = entry.get("CLINVAR_CLNSIG") {
            transcript.insert(
                "clinvar_clnsig",
                Bson::String(clnsig.to_lowercase()),
            );
        }

        if let Some(revstat) = entry.get("CLINVAR_CLNREVSTAT") {
            transcript.insert(
                "clinvar_revstat",
                Bson::String(revstat.to_lowercase()),
            );
        }
    }

    let clnsig = entry
        .get("CLIN_SIG")
        .or_else(|| entry.get("ClinVar_CLNSIG"));

    if let Some(clnsig) = clnsig.filter(|value| !value.is_empty()) {
        transcript.insert(
            "clnsig",
            Bson::Array(
                clnsig
                    .split('&')
                    .map(|value| Bson::String(value.to_string()))
                    .collect(),
            ),
        );
    }
}

/// Parse dbSNP identifiers from a VEP transcript entry.
///
/// Supports different VEP versions where dbSNP identifiers can be stored
/// in `EXISTING_VARIATION`, `RS_DBSNP150`, or `RS_DBSNP`. Only identifiers
/// starting with `rs` are retained.
pub fn parse_dbsnp(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let dbsnp_keys = [
        "EXISTING_VARIATION",
        "RS_DBSNP150",
        "RS_DBSNP",
    ];

    let mut dbsnp_ids: HashSet<String> = HashSet::new();

    for key in dbsnp_keys {
        if let Some(variant_ids) = entry.get(key) {
            for variant_id in variant_ids.split('&') {
                if variant_id.starts_with("rs") {
                    dbsnp_ids.insert(variant_id.to_string());
                }
            }
        }
    }

    transcript.insert(
        "dbsnp",
        Bson::Array(
            dbsnp_ids
                .into_iter()
                .map(Bson::String)
                .collect(),
        ),
    );
}

/// Parse COSMIC identifiers from a VEP transcript entry.
///
/// Extracts COSMIC identifiers from `EXISTING_VARIATION` (COSM/COSV
/// prefixes) and from the dedicated `COSMIC` field when available.
pub fn parse_cosmic(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let mut cosmic_ids: Vec<Bson> = Vec::new();

    if let Some(variant_ids) = entry.get("EXISTING_VARIATION") {
        for variant_id in variant_ids.split('&') {
            if variant_id.starts_with("COSM") || variant_id.starts_with("COSV") {
                cosmic_ids.push(Bson::String(variant_id.to_string()));
            }
        }
    }

    if let Some(cosmic_ids_entry) = entry.get("COSMIC") {
        for cosmic_id in cosmic_ids_entry.split('&') {
            cosmic_ids.push(Bson::String(cosmic_id.to_string()));
        }
    }

    transcript.insert(
        "cosmic",
        Bson::Array(cosmic_ids),
    );
}