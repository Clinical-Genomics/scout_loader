use rust_htslib::bcf::Record;
use mongodb::bson::{doc, Document};
use crate::parse::info::parse_info_string;

pub const ONC_CLNSIG: &[&str] = &[
    "Oncogenic",
    "Likely oncogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
];

/// Remove leading underscores and split grouped ClinVar annotations.
///
/// The input may contain comma or ampersand separated groups, with multiple
/// values inside each group separated by `/`.
/// Spaces are replaced with underscores.
///
/// # Arguments
///
/// * `value` - Raw annotation string.
///
/// # Returns
///
/// A vector of normalized annotation values.
pub fn split_groups(value: &str) -> Vec<String> {
    value
        .replace('&', ",")
        .split(',')
        .flat_map(|group| group.split('/'))
        .map(|item| {
            item.trim_start_matches('_')
                .replace(' ', "_")
        })
        .collect()
}

/// Capitalize the first character of a string.
///
/// Equivalent to Python's `str.capitalize()`.
fn capitalize(value: &str) -> String {
    let mut chars = value.chars();

    match chars.next() {
        Some(first) => {
            first.to_uppercase().collect::<String>() + chars.as_str()
        }
        None => String::new(),
    }
}

/// Collect somatic oncogenicity ClinVar classifications for a variant.
///
/// Parses ClinVar oncogenicity annotations from the VCF INFO fields:
///
/// - `ONC` - oncogenicity classification
/// - `ONCREVSTAT` - review status
/// - `ONCDN` - oncogenicity description/name
/// - `CLNVID` - ClinVar accession
///
/// Entries without a recognized oncogenicity classification are ignored.
///
/// Returns a list of BSON documents representing oncogenicity annotations.
pub fn parse_clnsig_onc(record: &Record) -> Vec<Document> {
    let Some(onc) = parse_info_string(record, b"ONC") else {
        return Vec::new();
    };

    let acc = parse_info_string(record, b"CLNVID")
        .and_then(|value| value.parse::<i64>().ok())
        .unwrap_or(0);

    let onc_sig_groups = split_groups(&onc.to_lowercase());

    let onc_revstat = split_groups(
        &parse_info_string(record, b"ONCREVSTAT")
            .unwrap_or_default()
            .to_lowercase(),
    )
    .join(",");

    let onc_dn_groups = split_groups(
        &parse_info_string(record, b"ONCDN")
            .unwrap_or_default(),
    );

    let mut onc_clnsig_accessions = Vec::new();

    for (i, onc_sig) in onc_sig_groups.iter().enumerate() {
        // Equivalent to:
        // if onc_sig.capitalize() not in ONC_CLNSIG:
        let normalized = capitalize(onc_sig);

        if !ONC_CLNSIG.contains(&normalized.as_str()) {
            continue;
        }

        let dn = onc_dn_groups
            .get(i)
            .map(|value| value.replace('|', ","))
            .unwrap_or_default();

        onc_clnsig_accessions.push(doc! {
            "accession": acc,
            "value": onc_sig,
            "revstat": onc_revstat.clone(),
            "dn": dn,
        });
    }

    onc_clnsig_accessions
}
