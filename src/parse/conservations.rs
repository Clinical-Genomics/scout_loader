use mongodb::bson::{Bson, Document};
use rust_htslib::bcf::Record;

const CONSERVATION: &[(&str, f64)] = &[("gerp", 2.0), ("phast", 0.8), ("phylop", 2.5)];

/// Return the minimum score required for a conservation field to be
/// classified as conserved.
fn conserved_min(field_key: &str) -> f64 {
    CONSERVATION
        .iter()
        .find(|(key, _)| *key == field_key)
        .map(|(_, threshold)| *threshold)
        .unwrap_or(0.0)
}

/// Format a prediction with Scout's rounded score representation.
fn format_conservation(score: f64, field_key: &str) -> Bson {
    let prediction = if score >= conserved_min(field_key) {
        "Conserved"
    } else {
        "NotConserved"
    };
    let rounded = format!("{score:.2}");
    // Python's round(score, 2) retains at least one decimal for floats.
    let rounded = rounded.strip_suffix('0').unwrap_or(&rounded);

    Bson::String(format!("{prediction} ({rounded})"))
}

/// Parse and add conservation predictions to a variant.
///
/// Conservation scores are read from variant INFO fields when available,
/// otherwise from the first parsed VEP transcript. Predictions are stored
/// separately for GERP, PhastCons, and PhyloP.
pub fn parse_conservations(
    record: &Record,
    parsed_transcripts: &[Document],
    variant: &mut Document,
) {
    let conservation_keys = [
        ("gerp", b"dbNSFP_GERP___RS".as_slice()),
        ("phast", b"dbNSFP_phastCons100way_vertebrate".as_slice()),
        ("phylop", b"dbNSFP_phyloP100way_vertebrate".as_slice()),
    ];

    for (field_key, info_key) in conservation_keys {
        let mut conservation = parse_conservation_info(record, info_key, field_key);

        if conservation.is_empty()
            && let Some(transcript) = parsed_transcripts.first()
        {
            conservation = parse_conservation_csq(transcript, field_key);
        }

        variant.insert(
            format!("{field_key}_conservation"),
            Bson::Array(conservation),
        );
    }
}

/// Parse conservation scores from a VCF INFO field.
///
/// Converts each score into a `Conserved` or `NotConserved` annotation
/// based on the minimum conservation threshold for the specified field.
///
/// Returns an empty vector if the INFO field is missing or contains no scores.
fn parse_conservation_info(record: &Record, info_key: &[u8], field_key: &str) -> Vec<Bson> {
    let Some(scores) = record.info(info_key).float().ok().flatten() else {
        return Vec::new();
    };

    scores
        .iter()
        .map(|&score| format_conservation(f64::from(score), field_key))
        .collect()
}

/// Parse conservation scores from a parsed VEP transcript.
///
/// The transcript field may contain multiple scores separated by `&`.
/// Each score is converted into a `Conserved` or `NotConserved` annotation
/// based on the minimum conservation threshold for the specified field.
///
/// Invalid or missing scores are ignored.
fn parse_conservation_csq(transcript: &Document, field_key: &str) -> Vec<Bson> {
    let Some(Bson::String(value)) = transcript.get(field_key) else {
        return Vec::new();
    };

    value
        .split('&')
        .filter_map(|score| score.parse::<f64>().ok())
        .map(|score| format_conservation(score, field_key))
        .collect()
}
