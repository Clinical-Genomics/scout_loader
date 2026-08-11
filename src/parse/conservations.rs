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

/// Parse conservation predictions from variant INFO fields or VEP transcripts.
///
/// Variant-level conservation annotations are preferred. If an annotation is
/// not available in INFO, the first parsed transcript is checked instead.
///
/// Returns conservation predictions for GERP, PhastCons, and PhyloP.
pub fn parse_conservations(record: &Record, parsed_transcripts: &[Document]) -> Document {
    let conservation_keys = [
        ("gerp", b"dbNSFP_GERP___RS".as_slice()),
        ("phast", b"dbNSFP_phastCons100way_vertebrate".as_slice()),
        ("phylop", b"dbNSFP_phyloP100way_vertebrate".as_slice()),
    ];

    let mut conservations = Document::new();

    for (field_key, info_key) in conservation_keys {
        let result = parse_conservation_info(record, info_key, field_key)
            .or_else(|| {
                parsed_transcripts
                    .first()
                    .map(|transcript| Bson::Array(parse_conservation_csq(transcript, field_key)))
            })
            .unwrap_or(Bson::Array(Vec::new()));

        conservations.insert(field_key, result);
    }

    conservations
}

/// Parse conservation scores from a VCF INFO field.
///
/// Classifies each score as `Conserved` or `NotConserved` based on the
/// minimum conservation threshold for the specified field.
fn parse_conservation_info(record: &Record, info_key: &[u8], field_key: &str) -> Option<Bson> {
    let scores = record.info(info_key).float().ok().flatten()?;

    let values = scores
        .iter()
        .map(|&score| {
            let score = f64::from(score);
            let label = if score >= conserved_min(field_key) {
                format!("Conserved ({score:.2})")
            } else {
                format!("NotConserved ({score:.2})")
            };

            Bson::String(label)
        })
        .collect();

    Some(Bson::Array(values))
}

/// Parse a conservation score from a parsed VEP transcript.
///
/// The transcript field may contain multiple scores separated by `&`.
/// Each score is converted into a `Conserved` or `NotConserved` annotation
/// based on the conservation threshold for the given field.
///
/// Invalid or missing scores are ignored.
fn parse_conservation_csq(transcript: &Document, field_key: &str) -> Vec<Bson> {
    let Some(Bson::String(value)) = transcript.get(field_key) else {
        return Vec::new();
    };

    value
        .split('&')
        .filter_map(|score| score.parse::<f64>().ok())
        .map(|score| {
            let label = if score >= conserved_min(field_key) {
                format!("Conserved ({:.2})", score)
            } else {
                format!("NotConserved ({:.2})", score)
            };

            Bson::String(label)
        })
        .collect()
}
