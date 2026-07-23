use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use rust_htslib::bcf::Record;
use mongodb::bson::{Bson, Document};
use std::collections::HashSet;
use crate::HashMap;
use crate::parse::genes::parse_genes;
use crate::models::gene::GeneAnnotation;
use crate::models::consequence::SO_TERMS;


/// Extract the VEP CSQ annotation fields from the VCF header.
///
/// The VEP annotations are stored in the CSQ INFO field. The field names are
/// defined after `Format:` in the CSQ description and are separated by `|`.
/// Returned field names are converted to uppercase.
///
/// Returns an empty vector if no CSQ INFO header is found.
pub fn parse_vep_header(header: &HeaderView) -> Vec<String> {
    for record in header.header_records() {
        if let HeaderRecord::Info { values, .. } = record {
            if let Some(id) = values.get("ID") {
                if id == "CSQ" {
                    if let Some(description) = values.get("Description") {
                        if let Some(format) = description.split("Format:").nth(1) {
                            return format
                                .trim()
                                .trim_matches('"')
                                .trim_end_matches('>')
                                .split('|')
                                .map(|field| field.to_uppercase())
                                .collect();
                        }
                    }
                }
            }
        }
    }

    Vec::new()
}

/// Extract string values from a BSON document field.
///
/// The field can either be stored as a BSON array of strings or as a single
/// string value. Any other BSON type, missing field, or non-string array
/// elements are ignored.
///
/// Args:
///     document: BSON document containing the field.
///     key: Name of the field to extract.
///
/// Returns:
///     A vector containing the extracted string values.
fn get_string_array(document: &Document, key: &str) -> Vec<String> {
    match document.get(key) {
        Some(Bson::Array(values)) => values
            .iter()
            .filter_map(|value| value.as_str().map(str::to_string))
            .collect(),

        Some(Bson::String(value)) => vec![value.clone()],

        _ => Vec::new(),
    }
}



/// Parse VEP CSQ annotations from a VCF record.
///
/// Extracts transcript annotations from the CSQ INFO field using the VEP
/// header, collects dbSNP and COSMIC identifiers, and builds gene-level
/// information separately from transcript-level information.
///
/// Gene identifiers are not stored in transcripts to avoid duplication.
/// They are passed separately to `parse_genes`.
///
/// Returns the parsed transcripts and gene annotations.
pub fn parse_vep_transcripts(
    record: &Record,
    vep_header: &[String],
    variant: &mut Document,
) -> (Vec<Document>, Vec<GeneAnnotation>) {
    let mut parsed_transcripts = Vec::new();
    let mut gene_annotations = Vec::new();

    let mut dbsnp_ids = HashSet::new();
    let mut cosmic_ids = HashSet::new();

    if !vep_header.is_empty() {
        if let Ok(Some(csq)) = record.info(b"CSQ").string() {
            let csq_string = csq
                .iter()
                .map(|value| String::from_utf8_lossy(value))
                .collect::<Vec<_>>()
                .join(",");
            for transcript_info in csq_string.split(',') {
                println!("tx info: {}", transcript_info);
                let raw_transcript: HashMap<String, String> = vep_header
                    .iter()
                    .zip(transcript_info.split('|'))
                    .map(|(key, value)| {
                        (key.clone(), value.to_string())
                    })
                    .collect();

                gene_annotations.push(GeneAnnotation {
                    hgnc_id: get_hgnc_id(&raw_transcript),
                    hgnc_symbol: raw_transcript
                        .get("SYMBOL")
                        .cloned(),
                });

                if let Some(transcript) = parse_vep_transcript(raw_transcript) {
                    if let Ok(values) = transcript.get_array("dbsnp") {
                        for value in values {
                            if let Bson::String(id) = value {
                                dbsnp_ids.insert(id.clone());
                            }
                        }
                    }

                    if let Ok(values) = transcript.get_array("cosmic") {
                        for value in values {
                            if let Bson::String(id) = value {
                                cosmic_ids.insert(id.clone());
                            }
                        }
                    }

                    parsed_transcripts.push(transcript);
                }
            }
        }
    }

    if !dbsnp_ids.is_empty() && !variant.contains_key("dbsnp_id") {
        variant.insert(
            "dbsnp_id",
            Bson::String(
                dbsnp_ids.into_iter().collect::<Vec<_>>().join(";"),
            ),
        );
    }

    if !cosmic_ids.is_empty() {
        variant.insert(
            "cosmic_ids",
            Bson::Array(
                cosmic_ids
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }

    (parsed_transcripts, gene_annotations)
}

/// Parse a single VEP transcript annotation.
///
/// Extracts transcript-specific information from a single CSQ entry.
/// Gene-level information such as HGNC ID and gene symbol is handled
/// separately by the gene parser.
///
/// Returns None if no transcript ID is available.
pub fn parse_vep_transcript(
    entry: HashMap<String, String>,
) -> Option<Document> {
    let transcript_id = entry
        .get("FEATURE")
        .map(|id| id.split(':').next().unwrap_or(""))
        .unwrap_or("");

    if transcript_id.is_empty() {
        return None;
    }

    let mut transcript = Document::new();

    transcript.insert(
        "transcript_id",
        transcript_id.to_string(),
    );

    transcript.insert(
        "protein_id",
        entry
            .get("ENSP")
            .cloned()
            .unwrap_or_default(),
    );

    let polyphen = get_prediction(&entry, &["POLYPHEN"]);

    transcript.insert(
        "polyphen_prediction",
        Bson::String(polyphen),
    );

    let sift = get_prediction(&entry, &["SIFT", "SIFT_PRED"]);
    transcript.insert(
        "sift_prediction",
        Bson::String(sift),
    );

    if let Some(value) = entry.get("REVEL_RANKSCORE").filter(|v| !v.is_empty()) {
        if let Some(rankscore) = get_highest_float_score_in_string(value) {
            transcript.insert(
                "revel_rankscore",
                Bson::Double(rankscore),
            );
        }
    }

    if let Some(value) = entry.get("REVEL_SCORE").filter(|v| !v.is_empty()) {
        if let Some(score) = get_highest_float_score_in_string(value) {
            transcript.insert(
                "revel_raw_score",
                Bson::Double(score),
            );
        }
    }

    parse_transcripts_spliceai(&mut transcript, &entry);


    transcript.insert(
        "swiss_prot",
        Bson::String(
            entry
                .get("SWISSPROT")
                .filter(|value| !value.is_empty())
                .cloned()
                .unwrap_or_else(|| "unknown".to_string()),
        ),
    );


    for (source, target) in [
        ("GERP++_RS", "gerp"),
        ("PHASTCONS100WAY_VERTEBRATE", "phast"),
        ("PHYLOP100WAY_VERTEBRATE", "phylop"),
    ] {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            transcript.insert(
                target,
                Bson::String(value.clone()),
            );
        }
    }

    parse_domains(&mut transcript, &entry);

    if let Some(coding_sequence) = get_sequence_aux(&entry, "HGVSC") {
        transcript.insert(
            "coding_sequence_name",
            Bson::String(coding_sequence),
        );
    }

    if let Some(coding_sequence) = get_sequence_aux(&entry, "HGVSP") {
        transcript.insert(
            "protein_sequence_name",
            Bson::String(coding_sequence),
        );
    }

    transcript.insert(
        "biotype",
        Bson::String(
            entry.get("BIOTYPE")
                .cloned()
                .unwrap_or_default(),
        ),
    );

    for (source, target) in [("EXON", "exon"),("INTRON", "intron")] {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            transcript.insert(
                target,
                Bson::String(value.clone()),
            );
        }
    }

    if let Some(strand) = get_strand(&entry) {
        transcript.insert(
            "strand",
            Bson::String(strand),
        );
    }

    let functional_annotations: Vec<String> = entry
        .get("CONSEQUENCE")
        .unwrap_or(&String::new())
        .split('&')
        .map(|x| x.to_string())
        .collect();

    let functional_annotations: Vec<String> = entry
        .get("CONSEQUENCE")
        .filter(|value| !value.is_empty())
        .map(|value| value.split('&').map(String::from).collect())
        .unwrap_or_default();

    if !functional_annotations.is_empty() {
        transcript.insert(
            "functional_annotations",
            Bson::Array(
                functional_annotations
                    .iter()
                    .map(|annotation| Bson::String(annotation.clone()))
                    .collect(),
            ),
        );

        let region_annotations = get_regional_annotation(&functional_annotations);

        transcript.insert(
            "region_annotations",
            Bson::Array(
                region_annotations
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }



    /*

    // Canonical transcript
    transcript.insert(
        "is_canonical",
        Bson::Boolean(
            entry.get("CANONICAL")
                .map(|x| x == "YES")
                .unwrap_or(false),
        ),
    );


    */

    Some(transcript)
}

/// Extract the HGNC identifier from a VEP annotation.
///
/// The HGNC_ID field may contain a prefix (for example, "HGNC:1234").
/// Only the numeric identifier is retained.
///
/// Returns `None` if the HGNC identifier is missing or invalid.
pub fn get_hgnc_id(entry: &HashMap<String, String>) -> Option<String> {
    entry
        .get("HGNC_ID")
        .and_then(|value| value.split(':').last())
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
}

/// Extract a prediction from VEP transcript annotations.
///
/// The prediction fields are typically formatted as
/// `prediction(score)`, for example `deleterious(0.01)`.
/// Returns only the prediction label. If none of the provided fields
/// are available or contain a value, `"unknown"` is returned.
fn get_prediction(entry: &HashMap<String, String>, fields: &[&str]) -> String {
    for field in fields {
        if let Some(value) = entry.get(*field).filter(|value| !value.is_empty()) {
            return value
                .split('(')
                .next()
                .unwrap_or("unknown")
                .to_string();
        }
    }

    "unknown".to_string()
}

/// Return the highest float value from a string with numbers possibly
/// separated by `&`.
///
/// Invalid values are ignored. Returns `None` if the input is empty or
/// no valid float values are found.
fn get_highest_float_score_in_string(value: &str) -> Option<f64> {
    value
        .split('&')
        .filter_map(|part| part.trim().parse::<f64>().ok())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
}


/// Parse SpliceAI annotations from a VEP transcript entry.
///
/// Extracts SpliceAI delta scores and positions from VEP CSQ fields.
/// The maximum delta score is stored together with its corresponding
/// position. Also stores a summary of all splice predictions.
fn parse_transcripts_spliceai(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let spliceai_positions = [
        ("SPLICEAI_PRED_DP_AG", "spliceai_dp_ag"),
        ("SPLICEAI_PRED_DP_AL", "spliceai_dp_al"),
        ("SPLICEAI_PRED_DP_DG", "spliceai_dp_dg"),
        ("SPLICEAI_PRED_DP_DL", "spliceai_dp_dl"),
    ];

    let spliceai_delta_scores = [
        ("SPLICEAI_PRED_DS_AG", "spliceai_ds_ag"),
        ("SPLICEAI_PRED_DS_AL", "spliceai_ds_al"),
        ("SPLICEAI_PRED_DS_DG", "spliceai_ds_dg"),
        ("SPLICEAI_PRED_DS_DL", "spliceai_ds_dl"),
    ];

    for (source, target) in spliceai_positions {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            if let Ok(position) = value.parse::<i32>() {
                transcript.insert(
                    target,
                    Bson::Int32(position),
                );
            }
        }
    }

    for (source, target) in spliceai_delta_scores {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            if let Ok(score) = value.parse::<f64>() {
                transcript.insert(
                    target,
                    Bson::Double(score),
                );
            }
        }
    }

    let spliceai_pairs = [
        ("spliceai_ds_ag", "spliceai_dp_ag"),
        ("spliceai_ds_al", "spliceai_dp_al"),
        ("spliceai_ds_dg", "spliceai_dp_dg"),
        ("spliceai_ds_dl", "spliceai_dp_dl"),
    ];

    let mut max_score: Option<f64> = None;
    let mut max_position: Option<i32> = None;
    let mut predictions = Vec::new();

    for (score_key, position_key) in spliceai_pairs {
        let score = transcript.get_f64(score_key).ok();
        let position = transcript.get_i32(position_key).ok();

        if let Some(score) = score {
            if max_score.map_or(true, |current| score > current) {
                max_score = Some(score);
                max_position = position;
            }
        }

        predictions.push(format!(
            "{} {} {} {}",
            score_key,
            score.map(|x| x.to_string()).unwrap_or("-".to_string()),
            position_key,
            position.map(|x| x.to_string()).unwrap_or("-".to_string()),
        ));
    }

    if let Some(score) = max_score {
        transcript.insert(
            "spliceai_score",
            Bson::Double(score),
        );

        if let Some(position) = max_position {
            transcript.insert(
                "spliceai_position",
                Bson::Int32(position),
            );
        }

        transcript.insert(
            "spliceai_prediction",
            Bson::Array(
                predictions
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }
}

/// Parse protein domain annotations from a VEP transcript entry.
///
/// Extracts supported protein domains from the VEP `DOMAINS` field.
fn parse_domains(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if let Some(domains) = entry.get("DOMAINS").filter(|v| !v.is_empty()) {
        for annotation in domains.split('&') {
            let parts: Vec<&str> = annotation.split(':').collect();

            if parts.len() < 2 {
                continue;
            }

            let domain_name = parts[0];
            let domain_id = parts[1];

            match domain_name {
                "Pfam_domain" => {
                    transcript.insert(
                        "pfam_domain",
                        Bson::String(domain_id.to_string()),
                    );
                }
                "PROSITE_profiles" => {
                    transcript.insert(
                        "prosite_profile",
                        Bson::String(domain_id.to_string()),
                    );
                }
                "SMART_domains" => {
                    transcript.insert(
                        "smart_domain",
                        Bson::String(domain_id.to_string()),
                    );
                }
                "hmmpanther" => {
                    transcript.insert(
                        "panther_domain",
                        Bson::String(domain_id.to_string()),
                    );
                }
                _ => {}
            }
        }
    }
}

/// Extract a sequence annotation from a VEP transcript entry.
///
/// The VEP fields `HGVSC` and `HGVSP` are formatted as
/// `transcript:sequence`. This function returns only the sequence part.
/// Returns `None` if the field is missing or does not contain `:`.
fn get_sequence_aux(entry: &HashMap<String, String>, name: &str) -> Option<String> {
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

/// Get strand information from a VEP transcript entry.
fn get_strand(entry: &HashMap<String, String>) -> Option<String> {
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
fn get_regional_annotation(functional_annotations: &[String]) -> Vec<String> {
    functional_annotations
        .iter()
        .filter_map(|annotation| {
            SO_TERMS
                .get(annotation.as_str())
                .map(|term| term.region.to_string())
        })
        .collect()
}