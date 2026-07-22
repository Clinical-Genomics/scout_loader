use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use rust_htslib::bcf::Record;
use mongodb::bson::{Bson, Document};
use std::collections::HashSet;
use crate::HashMap;
use crate::parse::info::parse_info_string;
use crate::parse::genes::parse_genes;
use crate::models::gene::GeneAnnotation;


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
        if let Some(csq) = parse_info_string(record, b"CSQ") {
            for transcript_info in csq.split(',') {
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

    transcript.insert(
        "exon",
        entry
            .get("EXON")
            .cloned()
            .unwrap_or_default(),
    );

    transcript.insert(
        "intron",
        entry
            .get("INTRON")
            .cloned()
            .unwrap_or_default(),
    );

    transcript.insert(
        "functional_annotations",
        Bson::Array(
            entry
                .get("Consequence")
                .unwrap_or(&String::new())
                .split('&')
                .map(|x| Bson::String(x.to_string()))
                .collect(),
        ),
    );

    transcript.insert(
        "is_canonical",
        entry
            .get("CANONICAL")
            .map(|x| x == "YES")
            .unwrap_or(false),
    );

    transcript.insert(
        "coding_sequence_name",
        entry
            .get("HGVSc")
            .cloned()
            .unwrap_or_default(),
    );

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



