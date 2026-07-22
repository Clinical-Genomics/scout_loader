use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use rust_htslib::bcf::Record;
use mongodb::bson::{Bson, Document};
use std::collections::HashSet;
use crate::HashMap;
use crate::parse::info::parse_info_string;
use crate::parse::genes::parse_genes;


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



/// Parse VEP transcript annotations and add gene information to a variant.
///
/// Extracts CSQ annotations from a VCF record using the VEP header, converts
/// raw VEP transcript entries into parsed transcript documents, collects dbSNP
/// and COSMIC identifiers, and derives gene information from the parsed
/// transcripts.
///
/// Fusion variants are not handled here because they have their own dedicated
/// gene/transcript parser.
///
/// # Arguments
///
/// * `record` - VCF record containing the variant annotations.
/// * `vep_header` - Parsed VEP CSQ header fields used to decode transcript
///   annotations.
/// * `variant` - BSON document representing the parsed variant. Gene, HGNC,
///   dbSNP and COSMIC information is added to this document.
///
/// # Returns
///
/// A vector of parsed transcript documents for downstream processing.
pub fn parse_vep_transcripts(
    record: &Record,
    vep_header: &[String],
    variant: &mut Document,
) -> Vec<Document> {
    let mut parsed_transcripts = Vec::new();

    let mut dbsnp_ids: HashSet<String> = HashSet::new();
    let mut cosmic_ids: HashSet<String> = HashSet::new();

    // Parse VEP CSQ annotations
    if !vep_header.is_empty() {
        if let Some(csq) = parse_info_string(record, b"CSQ") {
            let raw_transcripts = csq.split(',')
                .map(|transcript_info| {
                    vep_header
                        .iter()
                        .zip(transcript_info.split('|'))
                        .map(|(key, value)| {
                            (key.clone(), value.to_string())
                        })
                        .collect::<HashMap<String, String>>()
                });

            for raw_transcript in raw_transcripts {
                if let Some(transcript) = parse_vep_transcript(raw_transcript) {
                    for dbsnp in get_string_array(&transcript, "dbsnp") {
                        dbsnp_ids.insert(dbsnp);
                    }

                    for cosmic in get_string_array(&transcript, "cosmic") {
                        cosmic_ids.insert(cosmic);
                    }

                    parsed_transcripts.push(transcript);
                }
            }
        }
    }

    // COSMIC can also be added independently by VEP/bcftools annotate
    if let Some(cosmic_tag) = parse_info_string(record, b"COSMIC") {
        for cosmic_id in cosmic_tag.split('&') {
            cosmic_ids.insert(cosmic_id.to_string());
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

    // Derive genes from parsed transcripts
    let mut genes = parse_genes(&parsed_transcripts);

    let mut hgnc_ids: HashSet<String> = genes
        .iter()
        .filter_map(|gene| {
            gene.get_str("hgnc_id")
                .ok()
                .map(str::to_string)
        })
        .collect();

    // HGNC IDs annotated by Stranger for STR variants
    if let Some(str_hgnc_id) = parse_info_string(record, b"HGNCId") {
        hgnc_ids.insert(str_hgnc_id.clone());

        if genes.is_empty() {
            let mut gene = Document::new();
            gene.insert("hgnc_id", str_hgnc_id);
            genes.push(gene);
        }
    }

    variant.insert(
        "genes",
        Bson::Array(
            genes
                .into_iter()
                .map(Bson::Document)
                .collect(),
        ),
    );

    variant.insert(
        "hgnc_ids",
        Bson::Array(
            hgnc_ids
                .into_iter()
                .map(Bson::String)
                .collect(),
        ),
    );

    parsed_transcripts
}

/// Parse a single VEP transcript annotation.
///
/// Extracts the basic transcript, gene and consequence information from a
/// single VEP CSQ entry.
///
/// # Arguments
///
/// * `entry` - A single VEP CSQ annotation represented as key-value pairs.
///
/// # Returns
///
/// A parsed transcript document, or `None` if no transcript identifier
/// can be extracted.
pub fn parse_vep_transcript(
    entry: HashMap<String, String>,
) -> Option<Document> {
    let mut transcript = Document::new();

    // Transcript ID
    let transcript_id = entry
        .get("FEATURE")
        .map(|value| value.split(':').next().unwrap_or(""))
        .unwrap_or("");

    if transcript_id.is_empty() {
        return None;
    }

    transcript.insert(
        "transcript_id",
        transcript_id.to_string(),
    );

    // Gene information
    transcript.insert(
        "hgnc_id",
        get_hgnc_id(&entry),
    );

    transcript.insert(
        "hgnc_symbol",
        entry
            .get("SYMBOL")
            .cloned()
            .unwrap_or_default(),
    );

    // Consequences
    let consequences: Vec<Bson> = entry
        .get("CONSEQUENCE")
        .unwrap_or(&String::new())
        .split('&')
        .filter(|value| !value.is_empty())
        .map(|value| Bson::String(value.to_string()))
        .collect();

    transcript.insert(
        "functional_annotations",
        Bson::Array(consequences),
    );

    // Transcript annotations
    transcript.insert(
        "coding_sequence_name",
        entry
            .get("HGVSc")
            .cloned()
            .unwrap_or_default(),
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
        "is_canonical",
        entry
            .get("CANONICAL")
            .map(|value| value == "YES")
            .unwrap_or(false),
    );

    Some(transcript)
}

/// Extract the HGNC identifier from a VEP transcript annotation.
///
/// The HGNC identifier is extracted from the HGNC_ID field. If the field
/// contains a prefix (for example, "HGNC:1234"), only the numeric identifier
/// is retained.
///
/// # Arguments
///
/// * `entry` - A raw VEP transcript annotation represented as key-value pairs.
///
/// # Returns
///
/// A BSON integer containing the HGNC identifier, or BSON null if the field
/// is missing or cannot be parsed.
pub fn get_hgnc_id(entry: &HashMap<String, String>) -> Bson {
    match entry.get("HGNC_ID") {
        Some(value) => {
            value
                .split(':')
                .last()
                .and_then(|id| id.parse::<i32>().ok())
                .map(Bson::Int32)
                .unwrap_or(Bson::Null)
        }
        None => Bson::Null,
    }
}



