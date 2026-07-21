use mongodb::bson::{Bson, Document};
use rust_htslib::bcf::Record;

use crate::parse::info::{
    parse_info_float,
    parse_info_int,
    parse_info_string,
};

/// Add fusion-specific information from VCF INFO fields.
///
/// Parses fusion-related annotations and inserts them into the variant
/// document. Missing or placeholder values (e.g. `"nan"` or `"nan,nan"`)
/// are replaced with empty strings where appropriate. The `FOUND_DB` field
/// is converted into an array of database names. Gene and transcript
/// information for fusion partners is also added.
pub fn set_fusion_info(record: &Record, variant: &mut Document) {
    fn replace_nan(value: Option<String>, nan_value: &str) -> String {
        match value {
            Some(value) if value != nan_value => value,
            _ => String::new(),
        }
    }

    fn parse_found_db(value: Option<String>) -> Option<Vec<Bson>> {
        value.and_then(|value| {
            if value.is_empty() || value == "[]" {
                None
            } else {
                Some(
                    value
                        .split(',')
                        .map(|entry| Bson::String(entry.to_string()))
                        .collect(),
                )
            }
        })
    }

    if let Some(value) = parse_info_int(record, b"TOOL_HITS") {
        variant.insert("tool_hits", value);
    }

    if let Some(value) = parse_found_db(parse_info_string(record, b"FOUND_DB")) {
        variant.insert("found_db", Bson::Array(value));
    }

    if let Some(value) = parse_info_float(record, b"SCORE") {
        variant.insert("fusion_score", value);
    }

    variant.insert(
        "orientation",
        replace_nan(parse_info_string(record, b"ORIENTATION"), "nan,nan"),
    );

    variant.insert(
        "frame_status",
        replace_nan(parse_info_string(record, b"FRAME_STATUS"), "nan"),
    );

    set_fusion_genes(record, variant);
}

/// Add gene and transcript information for fusion variants.
///
/// Parses fusion partner annotations from VCF INFO fields and populates the
/// variant document with gene information, transcripts, HGNC IDs, and HGNC
/// symbols. Gene entries are added when either a gene symbol or an HGNC ID is
/// available. Transcript and exon information is retained even when HGNC IDs
/// are missing.
fn set_fusion_genes(record: &Record, variant: &mut Document) {
    let mut genes = Vec::new();
    let mut hgnc_ids = Vec::new();
    let mut hgnc_symbols = Vec::new();

    for suffix in ["A", "B"] {
        let gene = parse_info_string(
            record,
            format!("GENE{suffix}").as_bytes(),
        )
        .unwrap_or_default();

        let hgnc_id = parse_info_int(
            record,
            format!("HGNC_ID_{suffix}").as_bytes(),
        );

        let transcript_id = parse_info_string(
            record,
            format!("TRANSCRIPT_ID_{suffix}").as_bytes(),
        );

        let exon_number = parse_info_string(
            record,
            format!("EXON_NUMBER_{suffix}").as_bytes(),
        );

        // Keep the gene entry if we have at least a symbol or an identifier.
        if gene.is_empty() && hgnc_id.is_none() {
            continue;
        }

        if let Some(hgnc_id) = hgnc_id {
            hgnc_ids.push(Bson::Int32(hgnc_id));
        }

        if !gene.is_empty() {
            hgnc_symbols.push(Bson::String(gene.clone()));
        }

        let mut gene_doc = Document::new();

        if !gene.is_empty() {
            gene_doc.insert("hgnc_symbol", gene.clone());
        }

        if let Some(hgnc_id) = hgnc_id {
            gene_doc.insert("hgnc_id", hgnc_id);
        }

        let mut transcripts = Vec::new();

        if let Some(transcript_id) = transcript_id {
            let mut transcript = Document::new();

            transcript.insert("transcript_id", transcript_id);

            if let Some(hgnc_id) = hgnc_id {
                transcript.insert("hgnc_id", hgnc_id);
            }

            if !gene.is_empty() {
                transcript.insert("hgnc_symbol", gene.clone());
            }

            if let Some(exon_number) = exon_number {
                transcript.insert("exon", exon_number);
            }

            transcripts.push(Bson::Document(transcript));
        }

        gene_doc.insert(
            "transcripts",
            Bson::Array(transcripts),
        );

        genes.push(Bson::Document(gene_doc));
    }

    variant.insert("genes", Bson::Array(genes));
    variant.insert("hgnc_ids", Bson::Array(hgnc_ids));
    variant.insert("hgnc_symbols", Bson::Array(hgnc_symbols));
}