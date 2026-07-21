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
/// is converted into an array of database names, and `fusion_genes` is
/// populated with the pair of fusion partner genes.
pub fn set_fusion_info(
    record: &Record,
    variant: &mut Document,
) {
    
    fn replace_nan(value: Option<String>, nan_value: &str) -> String {
        match value {
            Some(v) if v != nan_value => v,
            _ => String::new(),
        }
    }

    fn set_found_db(found_db_info: Option<String>) -> Option<Vec<String>> {
        found_db_info.and_then(|value| {
            if value.is_empty() || value == "[]" {
                None
            } else {
                Some(value.split(',').map(|s| s.to_string()).collect())
            }
        })
    }

    let gene_a = parse_info_string(record, b"GENEA").unwrap_or_default();
    let gene_b = parse_info_string(record, b"GENEB").unwrap_or_default();

    variant.insert("gene_a", gene_a.clone());
    variant.insert("gene_b", gene_b.clone());

    if let Some(value) = parse_info_int(record, b"TOOL_HITS") {
        variant.insert("tool_hits", value);
    }

    if let Some(value) = set_found_db(parse_info_string(record, b"FOUND_DB")) {
        variant.insert(
            "found_db",
            Bson::Array(value.into_iter().map(Bson::String).collect()),
        );
    }

    if let Some(value) = parse_info_float(record, b"SCORE") {
        variant.insert("fusion_score", value);
    }

    if let Some(value) = parse_info_int(record, b"HGNC_ID_A") {
        variant.insert("hgnc_id_a", value);
    }

    if let Some(value) = parse_info_int(record, b"HGNC_ID_B") {
        variant.insert("hgnc_id_b", value);
    }

    variant.insert(
        "orientation",
        replace_nan(parse_info_string(record, b"ORIENTATION"), "nan,nan"),
    );

    variant.insert(
        "frame_status",
        replace_nan(parse_info_string(record, b"FRAME_STATUS"), "nan"),
    );

    variant.insert(
        "transcript_id_a",
        replace_nan(parse_info_string(record, b"TRANSCRIPT_ID_A"), "nan"),
    );

    variant.insert(
        "transcript_id_b",
        replace_nan(parse_info_string(record, b"TRANSCRIPT_ID_B"), "nan"),
    );

    if let Some(value) = replace_nan(parse_info_string(record, b"EXON_NUMBER_A"), "nan")
        .parse::<i32>()
        .ok()
    {
        variant.insert("exon_number_a", value.to_string());
    }

    if let Some(value) = replace_nan(parse_info_string(record, b"EXON_NUMBER_B"), "nan")
        .parse::<i32>()
        .ok()
    {
        variant.insert("exon_number_b", value.to_string());
    }

    variant.insert(
        "fusion_genes",
        Bson::Array(vec![Bson::String(gene_a), Bson::String(gene_b)]),
    );

}