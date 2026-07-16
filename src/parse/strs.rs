use mongodb::bson::Document;
use rust_htslib::bcf::Record;
use mongodb::bson::doc;

use crate::parse::info::{
    insert_info_float,
    insert_info_int,
    insert_info_string,
    parse_info_string
};

/// Add STR-specific annotations from VCF INFO fields to a MongoDB variant document.
///
/// Only fields present in the VCF record are added. Missing STR annotations
/// are ignored, allowing the same parser to handle STR records with different
/// annotation sets.
pub fn set_str_info(record: &Record, variant: &mut Document) {
    insert_info_float(record, variant, b"SweGenMean", "str_swegen_mean");
    insert_info_float(record, variant, b"SweGenStd", "str_swegen_std");

    insert_info_string(record, variant, b"REPID", "str_repid");
    insert_info_string(record, variant, b"TRID", "str_trid"); // Doesn't seem to be used downstream
    insert_info_string(record, variant, b"STRUC", "str_struc"); // Doesn't seem to be used downstream
    insert_info_string(record, variant, b"MOTIFS", "str_motifs"); // Doesn't seem to be used downstream
    insert_info_string(record, variant, b"PathologicStruc", "str_pathologic_struc"); // Doesn't seem to be used downstream

    insert_info_string(record, variant, b"RU", "str_ru");
    insert_info_string(record, variant, b"DisplayRU", "str_display_ru");

    insert_info_int(record, variant, b"REF", "str_ref");
    insert_info_int(record, variant, b"RL", "str_len");

    insert_info_string(record, variant, b"STR_STATUS", "str_status");

    insert_info_int(record, variant, b"STR_NORMAL_MAX", "str_normal_max");
    insert_info_int(record, variant, b"STR_PATHOLOGIC_MIN", "str_pathologic_min");

    insert_info_string(record, variant, b"Disease", "str_disease");
    insert_info_string(record, variant, b"InheritanceMode", "str_inheritance_mode");

    set_str_source(record, variant);
}

/// Add STR source annotation information to a MongoDB variant document.
///
/// The source information is read from the VCF INFO fields:
/// - `SourceDisplay` → `display`
/// - `Source` → `type`
/// - `SourceId` → `id`
///
/// The `str_source` subdocument is only added if at least one of these
/// fields is present in the VCF record.
///
/// The resulting document structure is:
///
/// ```text
/// "str_source": {
///     "display": "...",
///     "type": "...",
///     "id": "..."
/// }
/// ```
fn set_str_source(record: &Record, variant: &mut Document) {
    let display = parse_info_string(record, b"SourceDisplay");
    let source_type = parse_info_string(record, b"Source");
    let source_id = parse_info_string(record, b"SourceId");

    if display.is_some() || source_type.is_some() || source_id.is_some() {
        let source = doc! {
            "display": display,
            "type": source_type,
            "id": source_id,
        };

        variant.insert("str_source", source);
    }
}