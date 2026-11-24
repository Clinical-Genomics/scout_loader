use mongodb::bson::Document;
use mongodb::bson::doc;
use rust_htslib::bcf::Record;

use crate::parse::genotypes::get_pathologic_struc;
use crate::parse::info::{
    insert_info_float, insert_info_int, insert_info_string, insert_info_string_decoded,
    parse_info_string, percent_decode,
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
    // Doesn't seem to be used downstream. Stored as e.g. "[0, 1, 2, 6]" for both the
    // old String and the new Integer list format of PathologicStruc (stranger#124).
    if let Some(pathologic_struc) = get_pathologic_struc(record) {
        variant.insert("str_pathologic_struc", format!("{pathologic_struc:?}"));
    }

    insert_info_string(record, variant, b"RU", "str_ru");
    insert_info_string_decoded(record, variant, b"DisplayRU", "str_display_ru");

    insert_info_int(record, variant, b"REF", "str_ref");
    insert_info_int(record, variant, b"RL", "str_len");

    insert_info_string(record, variant, b"STR_STATUS", "str_status");

    insert_info_int(record, variant, b"STR_NORMAL_MAX", "str_normal_max");
    insert_info_int(record, variant, b"STR_PATHOLOGIC_MIN", "str_pathologic_min");

    // Catalog annotations added by stranger may be percent-encoded (stranger#124).
    insert_info_string_decoded(record, variant, b"Disease", "str_disease");
    insert_info_string_decoded(record, variant, b"InheritanceMode", "str_inheritance_mode");

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
    // Catalog annotations added by stranger may be percent-encoded (stranger#124).
    let decoded = |key: &[u8]| parse_info_string(record, key).map(|value| percent_decode(&value));

    let display = decoded(b"SourceDisplay");
    let source_type = decoded(b"Source");
    let source_id = decoded(b"SourceId");

    if display.is_some() || source_type.is_some() || source_id.is_some() {
        let source = doc! {
            "display": display,
            "type": source_type,
            "id": source_id,
        };

        variant.insert("str_source", source);
    }
}
