use mongodb::bson::Document;
use rust_htslib::bcf::Record;

use crate::parse::info::{
    insert_info_float,
    insert_info_int,
    insert_info_string,
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
    insert_info_string(record, variant, b"TRID", "str_trid");
    insert_info_string(record, variant, b"STRUC", "str_struc");
    insert_info_string(record, variant, b"MOTIFS", "str_motifs");
    insert_info_string(record, variant, b"PathologicStruc", "str_pathologic_struc");

    insert_info_string(record, variant, b"RU", "str_ru");
    insert_info_string(record, variant, b"DisplayRU", "str_display_ru");

    insert_info_int(record, variant, b"REF", "str_ref");
    insert_info_int(record, variant, b"RL", "str_len");

    insert_info_string(record, variant, b"STR_STATUS", "str_status");

    insert_info_int(record, variant, b"STR_NORMAL_MAX", "str_normal_max");
    insert_info_int(record, variant, b"STR_PATHOLOGIC_MIN", "str_pathologic_min");

    insert_info_string(record, variant, b"Disease", "str_disease");
    insert_info_string(record, variant, b"InheritanceMode", "str_inheritance_mode");
}