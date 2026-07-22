use rust_htslib::bcf::Record;
use mongodb::bson::{Document};

use crate::parse::info::{parse_info_string};

/// Parse Mitomap-associated diseases from the VCF INFO field.
///
/// If the `MitomapAssociatedDiseases` annotation is present and not equal to
/// `"."`, stores it in the variant document after replacing underscores with
/// spaces.
pub fn set_mitomap_associated_diseases(
    record: &Record,
    variant: &mut Document,
) {
    if let Some(value) = parse_info_string(record, b"MitomapAssociatedDiseases") {
        if value != "." {
            variant.insert(
                "mitomap_associated_diseases",
                value.replace('_', " "),
            );
        }
    }
}


/// Parse the HmtVar variant identifier from the VCF INFO field.
///
/// If the `HmtVar` annotation is present and not equal to `"."`, stores the
/// variant identifier in the variant document.
pub fn set_hmtvar(
    record: &Record,
    variant: &mut Document,
) {
    if let Some(value) = parse_info_string(record, b"HmtVar") {
        if value != "." {
            if let Ok(id) = value.parse::<i32>() {
                variant.insert("hmtvar_variant_id", id);
            }
        }
    }
}