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