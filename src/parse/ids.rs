use crate::models::variant::VariantIds;
use crate::utils::hash::generate_md5_key;

/// Constructs the identifiers associated with a variant.
///
/// # Arguments
///
/// * `chrom` - Variant chromosome.
/// * `pos` - Variant position.
/// * `reference` - Reference allele.
/// * `alternative` - Alternative allele.
/// * `case_id` - Unique case identifier.
/// * `variant_type` - Variant type (`clinical` or `research`).
///
/// # Returns
///
/// A `VariantIds` object containing the generated identifiers.
pub fn parse_ids(
    chrom: &str,
    pos: &u64,
    reference: &str,
    alternative: &str,
    case_id: &str,
    variant_type: &str,
) -> VariantIds {

    VariantIds {
        simple_id: parse_simple_id(chrom, pos, reference, alternative),
        variant_id: parse_variant_id(chrom, pos, reference, alternative, variant_type),
        display_name: parse_display_name(chrom, pos, reference, alternative, variant_type),
        document_id: parse_document_id(
            chrom,
            pos,
            reference,
            alternative,
            variant_type,
            case_id,
        ),
    }
}

/// Generates the simple identifier for a variant.
///
/// The simple identifier is a human-readable representation of the variant
/// location and alleles. It is **not** guaranteed to be unique.
pub fn parse_simple_id(
    chrom: &str,
    pos: &u64,
    reference: &str,
    alternative: &str,
) -> String {
    format!("{chrom}_{pos}_{reference}_{alternative}")
}

/// Generates the variant identifier for a variant.
///
/// The variant identifier uniquely identifies a variant within a specific
/// analysis type (`clinical` or `research`). It is generated as an MD5 hash
/// and is therefore not intended to be human-readable.
pub fn parse_variant_id(
    chrom: &str,
    pos: &u64,
    reference: &str,
    alternative: &str,
    variant_type: &str,
) -> String {
    generate_md5_key(&[
        chrom,
        &pos.to_string(),
        reference,
        alternative,
        variant_type,
    ])
}

/// Generates the display name for a variant.
///
/// The display name is a human-readable identifier used to display the
/// variant in Scout.
pub fn parse_display_name(
    chrom: &str,
    pos: &u64,
    reference: &str,
    alternative: &str,
    variant_type: &str,
) -> String {
    format!(
        "{chrom}_{pos}_{reference}_{alternative}_{variant_type}"
    )
}

/// Generates the unique document identifier for a variant.
///
/// The document identifier uniquely identifies a variant document in the
/// database. It is generated as an MD5 hash using the variant coordinates,
/// alleles, analysis type, and case identifier.
pub fn parse_document_id(
    chrom: &str,
    pos: &u64,
    reference: &str,
    alternative: &str,
    variant_type: &str,
    case_id: &str,
) -> String {
    generate_md5_key(&[
        chrom,
        &pos.to_string(),
        reference,
        alternative,
        variant_type,
        case_id,
    ])
}