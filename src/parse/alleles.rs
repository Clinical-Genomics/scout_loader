use rust_htslib::bcf::Record;
use crate::models::variant::VariantCategory;

/// Parses the reference and alternative alleles from a VCF record.
///
/// # Arguments
///
/// * `record` - A VCF record containing allele information.
///
/// # Returns
///
/// A tuple containing:
/// * reference allele (REF)
/// * alternative allele (ALT)
///
/// # Panics
///
/// Panics if the record does not contain a reference or alternative allele.
pub fn parse_alleles(record: &Record, category: VariantCategory,) -> (String, String) {
    let alleles = record.alleles();

    let reference = String::from_utf8_lossy(alleles[0]).to_string();

    let alternative = if alleles.len() > 1 {
        if alleles.len() > 2 {
            panic!("VCF records must be split and normalized.");
        }
        String::from_utf8_lossy(alleles[1]).to_string()
    } else if category == VariantCategory::Str {
        ".".to_string()
    } else {
        panic!("VCF record has no ALT allele.");
    };

    (reference, alternative)
}