use std::collections::HashMap;
use rust_htslib::bcf::header::HeaderView;
use mongodb::bson::{doc, Document};
use rust_htslib::bcf::Record;

use crate::models::sample::SampleInfo;

/// Validate that provided sample IDs match the samples present in the VCF header.
///
/// Checks that each requested sample exists at the expected VCF position.
pub fn validate_sample_mapping(
    header: &HeaderView,
    sample_mapping: &HashMap<String, SampleInfo>,
) -> Result<(), String> {
    let vcf_samples = header.samples();

    for (sample_id, sample_info) in sample_mapping {
        let index = sample_info.vcf_index;

        let vcf_sample = vcf_samples
            .get(index)
            .ok_or_else(|| {
                format!(
                    "Sample '{}' has VCF index {}, but the VCF contains only {} samples",
                    sample_id,
                    index,
                    vcf_samples.len()
                )
            })?;

        let vcf_sample = String::from_utf8_lossy(vcf_sample);

        if vcf_sample != *sample_id {
            return Err(format!(
                "Sample mismatch: CLI sample '{}' is assigned to VCF index {}, but VCF contains '{}'",
                sample_id,
                index,
                vcf_sample
            ));
        }
    }

    Ok(())
}


/// Parse genotype information for all selected samples.
///
/// Returns a MongoDB document for each sample.
pub fn parse_genotypes(
    record: &Record,
    sample_mapping: &HashMap<String, SampleInfo>,
) -> Vec<Document> {
    sample_mapping
        .iter()
        .map(|(sample_id, sample_info)| {
            parse_genotype(record, sample_id, sample_info)
        })
        .collect()
}

/// Parse genotype information for a single sample.
///
/// Extracts genotype-related FORMAT fields for the sample identified
/// by `sample_info.vcf_index` and returns a MongoDB document.
fn parse_genotype(
    record: &Record,
    sample_id: &str,
    sample_info: &SampleInfo,
) -> Document {
    let pos = sample_info.vcf_index;

    let mut genotype = doc! {
        "sample_id": sample_id,
        "display_name": &sample_info.display_name,
    };

    // GT
    // if let Some(gt) = parse_genotype_call(record, pos) {
    //     genotype.insert("genotype_call", gt);
    // }

    // STR-specific fields

    // SV-specific fields

    // MEI-specific fields

    // Derived fields

    genotype
}