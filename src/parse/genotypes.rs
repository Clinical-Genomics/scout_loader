use std::collections::HashMap;
use rust_htslib::bcf::header::HeaderView;

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