use std::collections::HashMap;
use rust_htslib::bcf::header::HeaderView;
use mongodb::bson::{doc, Document};
use rust_htslib::bcf::Record;
use rust_htslib::bcf::record::GenotypeAllele;

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


/// Parse genotype calls for selected samples.
pub fn parse_genotypes(
    record: &Record,
    sample_mapping: &HashMap<String, SampleInfo>,
) -> Vec<Document> {
    let mut genotypes = Vec::new();

    for (sample_id, sample_info) in sample_mapping {
        let pos = sample_info.vcf_index;

        genotypes.push(parse_genotype(
            record,
            sample_id,
            &sample_info.display_name,
            pos,
        ));
    }

    genotypes
}

/// Parse genotype information for a single sample.
///
/// Extracts the GT field for the sample identified by its VCF index
/// and returns a MongoDB document containing the genotype call.
fn parse_genotype(
    record: &Record,
    sample_id: &str,
    display_name: &str,
    pos: usize,
) -> Document {
    let mut gt_call = doc! {
        "sample_id": sample_id,
        "display_name": display_name,
    };

    let genotypes = record.genotypes().expect("Could not read genotypes");
    let genotype = genotypes.get(pos);

    let allele_1 = genotype_allele_to_string(genotype.get(0));
    let allele_2 = genotype_allele_to_string(genotype.get(1));

    let phase_sep = match genotype.get(1) {
        Some(GenotypeAllele::Phased(_)) | Some(GenotypeAllele::PhasedMissing) => "|",
        _ => "/",
    };

    gt_call.insert(
        "genotype_call",
        format!("{}{}{}", allele_1, phase_sep, allele_2),
    );

    // STR-specific fields

    // SV-specific fields

    // MEI-specific fields

    // Derived fields

    gt_call
}

/// Converts a genotype allele to its string representation.
///
/// Maps allele indexes to their VCF genotype representation.
/// Missing or unsupported alleles are represented as ".".
fn genotype_allele_to_string(allele: Option<&GenotypeAllele>) -> &'static str {
    match allele {
        Some(GenotypeAllele::Unphased(0)) | Some(GenotypeAllele::Phased(0)) => "0",
        Some(GenotypeAllele::Unphased(1)) | Some(GenotypeAllele::Phased(1)) => "1",
        Some(GenotypeAllele::UnphasedMissing)
        | Some(GenotypeAllele::PhasedMissing)
        | None => ".",
        _ => ".",
    }
}