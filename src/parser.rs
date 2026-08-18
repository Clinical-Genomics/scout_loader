use crate::models::build::GenomeBuild;
use crate::models::case::CaseConfig;
use crate::models::variant::{VariantCategory, VariantType};
use crate::parse::cytobands::set_cytobands;
use crate::parse::vcf::process_vcf;
use mongodb::bson::Document;
use std::collections::{HashMap, HashSet};

/// Parses and processes all clinical VCFs provided for a case.
///
/// Shared information such as the genome build and cytobands is prepared once
/// and reused for every VCF. Sample mappings are created separately for each
/// VCF because sample indices may differ between VCF files.
///
/// The variant type is derived from the corresponding VCF key in the case
/// configuration. Research VCFs are intentionally excluded.
pub fn parse(
    config: &CaseConfig,
    gene_to_panels: &HashMap<i32, HashSet<String>>,
    hgncid_to_gene: &HashMap<i32, Document>,
) -> Result<(), Box<dyn std::error::Error>> {
    let genome_build = GenomeBuild::from_str(&config.human_genome_build)
        .map_err(|_| format!("Invalid genome build: {}", config.human_genome_build))?;

    let cytobands = set_cytobands(genome_build.cytoband_path())
        .map_err(|error| format!("Could not load cytobands: {error}"))?;

    // Load only clinical variants for the time being
    let vcfs = [
        (&config.vcf_snv, VariantCategory::Snv),
        (&config.vcf_cancer, VariantCategory::Cancer),
        (&config.vcf_sv, VariantCategory::Sv),
        (&config.vcf_cancer_sv, VariantCategory::CancerSv),
        (&config.vcf_fusion, VariantCategory::Fusion),
        (&config.vcf_mei, VariantCategory::Mei),
        (&config.vcf_str, VariantCategory::Str),
    ];

    let variant_type =
        VariantType::from_str("clinical").map_err(|_| "Invalid variant type: clinical")?;

    for (vcf, category) in vcfs {
        if let Some(vcf) = vcf {
            process_vcf(
                vcf.to_str().ok_or("Invalid VCF path")?,
                category,
                variant_type,
                &config.family,
                &cytobands,
                &config.samples,
                gene_to_panels,
                hgncid_to_gene,
            );
        }
    }

    Ok(())
}
