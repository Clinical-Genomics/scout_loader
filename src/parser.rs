use crate::loader::Loader;
use crate::models::build::GenomeBuild;
use crate::models::case::CaseConfig;
use crate::models::variant::VariantAnnotations;
use crate::models::variant::{VariantCategory, VariantType};
use crate::parse::cytobands::set_cytobands;
use crate::parse::vcf::process_vcf;
use crate::updater;
use std::str::FromStr;

/// Parses and processes all clinical VCFs provided for a case.
///
/// Shared information such as the genome build and cytobands is prepared once
/// and reused for every VCF. Sample mappings are created separately for each
/// VCF because sample indices may differ between VCF files.
///
/// The variant type is derived from the corresponding VCF key in the case
/// configuration. Research VCFs are intentionally excluded.
///
/// After each VCF is processed, the number of newly inserted variants is
/// checked. If variants were inserted, their `variant_rank` values are
/// updated for the corresponding variant category and type.
///
/// Returns the total number of variants inserted across all VCFs for the case.
pub async fn parse(
    config: &CaseConfig,
    mut annotations: VariantAnnotations<'_>,
    loader: &Loader,
) -> Result<usize, Box<dyn std::error::Error>> {
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

    let mut inserted_variants = 0;
    let mut total_inserted_variants = 0;

    for (vcf, category) in vcfs {
        if let Some(vcf) = vcf {
            annotations.managed_variant_ids = loader
                .get_managed_variant_ids(&category.to_string(), &config.human_genome_build)
                .await?;

            inserted_variants = process_vcf(
                vcf.to_str().ok_or("Invalid VCF path")?,
                category,
                variant_type,
                config,
                &cytobands,
                &annotations,
                loader,
            )
            .await?;
        }

        println!("{category:?}: {inserted_variants} variants added");
        total_inserted_variants += inserted_variants;

        if inserted_variants > 0 {
            updater::update_variant_rank(loader, &config.family, variant_type, category).await?;
        }

        inserted_variants = 0;
    }

    Ok(total_inserted_variants)
}
