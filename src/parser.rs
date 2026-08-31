use crate::loader::Loader;
use crate::models::build::GenomeBuild;
use crate::models::case::CaseConfig;
use crate::models::variant::VariantAnnotations;
use crate::models::variant::{VariantCategory, VariantType};
use crate::parse::cytobands::set_cytobands;
use crate::parse::vcf::process_vcf;
use crate::updater;
use std::str::FromStr;

/// Parse and process the selected VCFs provided for a case.
///
/// By default, all available clinical VCFs are processed. When `research` is
/// true, all available research VCFs are processed instead.
///
/// If `categories` is provided, only the requested categories are processed.
/// The `research` flag determines whether the clinical or research VCF is
/// selected for each category.
///
/// Shared information such as the genome build and cytobands is prepared once
/// and reused for every VCF. Sample mappings are created separately for each
/// VCF because sample indices may differ between VCF files.
///
/// After each VCF is processed, the number of newly inserted variants is
/// checked. If variants were inserted, their `variant_rank` values are
/// updated for the corresponding variant category and type.
///
/// Returns the total number of variants inserted across all selected VCFs.
pub async fn parse(
    config: &CaseConfig,
    mut annotations: VariantAnnotations<'_>,
    loader: &Loader,
    categories: Option<&[String]>,
    research: bool,
) -> Result<usize, Box<dyn std::error::Error>> {
    let genome_build = GenomeBuild::from_str(&config.human_genome_build)
        .map_err(|_| format!("Invalid genome build: {}", config.human_genome_build))?;

    let cytobands = set_cytobands(genome_build.cytoband_path())
        .map_err(|error| format!("Could not load cytobands: {error}"))?;

    let variant_type = VariantType::from_str(if research { "research" } else { "clinical" })
        .map_err(|_| "Invalid variant type")?;

    let vcfs = [
        (
            &config.vcf_snv,
            &config.vcf_snv_research,
            VariantCategory::Snv,
            "snv",
        ),
        (
            &config.vcf_cancer,
            &config.vcf_cancer_research,
            VariantCategory::Cancer,
            "cancer",
        ),
        (
            &config.vcf_sv,
            &config.vcf_sv_research,
            VariantCategory::Sv,
            "sv",
        ),
        (
            &config.vcf_cancer_sv,
            &config.vcf_cancer_sv_research,
            VariantCategory::CancerSv,
            "cancer_sv",
        ),
        (
            &config.vcf_fusion,
            &config.vcf_fusion_research,
            VariantCategory::Fusion,
            "fusion",
        ),
        (
            &config.vcf_mei,
            &config.vcf_mei_research,
            VariantCategory::Mei,
            "mei",
        ),
        (
            &config.vcf_str,
            &config.vcf_str_research,
            VariantCategory::Str,
            "str",
        ),
    ];

    let mut total_inserted_variants = 0;

    for (clinical_vcf, research_vcf, category, category_name) in vcfs {
        // If categories were explicitly requested, skip everything else.
        if let Some(categories) = categories
            && !categories
                .iter()
                .any(|requested| requested == category_name)
        {
            continue;
        }

        // Select either the clinical or research VCF.
        let vcf = if research { research_vcf } else { clinical_vcf };

        let Some(vcf) = vcf else {
            continue;
        };

        annotations.managed_variant_ids = loader
            .get_managed_variant_ids(&category.to_string(), &config.human_genome_build)
            .await?;

        let inserted_variants = process_vcf(
            vcf.to_str().ok_or("Invalid VCF path")?,
            category,
            variant_type,
            config,
            &cytobands,
            &annotations,
            loader,
        )
        .await?;

        println!("{category:?}: {inserted_variants} variants added");

        total_inserted_variants += inserted_variants;

        if inserted_variants > 0 {
            updater::update_variant_rank(loader, &config.family, variant_type, category).await?;
        }
    }

    Ok(total_inserted_variants)
}
