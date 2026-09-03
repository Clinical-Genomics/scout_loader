use crate::loader::Loader;
use crate::models::build::GenomeBuild;
use crate::models::case::CaseConfig;
use crate::models::variant::VariantAnnotations;
use crate::models::variant::{VariantCategory, VariantType};
use crate::parse::cytobands::set_cytobands;
use crate::parse::vcf::process_vcf;
use crate::updater;
use std::path::PathBuf;
use std::str::FromStr;

/// Select VCFs to process based on variant categories and data type.
///
/// By default, clinical VCFs are selected. When research is true, only
/// research VCFs are considered instead. If categories is provided, only
/// VCFs belonging to the requested variant categories are selected.
///
/// VCFs that are not specified in the case configuration are skipped.
///
/// Returns the selected VCF paths together with their corresponding variant
/// categories.
pub fn select_vcfs<'a>(
    config: &'a CaseConfig,
    categories: Option<&[String]>,
    research: bool,
) -> Vec<(&'a PathBuf, VariantCategory)> {
    let vcfs = [
        (
            if research {
                &config.vcf_snv_research
            } else {
                &config.vcf_snv
            },
            VariantCategory::Snv,
            "snv",
        ),
        (
            if research {
                &config.vcf_sv_research
            } else {
                &config.vcf_sv
            },
            VariantCategory::Sv,
            "sv",
        ),
        // ...
    ];

    vcfs.into_iter()
        .filter_map(|(vcf, category, name)| {
            if let Some(categories) = categories
                && !categories.iter().any(|requested| requested == name)
            {
                return None;
            }

            vcf.as_ref().map(|path| (path, category))
        })
        .collect()
}

/// Parses and processes the selected VCFs provided for a case.
///
/// Shared information such as the genome build and cytobands is prepared once
/// and reused for every VCF. Sample mappings are created separately for each
/// VCF because sample indices may differ between VCF files.
///
/// VCF selection is handled by `select_vcfs`, which determines whether
/// clinical or research VCFs should be loaded and optionally filters them by
/// variant category.
///
/// After each VCF is processed, the number of newly inserted variants is
/// checked. If variants were inserted, their `variant_rank` values are
/// updated for the corresponding variant category and type.
///
/// If an error occurs while processing a VCF or updating variant ranks, all
/// variants loaded for the case are removed before the error is returned.
///
/// Returns the total number of variants inserted across all selected VCFs
/// for the case.
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

    let vcfs = select_vcfs(config, categories, research);

    let variant_type = if research {
        VariantType::Research
    } else {
        VariantType::Clinical
    };

    let mut total_inserted_variants = 0;

    for (vcf, category) in vcfs {
        annotations.managed_variant_ids = loader
            .get_managed_variant_ids(&category.to_string(), &config.human_genome_build)
            .await?;

        let inserted_variants = match process_vcf(
            vcf.to_str().ok_or("Invalid VCF path")?,
            category,
            variant_type,
            config,
            &cytobands,
            &annotations,
            loader,
        )
        .await
        {
            Ok(count) => count,
            Err(error) => {
                loader.delete_case_variants(&config.family).await?;
                return Err(error);
            }
        };

        println!("{variant_type:?} {category:?}: {inserted_variants} variants added");

        total_inserted_variants += inserted_variants;

        if inserted_variants > 0 {
            if let Err(error) =
                updater::update_variant_rank(loader, &config.family, variant_type, category).await
            {
                loader.delete_case_variants(&config.family).await?;
                return Err(error);
            }
        }
    }

    Ok(total_inserted_variants)
}
