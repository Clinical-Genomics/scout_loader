use crate::Args;
use crate::models::build::GenomeBuild;
use crate::models::sample::SampleInfo;
use crate::models::variant::{VariantCategory, VariantType};
use crate::parse::cytobands::set_cytobands;
use crate::parse::vcf::process_vcf;
use std::collections::HashMap;

/// Parse command-line sample mappings in the format:
///
/// `SAMPLE_ID:DISPLAY_NAME:VCF_INDEX`
///
/// Example CLI usage:
///
/// `--samples ADM1059A1:NA12877:0 ADM1059A2:NA12882:1 ..
///
/// Creates `SampleInfo` entries containing the display name and VCF sample index.
pub fn parse_sample_mapping(
    samples: Option<Vec<String>>,
) -> Result<HashMap<String, SampleInfo>, String> {
    let mut mapping = HashMap::new();

    let Some(samples) = samples else {
        return Ok(mapping);
    };

    for sample in samples {
        let parts: Vec<&str> = sample.split(':').collect();

        if parts.len() != 3 {
            return Err(format!(
                "Invalid sample '{}'. Expected SAMPLE_ID:DISPLAY_NAME:VCF_POSITION",
                sample
            ));
        }

        let sample_id = parts[0].to_string();
        let display_name = parts[1].to_string();

        let vcf_index = parts[2].parse::<usize>().map_err(|_| {
            format!(
                "Invalid VCF position '{}' for sample '{}'",
                parts[2], sample_id
            )
        })?;

        mapping.insert(
            sample_id,
            SampleInfo {
                display_name,
                vcf_index,
            },
        );
    }

    Ok(mapping)
}

/// Parses and processes all VCFs provided for a case.
///
/// Shared information such as sample mappings, genome build, cytobands,
/// and variant type is prepared once and reused for every VCF.
pub fn parse(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let sample_mapping = parse_sample_mapping(args.samples)
        .map_err(|error| format!("Error parsing samples: {error}"))?;

    let variant_type = VariantType::from_str(&args.variant_type)
        .map_err(|_| format!("Invalid variant type: {}", args.variant_type))?;

    let genome_build = GenomeBuild::from_str(&args.genome_build)
        .map_err(|_| format!("Invalid genome build: {}", args.genome_build))?;

    let cytobands = set_cytobands(genome_build.cytoband_path())
        .map_err(|error| format!("Could not load cytobands: {error}"))?;

    let vcfs = [
        (args.snv, VariantCategory::Snv),
        (args.cancer, VariantCategory::Cancer),
        (args.sv, VariantCategory::Sv),
        (args.cancer_sv, VariantCategory::CancerSv),
        (args.fusion, VariantCategory::Fusion),
        (args.mei, VariantCategory::Mei),
        (args.str, VariantCategory::Str),
    ];

    for (vcf, category) in vcfs {
        if let Some(vcf) = vcf {
            process_vcf(
                vcf.to_str().ok_or("Invalid VCF path")?,
                category,
                variant_type,
                &args.case_id,
                &cytobands,
                &sample_mapping,
            );
        }
    }

    Ok(())
}
