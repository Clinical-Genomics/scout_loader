mod parse;
mod models;
mod utils;
use clap::Parser;
use std::collections::HashMap;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;
use crate::models::build::GenomeBuild;
use crate::models::sample::SampleInfo;
use parse::vcf::process_vcf;
use parse::cytobands::set_cytobands;


/// Parse command-line sample mappings in the format:
///
/// `SAMPLE_ID:DISPLAY_NAME:VCF_INDEX`
///
/// Example CLI usage:
///
/// `--samples ADM1059A1:NA12881:0 ADM1059A2:NA12882:1`
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

        let vcf_index = parts[2]
            .parse::<usize>()
            .map_err(|_| {
                format!(
                    "Invalid VCF position '{}' for sample '{}'",
                    parts[2],
                    sample_id
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


#[derive(Parser)]
struct Args {
    /// Path to the VCF containing the variants
    #[arg(long)]
    vcf: String,

    /// 'snv', 'sv', 'str', ..
    #[arg(long)]
    category: String,

    /// 'clinical' or 'research'
    #[arg(long)]
    variant_type: String,

    /// Case id 
    #[arg(long)]
    case_id: String,

    /// 'GRCh37' or 'GRCh38'
    #[arg(long)]
    genome_build: String,

    // Sample IDs to extract genotypes for, formatted as:
    // 'SAMPLE_ID:DISPLAY_NAME:VCF_POSITION'
    //
    // Example:
    // ADM1059A1:NA12881:0 ADM1059A2:NA12882:1
    #[arg(short, long, num_args = 1..)]
    samples: Option<Vec<String>>,
}

/// Reads variants from a VCF file and converts each record into a Variant.
///
/// The parser behavior depends on the provided variant category
/// (e.g. "snv", "sv", "str").
///
/// # Arguments
///
/// * `path` - Path to the input VCF file.
/// * `category` - Variant category used to select the appropriate parser.
/// * `variant_type` - Variant type (clinical or research).
/// * `case_id` - Case _id to be saved in the documents
fn main() {
    let args = Args::parse();

    let sample_mapping = match parse_sample_mapping(args.samples) {
        Ok(mapping) => mapping,
        Err(e) => {
            eprintln!("Error parsing samples: {}", e);
            return;
        }
    };

    let category = VariantCategory::from_str(&args.category)
        .expect("Invalid category");
    let variant_type = VariantType::from_str(&args.variant_type)
        .expect("Invalid variant type");
    let genome_build = GenomeBuild::from_str(&args.genome_build)
        .expect("Invalid genome build");

    let cytobands = match set_cytobands(genome_build.cytoband_path()) {
        Ok(value) => value,
        Err(error) => {
            eprintln!("Could not load cytobands: {}", error);
            return;
        }
    };

    process_vcf(&args.vcf, category, variant_type, &args.case_id, &cytobands, &sample_mapping);
}