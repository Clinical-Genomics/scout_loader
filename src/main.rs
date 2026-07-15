mod parse;
mod models;
mod utils;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;
use crate::models::build::GenomeBuild;
use parse::vcf::process_vcf;
use parse::cytobands::set_cytobands;

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(long)]
    vcf: String,

    #[arg(long)]
    category: String,

    #[arg(long)]
    variant_type: String,

    #[arg(long)]
    case_id: String,

    #[arg(long)]
    genome_build: String,
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

    process_vcf(&args.vcf, category, variant_type, &args.case_id, genome_build.as_str(), &cytobands);
}