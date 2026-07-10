mod parse;
mod models;
use crate::models::variant::VariantCategory;
use parse::vcf::process_vcf;

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(long)]
    vcf: String,

    #[arg(long)]
    category: String,
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
fn main() {
    let args = Args::parse();

    let category = VariantCategory::from_str(&args.category)
        .expect("Invalid category");

    process_vcf(&args.vcf, category);
}