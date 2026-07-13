mod parse;
mod models;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;
use parse::vcf::process_vcf;

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

    process_vcf(&args.vcf, category, variant_type, &args.case_id);
}