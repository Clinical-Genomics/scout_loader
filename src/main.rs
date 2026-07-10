mod parse;
mod models;
use parse::vcf::process_vcf;

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(long)]
    vcf: String,

    #[arg(long)]
    category: String,
}

fn main() {
    let args = Args::parse();

    println!("VCF: {}", args.vcf);
    println!("Category: {}", args.category);

    process_vcf(&args.vcf, &args.category);
}