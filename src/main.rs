mod parse;
use parse::vcf_parser::process_vcf;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <vcf_file>", args[0]);
        std::process::exit(1);
    }

    let vcf_path = &args[1];

    process_vcf(vcf_path);
}