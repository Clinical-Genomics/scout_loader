use clap::Parser;
use std::collections::HashMap;
use std::path::PathBuf;

use loader::Loader;

mod config;
mod loader;
mod models;
mod parse;
mod parser;
mod utils;

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Args {
    /// Case ID.
    #[arg(long)]
    pub case_id: String,

    /// Sample mappings in the format SAMPLE_ID:DISPLAY_NAME:VCF_POSITION.
    ///
    /// Example: SAMPLE1:NA12877:0 SAMPLE2:NA12878:1
    #[arg(short, long, num_args = 1..)]
    pub samples: Option<Vec<String>>,

    /// Variant type, e.g. clinical or research.
    #[arg(long)]
    pub variant_type: String,

    /// Genome build, e.g. 37 or 38.
    #[arg(long)]
    pub genome_build: String,

    /// SNV VCF.
    #[arg(long)]
    pub snv: Option<PathBuf>,

    /// Cancer VCF.
    #[arg(long)]
    pub cancer: Option<PathBuf>,

    /// Structural variant VCF.
    #[arg(long)]
    pub sv: Option<PathBuf>,

    /// Cancer structural variant VCF.
    #[arg(long = "cancer-sv")]
    pub cancer_sv: Option<PathBuf>,

    /// Fusion VCF.
    #[arg(long)]
    pub fusion: Option<PathBuf>,

    /// Mobile element insertion VCF.
    #[arg(long)]
    pub mei: Option<PathBuf>,

    /// STR VCF.
    #[arg(long)]
    pub str: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let loader = Loader::new("config.toml")?;

    // Connect to the database and retrieve gene_to_panels and hgncid_to_gene

    // Orchestrate loading
    parser::parse(args)?;

    Ok(())
}
