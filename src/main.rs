use clap::Parser;
use std::collections::HashMap;

use crate::models::case::CaseConfig;
use crate::parser::parse;
use loader::Loader;
use std::fs;

mod config;
mod loader;
mod models;
mod parse;
mod parser;
mod utils;

#[derive(Parser, Debug)]
struct Args {
    /// Path to the case configuration YAML file.
    #[arg(long = "case-config")]
    case_config: String,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let yaml = fs::read_to_string(&args.case_config)?;
    let config: CaseConfig = serde_yaml::from_str(&yaml)?;

    println!("Case: {}", config.family);
    println!("Genome build: {}", config.human_genome_build);
    println!("Gene panels: {:?}", config.gene_panels);

    for sample in &config.samples {
        println!("Sample: {} ({:?})", sample.sample_id, sample.sample_name);
    }

    let loader = Loader::new("config.toml").await?;

    let panel_ids = config.gene_panels.as_deref().unwrap_or_default();

    let gene_to_panels = loader.gene_to_panels(panel_ids).await?;

    let hgncid_to_gene = loader.hgncid_to_gene(&config.human_genome_build).await?;

    println!(
        "Number of genes for genome build {}: {}",
        config.human_genome_build,
        hgncid_to_gene.len()
    );

    parse(&config, &gene_to_panels, &hgncid_to_gene)?;

    Ok(())
}
