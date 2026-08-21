use crate::models::case::CaseConfig;
use crate::parser::parse;
use clap::Parser;
use loader::Loader;
use mongodb::bson::Document;
use std::collections::{HashMap, HashSet};
use std::fs;

mod build;
mod config;
mod loader;
mod models;
mod parse;
mod parser;
mod utils;

pub struct VariantAnnotations<'a> {
    pub gene_to_panels: &'a HashMap<i32, HashSet<String>>,
    pub hgncid_to_gene: &'a HashMap<i32, Document>,
}

#[derive(Parser, Debug)]
struct Args {
    /// Path to the case configuration YAML file.
    #[arg(long = "case-config")]
    case_config: String,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let loader = Loader::new("config.toml").await?;

    let args = Args::parse();

    let yaml = fs::read_to_string(&args.case_config)?;
    let config: CaseConfig = serde_yaml::from_str(&yaml)?;

    println!("Case: {}", config.family);
    println!("Genome build: {}", config.human_genome_build);
    println!("Gene panels: {:?}", config.gene_panels);
    println!("Institute: {:?}", config.owner);

    if !loader.institute_exists(&config.owner).await? {
        return Err(format!("Institute '{}' does not exist in database", config.owner).into());
    }

    for sample in &config.samples {
        println!("Sample: {} ({:?})", sample.sample_id, sample.sample_name);
    }

    let panel_ids = config.gene_panels.as_deref().unwrap_or_default();

    let gene_to_panels = loader.gene_to_panels(panel_ids).await?;

    let hgncid_to_gene = loader.hgncid_to_gene(&config.human_genome_build).await?;

    println!(
        "Number of genes for genome build {}: {}",
        config.human_genome_build,
        hgncid_to_gene.len()
    );

    parse(
        &config,
        VariantAnnotations {
            gene_to_panels: &gene_to_panels,
            hgncid_to_gene: &hgncid_to_gene,
        },
    )?;

    Ok(())
}
