use clap::Parser;
use scout_loader::loader::Loader;
use scout_loader::models::case::CaseConfig;
use scout_loader::models::variant::VariantAnnotations;
use scout_loader::parser::parse;
use std::fs;

#[derive(Parser, Debug)]
struct Args {
    /// Path to the case configuration YAML file.
    #[arg(long = "case-config")]
    case_config: String,

    /// Path to the MongoDB configuration TOML file.
    #[arg(long)]
    config: Option<String>,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let config_path = args.config.as_deref().unwrap_or("config.toml");

    let loader = Loader::new(config_path).await?;

    let yaml = fs::read_to_string(&args.case_config)?;
    let config: CaseConfig = serde_yaml::from_str(&yaml)?;

    println!("Case: {}", config.family);
    println!("Genome build: {}", config.human_genome_build);
    println!("Gene panels: {:?}", config.gene_panels);
    println!("Institute: {:?}", config.owner);

    validate_institute(&loader, &config.owner).await?;

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

    let total_inserted_variants = parse(
        &config,
        VariantAnnotations {
            gene_to_panels: &gene_to_panels,
            hgncid_to_gene: &hgncid_to_gene,
        },
        &loader,
    )
    .await?;

    println!(
        "Total variants added for case {}: {}",
        config.family, total_inserted_variants
    );

    Ok(())
}

/// Validates that the specified institute exists in the database.
///
/// In the test environment (`TEST_ENV`), the database check is skipped.
/// Returns an error if the institute does not exist.
async fn validate_institute(
    loader: &Loader,
    institute_id: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if std::env::var("TEST_ENV").is_ok() {
        return Ok(());
    }

    if !loader.institute_exists(institute_id).await? {
        return Err(format!("Institute '{}' does not exist in database", institute_id).into());
    }

    Ok(())
}
