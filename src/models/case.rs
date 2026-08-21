use serde::Deserialize;
use std::path::PathBuf;

/// Configuration required to load a Scout case from a YAML file.
///
/// The YAML file may contain many additional Scout configuration fields,
/// but only the fields defined here are deserialized. Research VCFs are
/// intentionally excluded and will not be loaded.
#[derive(Debug, Deserialize)]
pub struct CaseConfig {
    pub owner: String,

    pub family: String,

    pub human_genome_build: String,

    pub gene_panels: Option<Vec<String>>,

    pub samples: Vec<SampleConfig>,

    pub vcf_snv: Option<PathBuf>,

    pub vcf_sv: Option<PathBuf>,

    pub vcf_str: Option<PathBuf>,

    pub vcf_mei: Option<PathBuf>,

    pub vcf_cancer: Option<PathBuf>,

    pub vcf_cancer_sv: Option<PathBuf>,

    pub vcf_fusion: Option<PathBuf>,
}

/// Sample information extracted from the case YAML configuration.
///
/// This describes how a sample is identified in the case configuration.
/// The position of the sample in a VCF is determined separately when the
/// VCF is parsed and represented by [`SampleInfo`].
#[derive(Debug, Deserialize)]
pub struct SampleConfig {
    /// Sample identifier used in the case configuration.
    pub sample_id: String,

    /// Human-readable sample name, if provided.
    pub sample_name: Option<String>,
}
