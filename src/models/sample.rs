#[derive(Debug, Clone)]
pub struct SampleInfo {
    /// Human-readable display name for the sample.
    pub display_name: String,

    /// Zero-based index of the sample in the VCF header.
    pub vcf_index: usize,
}