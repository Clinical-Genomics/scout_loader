#[derive(Debug)]

/// Represents a genomic variant (snv).
pub struct Variant {
    pub chromosome: String,
    pub position: u64,
    pub end: u64,
}
