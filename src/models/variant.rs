#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub enum VariantCategory {
    Snv,
    Sv,
    Str,
}

impl VariantCategory {
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "snv" => Ok(VariantCategory::Snv),
            "sv" => Ok(VariantCategory::Sv),
            "str" => Ok(VariantCategory::Str),
            _ => Err(format!("Unknown variant category: {}", s)),
        }
    }
}

/// Represents a genomic variant (snv).
#[derive(Debug)]
pub struct Variant {
    pub chromosome: String,
    pub position: u64,
    pub end: u64,
    pub reference: String,
    pub alternative: String
}
