use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Clinical,
    Research
}

impl VariantType {
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "clinical" => Ok(VariantType::Clinical),
            "research" => Ok(VariantType::Research),
            _ => Err(format!("Unknown variant type: {}", s)),
        }
    }    
}

impl fmt::Display for VariantType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VariantType::Clinical => write!(f, "clinical"),
            VariantType::Research => write!(f, "research"),
        }
    }
}

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

/// Represents the identifiers associated with a variant.
pub struct VariantIds {
    pub simple_id: String,
    pub variant_id: String,
    pub display_name: String,
    pub document_id: String,
}

#[derive(Debug)]
pub struct Compound {
    pub display_name: String,
    pub variant: String,
    pub score: f64,
}

/// Represents a genomic variant (snv).
#[derive(Debug)]
pub struct Variant {
    pub simple_id: String,
    pub variant_id: String,
    pub display_name: String,
    pub document_id: String,
    pub case_id: String,
    pub rank_score: Option<i32>,
    pub norm_rank_score: Option<i32>,
    pub r#type: String,
    pub chromosome: String,
    pub position: u64,
    pub end: u64,
    pub reference: String,
    pub alternative: String,
    pub filters: Vec<String>,
    pub quality: f32,
    pub compounds: Vec<Compound>,
}
