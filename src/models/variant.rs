use std::fmt;
use serde::Serialize;

#[derive(Debug)]
pub struct Coordinates {
    pub chromosome: String,
    pub position: u64,
    pub end: u64,
    pub end_chrom: String,
    pub length: i64,
    pub sub_category: String,
    pub mate_id: Option<String>,
    pub cytoband_start: Option<String>,
    pub cytoband_end: Option<String>

}

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
    Cancer,
    Sv,
    CancerSv,
    Fusion,
    Mei,
    Str,
}

impl VariantCategory {
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "snv" => Ok(VariantCategory::Snv),
            "cancer" => Ok(VariantCategory::Cancer),
            "sv" => Ok(VariantCategory::Sv),
            "cancer_sv" => Ok(VariantCategory::CancerSv),
            "fusion" => Ok(VariantCategory::Fusion),
            "mei" => Ok(VariantCategory::Mei),
            "str" => Ok(VariantCategory::Str),
            _ => Err(format!("Unknown variant category: {}", s)),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            VariantCategory::Snv => "snv",
            VariantCategory::Cancer => "cancer",
            VariantCategory::Sv => "sv",
            VariantCategory::CancerSv => "cancer_sv",
            VariantCategory::Fusion => "fusion",
            VariantCategory::Mei => "mei",
            VariantCategory::Str => "str",
        }
    }
}

impl fmt::Display for VariantCategory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

// Represents the identifiers associated with a variant.
pub struct VariantIds {
    pub simple_id: String,
    pub variant_id: String,
    pub display_name: String,
    pub document_id: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct Compound {
    pub display_name: String,
    pub variant: String,
    pub score: f64,
}