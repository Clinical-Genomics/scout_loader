use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenomeBuild {
    Grch37,
    Grch38,
}

impl FromStr for GenomeBuild {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "grch37" | "37" | "hg19" => Ok(GenomeBuild::Grch37),
            "grch38" | "38" | "hg38" => Ok(GenomeBuild::Grch38),
            _ => Err(format!("Unknown genome build: {s}")),
        }
    }
}

impl GenomeBuild {
    pub fn cytoband_path(&self) -> &'static str {
        match self {
            GenomeBuild::Grch37 => "resources/cytoBand_hg19.txt.gz",
            GenomeBuild::Grch38 => "resources/cytoBand_hg38.txt.gz",
        }
    }
}
