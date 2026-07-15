use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::models::cytoband::Cytoband;

/// Loads cytoband annotations from a gzipped TSV file.
///
/// The cytoband file is expected to contain tab-separated entries with the
/// following format:
///
/// ```text
/// chromosome    start    end    cytoband_name    stain
/// chr1          0        2300000 p36.33          gneg
/// ```
///
/// Chromosome names are normalized by removing a leading `chr` prefix, so
/// both `chr1` and `1` are stored as `1`.
///
/// The annotations are stored in a chromosome-indexed map, where each
/// chromosome contains a list of genomic intervals associated with cytoband
/// names.
///
/// # Arguments
///
/// * `path` - Path to the gzipped cytoband TSV file.
///
/// # Returns
///
/// Returns a `HashMap` where:
///
/// * the key is the normalized chromosome name
/// * the value is a vector of cytoband intervals for that chromosome
///
/// # Errors
///
/// Returns an error if:
///
/// * the cytoband file cannot be opened
/// * the gzip stream cannot be read
/// * a start or end coordinate cannot be parsed as an integer
pub fn set_cytobands(
    path: &str,
) -> Result<HashMap<String, Vec<Cytoband>>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;

    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut cytobands: HashMap<String, Vec<Cytoband>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0]
            .trim_start_matches("chr")
            .to_string();

        let start: u64 = fields[1].parse()?;
        let end: u64 = fields[2].parse()?;

        cytobands
            .entry(chrom)
            .or_default()
            .push(Cytoband {
                start,
                end,
                name: fields[3].to_string(),
            });
    }

    Ok(cytobands)
}