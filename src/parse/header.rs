use mongodb::bson::Document;
use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Extract the VEP CSQ annotation fields from the VCF header.
///
/// The VEP annotations are stored in the CSQ INFO field. The field names are
/// defined after `Format:` in the CSQ description and are separated by `|`.
/// Returned field names are converted to uppercase.
///
/// Returns an empty vector if no CSQ INFO header is found.
pub fn parse_vep_header(header: &HeaderView) -> Vec<String> {
    for record in header.header_records() {
        if let HeaderRecord::Info { values, .. } = record {
            if let Some(id) = values.get("ID") {
                if id == "CSQ" {
                    if let Some(description) = values.get("Description") {
                        if let Some(format) = description.split("Format:").nth(1) {
                            return format
                                .trim()
                                .trim_matches('"')
                                .trim_end_matches('>')
                                .split('|')
                                .map(|field| field.to_uppercase())
                                .collect();
                        }
                    }
                }
            }
        }
    }

    Vec::new()
}

/// Parse LoqusDB archive metadata from a VCF header.
///
/// Extracts metadata written by loqusdb in the VCF header:
///
/// - `INFO/Obs` description -> `Description`
/// - `##NrCases` value -> `NrCases`
/// - `##Software=<ID=loqusdb,...>` date -> `Date`
///
/// Returns `None` if no LoqusDB metadata is found.
///
/// # Arguments
///
/// * `path` - Path to the (possibly gzipped) VCF file.
///
/// # Returns
///
/// A BSON document containing local archive metadata, or `None` if no
/// metadata could be extracted.
pub fn parse_local_archive_header(path: &str) -> Option<Document> {
    let file = File::open(path).ok()?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut local_archive_info = Document::new();

    for line in reader.lines().flatten() {
        if line.starts_with("##INFO=<ID=Obs") {
            if let Some(description) = line.split("Description=\"").nth(1) {
                if let Some(description) = description.split('"').next() {
                    local_archive_info.insert("Description", description.to_string());
                }
            }
        }

        if line.starts_with("##NrCases=") {
            if let Some(value) = line.strip_prefix("##NrCases=") {
                if let Ok(value) = value.parse::<i32>() {
                    local_archive_info.insert("NrCases", value);
                }
            }
        }

        if line.starts_with("##Software=<ID=loqusdb") {
            if let Some(date) = line.split("Date=\"").nth(1) {
                if let Some(date) = date.split('"').next() {
                    local_archive_info.insert("Date", date.to_string());
                }
            }
        }

        // Stop reading once the header is finished.
        if line.starts_with("#CHROM") {
            break;
        }
    }

    (!local_archive_info.is_empty()).then_some(local_archive_info)
}