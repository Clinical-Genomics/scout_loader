use rust_htslib::bcf::header::{HeaderRecord, HeaderView};

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