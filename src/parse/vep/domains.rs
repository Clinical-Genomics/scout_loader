use crate::parse::vep::transcripts::VepFields;
use mongodb::bson::{Bson, Document};

/// Parse protein domain annotations from a VEP transcript entry.
///
/// Extracts supported protein domains from the VEP `DOMAINS` field.
pub fn parse_domains(transcript: &mut Document, entry: &VepFields<'_>) {
    if let Some(domains) = entry.get("DOMAINS").filter(|v| !v.is_empty()) {
        for annotation in domains.split('&') {
            let Some((domain_name, domain_id)) = annotation.split_once(':') else {
                continue;
            };

            match domain_name {
                "Pfam_domain" => {
                    transcript.insert("pfam_domain", Bson::String(domain_id.to_string()));
                }
                "PROSITE_profiles" => {
                    transcript.insert("prosite_profile", Bson::String(domain_id.to_string()));
                }
                "SMART_domains" => {
                    transcript.insert("smart_domain", Bson::String(domain_id.to_string()));
                }
                "hmmpanther" => {
                    transcript.insert("panther_domain", Bson::String(domain_id.to_string()));
                }
                _ => {}
            }
        }
    }
}
