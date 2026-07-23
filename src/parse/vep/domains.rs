/// Parse protein domain annotations from a VEP transcript entry.
///
/// Extracts supported protein domains from the VEP `DOMAINS` field.
fn parse_domains(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if let Some(domains) = entry.get("DOMAINS").filter(|v| !v.is_empty()) {
        for annotation in domains.split('&') {
            let parts: Vec<&str> = annotation.split(':').collect();

            if parts.len() < 2 {
                continue;
            }

            let domain_name = parts[0];
            let domain_id = parts[1];

            match domain_name {
                "Pfam_domain" => {
                    transcript.insert(
                        "pfam_domain",
                        Bson::String(domain_id.to_string()),
                    );
                }
                "PROSITE_profiles" => {
                    transcript.insert(
                        "prosite_profile",
                        Bson::String(domain_id.to_string()),
                    );
                }
                "SMART_domains" => {
                    transcript.insert(
                        "smart_domain",
                        Bson::String(domain_id.to_string()),
                    );
                }
                "hmmpanther" => {
                    transcript.insert(
                        "panther_domain",
                        Bson::String(domain_id.to_string()),
                    );
                }
                _ => {}
            }
        }
    }
}