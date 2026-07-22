use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use rust_htslib::bcf::Record;
use mongodb::bson::{Document};


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

/// Add gene and transcript information from VEP annotations.
///
/// Parses CSQ annotations using the VEP header extracted from the VCF header,
/// extracts transcript and gene information, and collects dbSNP and COSMIC
/// identifiers. Fusion variants are skipped because they have their own
/// dedicated gene/transcript parser.
///
/// Returns the parsed transcripts for later processing.
pub fn parse_vep_transcripts(
    record: &Record,
    vep_header: &[String],
    variant: &mut Document,
) -> Vec<Document> {
    
    let mut parsed_transcripts = Vec::new();

    /*
    let mut dbsnp_ids = HashSet::new();
    let mut cosmic_ids = HashSet::new();

    if !vep_header.is_empty() {
        if let Some(csq) = parse_info_string(record, b"CSQ") {
            for transcript_info in csq.split(',') {
                let raw_transcript: HashMap<String, String> = vep_header
                    .iter()
                    .zip(transcript_info.split('|'))
                    .map(|(key, value)| (key.clone(), value.to_string()))
                    .collect();

                if let Some(transcript) = parse_transcript(raw_transcript) {
                    collect_ids(
                        &transcript,
                        &mut dbsnp_ids,
                        &mut cosmic_ids,
                    );

                    parsed_transcripts.push(transcript);
                }
            }
        }
    }

    // COSMIC can also be annotated outside VEP
    if let Some(cosmic_tag) = parse_info_string(record, b"COSMIC") {
        for cosmic_id in cosmic_tag.split('&') {
            cosmic_ids.insert(cosmic_id.to_string());
        }
    }

    if !dbsnp_ids.is_empty() && !variant.contains_key("dbsnp_id") {
        variant.insert(
            "dbsnp_id",
            dbsnp_ids.into_iter().collect::<Vec<_>>().join(";"),
        );
    }

    if !cosmic_ids.is_empty() {
        variant.insert(
            "cosmic_ids",
            cosmic_ids.into_iter().collect::<Vec<_>>(),
        );
    }

    let genes = parse_genes(&parsed_transcripts);

    let mut hgnc_ids = HashSet::new();

    for gene in &genes {
        if let Some(id) = gene.get("hgnc_id") {
            hgnc_ids.insert(id.clone());
        }
    }

    // Stranger STR annotations
    if let Some(str_hgnc_id) = parse_info_string(record, b"HGNCId") {
        hgnc_ids.insert(str_hgnc_id.clone());

        if genes.is_empty() {
            variant.insert(
                "genes",
                vec![doc! {
                    "hgnc_id": str_hgnc_id
                }],
            );
        }
    } else {
        variant.insert("genes", genes);
    }

    variant.insert(
        "hgnc_ids",
        hgnc_ids.into_iter().collect::<Vec<_>>(),
    );
	*/

    parsed_transcripts
	
}

