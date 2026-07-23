use rust_htslib::bcf::Record;
use mongodb::bson::{Bson, Document};
use std::collections::HashSet;
use crate::HashMap;
use crate::models::gene::GeneAnnotation;
use crate::parse::vep::annotations::{get_hgnc_id, get_strand, get_regional_annotation, parse_mane_annotations, parse_superdups_fracmatch, parse_clinvar_annotations, parse_dbsnp, parse_cosmic};
use crate::parse::vep::predictors::{get_prediction, parse_transcripts_spliceai};
use crate::parse::vep::utils::{get_highest_float_score_in_string, get_sequence_aux};
use crate::parse::vep::domains::parse_domains;
use crate::parse::vep::scores::parse_cadd;
use crate::parse::vep::frequencies::{parse_mt_frequencies, parse_variant_frequencies};
use crate::parse::info::parse_info_string;
use crate::parse::vep::genes::parse_genes;

/// Parse VEP CSQ annotations from a VCF record.
///
/// Extracts transcript annotations from the CSQ INFO field using the VEP
/// header, collects dbSNP and COSMIC identifiers, and builds gene-level
/// information separately from transcript-level information.
///
/// Gene identifiers are not stored in transcripts to avoid duplication.
/// They are passed separately to `parse_genes`.
///
/// Returns the parsed transcripts and gene annotations.
pub fn parse_vep_transcripts(
    record: &Record,
    vep_header: &[String],
    variant: &mut Document,
) -> Vec<Document> {
    let mut parsed_transcripts = Vec::new();

    let mut dbsnp_ids = HashSet::new();
    let mut cosmic_ids = HashSet::new();

    if vep_header.is_empty() {
        return parsed_transcripts;
    }


    if let Ok(Some(csq)) = record.info(b"CSQ").string() {
        let csq_string = csq
            .iter()
            .map(|value| String::from_utf8_lossy(value))
            .collect::<Vec<_>>()
            .join(",");
        for transcript_info in csq_string.split(',') {
            let raw_transcript: HashMap<String, String> = vep_header
                .iter()
                .zip(transcript_info.split('|'))
                .map(|(key, value)| {
                    (key.clone(), value.to_string())
                })
                .collect();

            if let Some(transcript) = parse_vep_transcript(raw_transcript) {
                if let Ok(values) = transcript.get_array("dbsnp") {
                    for value in values {
                        if let Bson::String(id) = value {
                            dbsnp_ids.insert(id.clone());
                        }
                    }
                }

                if let Ok(values) = transcript.get_array("cosmic") {
                    for value in values {
                        if let Bson::String(id) = value {
                            cosmic_ids.insert(id.clone());
                        }
                    }
                }

                parsed_transcripts.push(transcript);
            }
        }

        // The COSMIC INFO tag may be added by VEP and/or bcftools annotate.
        if let Some(cosmic_tag) = parse_info_string(record, b"COSMIC") {
            cosmic_ids.extend(
                cosmic_tag
                    .split('&')
                    .filter(|id| !id.is_empty())
                    .map(str::to_string),
            );
        }

        if !cosmic_ids.is_empty() {
            variant.insert(
                "cosmic_ids",
                Bson::Array(
                    cosmic_ids
                        .iter()
                        .cloned()
                        .map(Bson::String)
                        .collect(),
                ),
            );
        }

        if !dbsnp_ids.is_empty() && !variant.contains_key("dbsnp_id") {
            variant.insert(
                "dbsnp_id",
                Bson::String(
                    dbsnp_ids.iter().cloned().collect::<Vec<_>>().join(";"),
                ),
            );
        }


        let genes = parse_genes(&parsed_transcripts);
        variant.insert(
            "genes",
            Bson::Array(
                genes
                    .into_iter()
                    .map(Bson::Document)
                    .collect(),
            ),
        );

    }

     parsed_transcripts
}

/// Parse a single VEP transcript annotation.
///
/// Extracts transcript-specific information from a single CSQ entry.
/// Gene-level information such as HGNC ID and gene symbol is handled
/// separately by the gene parser.
///
/// Returns None if no transcript ID is available.
pub fn parse_vep_transcript(
    entry: HashMap<String, String>,
) -> Option<Document> {

    let transcript_id = entry
        .get("FEATURE")
        .map(|id| id.split(':').next().unwrap_or(""))
        .unwrap_or("");

    if transcript_id.is_empty() {
        return None;
    }
    
    let mut transcript = Document::new();

    transcript.insert(
        "transcript_id",
        transcript_id.to_string(),
    );

    transcript.insert(
        "protein_id",
        entry
            .get("ENSP")
            .cloned()
            .unwrap_or_default(),
    );

    let polyphen = get_prediction(&entry, &["POLYPHEN"]);

    transcript.insert(
        "polyphen_prediction",
        Bson::String(polyphen),
    );

    let sift = get_prediction(&entry, &["SIFT", "SIFT_PRED"]);
    transcript.insert(
        "sift_prediction",
        Bson::String(sift),
    );

    if let Some(value) = entry.get("REVEL_RANKSCORE").filter(|v| !v.is_empty()) {
        if let Some(rankscore) = get_highest_float_score_in_string(value) {
            transcript.insert(
                "revel_rankscore",
                Bson::Double(rankscore),
            );
        }
    }

    if let Some(value) = entry.get("REVEL_SCORE").filter(|v| !v.is_empty()) {
        if let Some(score) = get_highest_float_score_in_string(value) {
            transcript.insert(
                "revel_raw_score",
                Bson::Double(score),
            );
        }
    }

    parse_transcripts_spliceai(&mut transcript, &entry);

    transcript.insert(
        "swiss_prot",
        Bson::String(
            entry
                .get("SWISSPROT")
                .filter(|value| !value.is_empty())
                .cloned()
                .unwrap_or_else(|| "unknown".to_string()),
        ),
    );


    for (source, target) in [
        ("GERP++_RS", "gerp"),
        ("PHASTCONS100WAY_VERTEBRATE", "phast"),
        ("PHYLOP100WAY_VERTEBRATE", "phylop"),
    ] {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            transcript.insert(
                target,
                Bson::String(value.clone()),
            );
        }
    }

    parse_domains(&mut transcript, &entry);

    if let Some(coding_sequence) = get_sequence_aux(&entry, "HGVSC") {
        transcript.insert(
            "coding_sequence_name",
            Bson::String(coding_sequence),
        );
    }

    if let Some(coding_sequence) = get_sequence_aux(&entry, "HGVSP") {
        transcript.insert(
            "protein_sequence_name",
            Bson::String(coding_sequence),
        );
    }

    transcript.insert(
        "biotype",
        Bson::String(
            entry.get("BIOTYPE")
                .cloned()
                .unwrap_or_default(),
        ),
    );

    for (source, target) in [("EXON", "exon"),("INTRON", "intron")] {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            transcript.insert(
                target,
                Bson::String(value.clone()),
            );
        }
    }

    if let Some(strand) = get_strand(&entry) {
        transcript.insert(
            "strand",
            Bson::String(strand),
        );
    }


    let functional_annotations: Vec<String> = entry
        .get("CONSEQUENCE")
        .filter(|value| !value.is_empty())
        .map(|value| value.split('&').map(String::from).collect())
        .unwrap_or_default();

    if !functional_annotations.is_empty() {
        transcript.insert(
            "functional_annotations",
            Bson::Array(
                functional_annotations
                    .iter()
                    .map(|annotation| Bson::String(annotation.clone()))
                    .collect(),
            ),
        );

        let region_annotations = get_regional_annotation(&functional_annotations);

        transcript.insert(
            "region_annotations",
            Bson::Array(
                region_annotations
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }

    transcript.insert(
        "is_canonical",
        Bson::Boolean(
            entry.get("CANONICAL")
                .map(|value| value == "YES")
                .unwrap_or(false),
        ),
    );

    parse_mane_annotations(&mut transcript, &entry);

    parse_cadd(&mut transcript, &entry);

    parse_superdups_fracmatch(&mut transcript, &entry);

    parse_mt_frequencies(&mut transcript, &entry);

    parse_variant_frequencies(&mut transcript, &entry);

    parse_clinvar_annotations(&mut transcript, &entry);

    parse_dbsnp(&mut transcript, &entry);

    parse_cosmic(&mut transcript, &entry);

    Some(transcript)
}