use crate::parse::info::parse_info_string;
use crate::parse::vep::annotations::{
    get_regional_annotation, get_strand, parse_clinvar_annotations, parse_cosmic, parse_dbsnp,
    parse_mane_annotations, parse_superdups_fracmatch,
};
use crate::parse::vep::domains::parse_domains;
use crate::parse::vep::frequencies::{parse_mt_frequencies, parse_variant_frequencies};
use crate::parse::vep::predictors::{get_prediction, parse_transcripts_spliceai};
use crate::parse::vep::scores::parse_cadd;
use crate::parse::vep::utils::{get_highest_float_score_in_string, get_sequence_aux};
use mongodb::bson::{Bson, Document};
use rust_htslib::bcf::Record;
use std::collections::HashMap;
use std::collections::HashSet;

pub struct VepFields<'a> {
    indices: &'a HashMap<String, usize>,
    values: &'a [&'a str],
}

impl<'a> VepFields<'a> {
    pub(crate) fn new(indices: &'a HashMap<String, usize>, values: &'a [&'a str]) -> Self {
        Self { indices, values }
    }

    pub(crate) fn get(&self, key: &str) -> Option<&'a str> {
        self.indices
            .get(key)
            .and_then(|&index| self.values.get(index).copied())
    }

    pub(crate) fn contains_key(&self, key: &str) -> bool {
        self.indices
            .get(key)
            .is_some_and(|&index| index < self.values.len())
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (&str, &str)> + '_ {
        self.indices
            .iter()
            .filter_map(|(key, &index)| self.values.get(index).map(|value| (key.as_str(), *value)))
    }
}

pub fn parse_vep_transcripts(
    record: &Record,
    field_indices: &HashMap<String, usize>,
    variant: &mut Document,
) -> Vec<Document> {
    let mut parsed_transcripts = Vec::new();

    if field_indices.is_empty() {
        return parsed_transcripts;
    }

    let mut dbsnp_ids = HashSet::new();
    let mut cosmic_ids = HashSet::new();

    if let Ok(Some(csq)) = record.info(b"CSQ").string() {
        for value in csq.iter() {
            let csq_string = String::from_utf8_lossy(value);
            let mut field_values: Vec<&str> = Vec::new();

            for transcript_info in csq_string.split(',') {
                field_values.clear();
                field_values.extend(transcript_info.split('|'));

                if let Some(transcript) =
                    parse_vep_transcript_from_values(field_indices, &field_values)
                {
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
        }

        // COSMIC identifiers may also come from VEP or bcftools annotate.
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
                Bson::Array(cosmic_ids.into_iter().map(Bson::String).collect()),
            );
        }

        if !dbsnp_ids.is_empty() && !variant.contains_key("dbsnp_id") {
            variant.insert(
                "dbsnp_id",
                Bson::String(dbsnp_ids.into_iter().collect::<Vec<_>>().join(";")),
            );
        }
    }

    parsed_transcripts
}

fn parse_vep_transcript_from_values(
    field_indices: &HashMap<String, usize>,
    field_values: &[&str],
) -> Option<Document> {
    let fields = VepFields::new(field_indices, field_values);
    parse_vep_transcript(&fields)
}

/// Parse a single VEP transcript annotation.
///
/// Extracts transcript-specific information from a single CSQ entry.
/// Gene-level information such as HGNC ID and gene symbol is handled
/// separately by the gene parser.
///
/// Returns None if no transcript ID is available.
pub fn parse_vep_transcript(entry: &VepFields<'_>) -> Option<Document> {
    let transcript_id = entry
        .get("FEATURE")
        .map(|id| id.split(':').next().unwrap_or(""))
        .unwrap_or("");

    if transcript_id.is_empty() {
        return None;
    }

    let mut transcript = Document::new();

    if let Some(hgnc_id) = entry.get("HGNC_ID")
        && !hgnc_id.is_empty()
    {
        let hgnc_id = hgnc_id
            .split(':')
            .next_back()
            .and_then(|value| value.parse::<i32>().ok());

        if let Some(hgnc_id) = hgnc_id {
            transcript.insert("hgnc_id", Bson::Int32(hgnc_id));
        }
    }

    if let Some(symbol) = entry.get("SYMBOL")
        && !symbol.is_empty()
    {
        transcript.insert("hgnc_symbol", Bson::String(symbol.to_string()));
    }

    transcript.insert("transcript_id", transcript_id.to_string());
    transcript.insert("protein_id", entry.get("ENSP").unwrap_or("").to_string());

    let polyphen = get_prediction(entry, &["POLYPHEN"]);
    transcript.insert("polyphen_prediction", Bson::String(polyphen));

    let sift = get_prediction(entry, &["SIFT", "SIFT_PRED"]);
    transcript.insert("sift_prediction", Bson::String(sift));

    if let Some(value) = entry.get("REVEL_RANKSCORE").filter(|v| !v.is_empty())
        && let Some(rankscore) = get_highest_float_score_in_string(value)
    {
        transcript.insert("revel_rankscore", Bson::Double(rankscore));
    }

    if let Some(value) = entry.get("REVEL_SCORE").filter(|v| !v.is_empty())
        && let Some(score) = get_highest_float_score_in_string(value)
    {
        transcript.insert("revel_raw_score", Bson::Double(score));
    }

    parse_transcripts_spliceai(&mut transcript, entry);

    transcript.insert(
        "swiss_prot",
        Bson::String(
            entry
                .get("SWISSPROT")
                .filter(|value| !value.is_empty())
                .unwrap_or("unknown")
                .to_string(),
        ),
    );

    for (source, target) in [
        ("GERP++_RS", "gerp"),
        ("PHASTCONS100WAY_VERTEBRATE", "phast"),
        ("PHYLOP100WAY_VERTEBRATE", "phylop"),
    ] {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            transcript.insert(target, Bson::String(value.to_string()));
        }
    }

    parse_domains(&mut transcript, entry);

    if let Some(coding_sequence) = get_sequence_aux(entry, "HGVSC") {
        transcript.insert("coding_sequence_name", Bson::String(coding_sequence));
    }

    if let Some(protein_sequence) = get_sequence_aux(entry, "HGVSP") {
        transcript.insert("protein_sequence_name", Bson::String(protein_sequence));
    }

    transcript.insert("biotype", entry.get("BIOTYPE").unwrap_or("").to_string());

    for (source, target) in [("EXON", "exon"), ("INTRON", "intron")] {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            transcript.insert(target, Bson::String(value.to_string()));
        }
    }

    if let Some(strand) = get_strand(entry) {
        transcript.insert("strand", Bson::String(strand));
    }

    let functional_annotations: Vec<String> = entry
        .get("CONSEQUENCE")
        .filter(|value| !value.is_empty())
        .map(|value| value.split('&').map(str::to_string).collect())
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
            Bson::Array(region_annotations.into_iter().map(Bson::String).collect()),
        );
    }

    transcript.insert(
        "is_canonical",
        Bson::Boolean(entry.get("CANONICAL") == Some("YES")),
    );

    parse_mane_annotations(&mut transcript, entry);
    parse_cadd(&mut transcript, entry);
    parse_superdups_fracmatch(&mut transcript, entry);
    parse_mt_frequencies(&mut transcript, entry);
    parse_variant_frequencies(&mut transcript, entry);
    parse_clinvar_annotations(&mut transcript, entry);
    parse_dbsnp(&mut transcript, entry);
    parse_cosmic(&mut transcript, entry);

    Some(transcript)
}
