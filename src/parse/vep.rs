use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use rust_htslib::bcf::Record;
use mongodb::bson::{Bson, Document};
use std::collections::HashSet;
use crate::HashMap;
use crate::models::gene::GeneAnnotation;
use crate::models::consequence::SO_TERMS;


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
) -> (Vec<Document>, Vec<GeneAnnotation>) {
    let mut parsed_transcripts = Vec::new();
    let mut gene_annotations = Vec::new();

    let mut dbsnp_ids = HashSet::new();
    let mut cosmic_ids = HashSet::new();

    if !vep_header.is_empty() {
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

                gene_annotations.push(GeneAnnotation {
                    hgnc_id: get_hgnc_id(&raw_transcript),
                    hgnc_symbol: raw_transcript
                        .get("SYMBOL")
                        .cloned(),
                });

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
        }
    }

    if !dbsnp_ids.is_empty() && !variant.contains_key("dbsnp_id") {
        variant.insert(
            "dbsnp_id",
            Bson::String(
                dbsnp_ids.into_iter().collect::<Vec<_>>().join(";"),
            ),
        );
    }

    if !cosmic_ids.is_empty() {
        variant.insert(
            "cosmic_ids",
            Bson::Array(
                cosmic_ids
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }

    (parsed_transcripts, gene_annotations)
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

/// Extract the HGNC identifier from a VEP annotation.
///
/// The HGNC_ID field may contain a prefix (for example, "HGNC:1234").
/// Only the numeric identifier is retained.
///
/// Returns `None` if the HGNC identifier is missing or invalid.
pub fn get_hgnc_id(entry: &HashMap<String, String>) -> Option<String> {
    entry
        .get("HGNC_ID")
        .and_then(|value| value.split(':').last())
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
}

/// Extract a prediction from VEP transcript annotations.
///
/// The prediction fields are typically formatted as
/// `prediction(score)`, for example `deleterious(0.01)`.
/// Returns only the prediction label. If none of the provided fields
/// are available or contain a value, `"unknown"` is returned.
fn get_prediction(entry: &HashMap<String, String>, fields: &[&str]) -> String {
    for field in fields {
        if let Some(value) = entry.get(*field).filter(|value| !value.is_empty()) {
            return value
                .split('(')
                .next()
                .unwrap_or("unknown")
                .to_string();
        }
    }

    "unknown".to_string()
}

/// Return the highest float value from a string with numbers possibly
/// separated by `&`.
///
/// Invalid values are ignored. Returns `None` if the input is empty or
/// no valid float values are found.
fn get_highest_float_score_in_string(value: &str) -> Option<f64> {
    value
        .split('&')
        .filter_map(|part| part.trim().parse::<f64>().ok())
        .max_by(|a, b| a.partial_cmp(b).unwrap())
}


/// Parse SpliceAI annotations from a VEP transcript entry.
///
/// Extracts SpliceAI delta scores and positions from VEP CSQ fields.
/// The maximum delta score is stored together with its corresponding
/// position. Also stores a summary of all splice predictions.
fn parse_transcripts_spliceai(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let spliceai_positions = [
        ("SPLICEAI_PRED_DP_AG", "spliceai_dp_ag"),
        ("SPLICEAI_PRED_DP_AL", "spliceai_dp_al"),
        ("SPLICEAI_PRED_DP_DG", "spliceai_dp_dg"),
        ("SPLICEAI_PRED_DP_DL", "spliceai_dp_dl"),
    ];

    let spliceai_delta_scores = [
        ("SPLICEAI_PRED_DS_AG", "spliceai_ds_ag"),
        ("SPLICEAI_PRED_DS_AL", "spliceai_ds_al"),
        ("SPLICEAI_PRED_DS_DG", "spliceai_ds_dg"),
        ("SPLICEAI_PRED_DS_DL", "spliceai_ds_dl"),
    ];

    for (source, target) in spliceai_positions {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            if let Ok(position) = value.parse::<i32>() {
                transcript.insert(
                    target,
                    Bson::Int32(position),
                );
            }
        }
    }

    for (source, target) in spliceai_delta_scores {
        if let Some(value) = entry.get(source).filter(|v| !v.is_empty()) {
            if let Ok(score) = value.parse::<f64>() {
                transcript.insert(
                    target,
                    Bson::Double(score),
                );
            }
        }
    }

    let spliceai_pairs = [
        ("spliceai_ds_ag", "spliceai_dp_ag"),
        ("spliceai_ds_al", "spliceai_dp_al"),
        ("spliceai_ds_dg", "spliceai_dp_dg"),
        ("spliceai_ds_dl", "spliceai_dp_dl"),
    ];

    let mut max_score: Option<f64> = None;
    let mut max_position: Option<i32> = None;
    let mut predictions = Vec::new();

    for (score_key, position_key) in spliceai_pairs {
        let score = transcript.get_f64(score_key).ok();
        let position = transcript.get_i32(position_key).ok();

        if let Some(score) = score {
            if max_score.map_or(true, |current| score > current) {
                max_score = Some(score);
                max_position = position;
            }
        }

        predictions.push(format!(
            "{} {} {} {}",
            score_key,
            score.map(|x| x.to_string()).unwrap_or("-".to_string()),
            position_key,
            position.map(|x| x.to_string()).unwrap_or("-".to_string()),
        ));
    }

    if let Some(score) = max_score {
        transcript.insert(
            "spliceai_score",
            Bson::Double(score),
        );

        if let Some(position) = max_position {
            transcript.insert(
                "spliceai_position",
                Bson::Int32(position),
            );
        }

        transcript.insert(
            "spliceai_prediction",
            Bson::Array(
                predictions
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }
}

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

/// Extract a sequence annotation from a VEP transcript entry.
///
/// The VEP fields `HGVSC` and `HGVSP` are formatted as
/// `transcript:sequence`. This function returns only the sequence part.
/// Returns `None` if the field is missing or does not contain `:`.
fn get_sequence_aux(entry: &HashMap<String, String>, name: &str) -> Option<String> {
    let sequence_entry = entry
        .get(name)?
        .split(':')
        .collect::<Vec<&str>>();

    if sequence_entry.len() > 1 {
        Some(sequence_entry.last()?.to_string())
    } else {
        None
    }
}

/// Get strand information from a VEP transcript entry.
fn get_strand(entry: &HashMap<String, String>) -> Option<String> {
    match entry.get("STRAND").map(String::as_str) {
        Some("1") => Some("+".to_string()),
        Some("-1") => Some("-".to_string()),
        _ => None,
    }
}

/// Get regional annotations from Sequence Ontology consequence terms.
///
/// Maps functional annotations (SO consequence terms) to their broader
/// functional regions using the predefined Sequence Ontology term mapping.
///
/// Unknown annotations are ignored.
///
/// # Arguments
///
/// * `functional_annotations` - A list of Sequence Ontology consequence terms
///   extracted from a VEP transcript annotation.
///
/// # Returns
///
/// A vector containing the regional annotations corresponding to the provided
/// functional annotations.
fn get_regional_annotation(functional_annotations: &[String]) -> Vec<String> {
    functional_annotations
        .iter()
        .filter_map(|annotation| {
            SO_TERMS
                .get(annotation.as_str())
                .map(|term| term.region.to_string())
        })
        .collect()
}


/// Parse MANE transcript annotations from a VEP transcript entry.
///
/// Extracts MANE Select and MANE Plus Clinical transcript identifiers from
/// VEP v103/MANE v0.92 annotations. Falls back to the older `MANE` field
/// used by previous VEP versions.
fn parse_mane_annotations(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if entry.contains_key("MANE_SELECT") {
        if let Some(mane_select) = entry.get("MANE_SELECT").filter(|v| !v.is_empty()) {
            transcript.insert(
                "mane_select_transcript",
                Bson::String(mane_select.clone()),
            );
        }

        if let Some(mane_plus_clinical) = entry
            .get("MANE_PLUS_CLINICAL")
            .filter(|v| !v.is_empty())
        {
            transcript.insert(
                "mane_plus_clinical_transcript",
                Bson::String(mane_plus_clinical.clone()),
            );
        }
    } else if entry.contains_key("MANE") {
        if let Some(mane) = entry.get("MANE").filter(|v| !v.is_empty()) {
            transcript.insert(
                "mane_select_transcript",
                Bson::String(mane.clone()),
            );
        }
    }
}


/// Parse transcript-level CADD Phred score.
fn parse_cadd(transcript: &mut Document, entry: &HashMap<String, String>) {
    if let Some(cadd_phred) = entry.get("CADD_PHRED") {
        if let Ok(value) = cadd_phred.parse::<f64>() {
            transcript.insert(
                "cadd",
                Bson::Double(value),
            );
        }
    }
}

/// Parse genomic superdups fractional match values from a VEP transcript entry.
fn parse_superdups_fracmatch(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if let Some(superdups_fractmatch) = entry
        .get("GENOMIC_SUPERDUPS_FRAC_MATCH")
        .filter(|value| !value.is_empty())
    {
        let values: Vec<Bson> = superdups_fractmatch
            .split('&')
            .filter_map(|fractmatch| {
                fractmatch.parse::<f64>().ok().map(Bson::Double)
            })
            .collect();

        if !values.is_empty() {
            transcript.insert(
                "superdups_fracmatch",
                Bson::Array(values),
            );
        }
    }
}

/// Parse mitochondrial gnomAD allele frequencies from a VEP transcript.
///
/// Extracts homoplasmic and heteroplasmic mitochondrial allele frequencies
/// when available.
fn parse_mt_frequencies(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if let Some(value) = entry
        .get("GNOMAD_MT_AF_HOM")
        .filter(|value| !value.is_empty())
    {
        if let Ok(value) = value.parse::<f64>() {
            transcript.insert(
                "gnomad_mt_homoplasmic",
                Bson::Double(value),
            );
        }
    }

    if let Some(value) = entry
        .get("GNOMAD_MT_AF_HET")
        .filter(|value| !value.is_empty())
    {
        if let Ok(value) = value.parse::<f64>() {
            transcript.insert(
                "gnomad_mt_heteroplasmic",
                Bson::Double(value),
            );
        }
    }
}

/// Parse variant population frequencies from VEP transcript annotations.
///
/// Supports VEP v90+ frequency fields, including 1000 Genomes, gnomAD,
/// and ExAC frequencies. Stores specific allele frequencies and calculates
/// maximum frequencies across population-specific annotations.
fn parse_variant_frequencies(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let gnomad_csq_keys = [
        "GNOMAD_AF",
        "GNOMADG_AF",
        "GNOMAD_EXOMES_AF",
        "GNOMAD_EXOMES_AF",
    ];

    let thousand_genomes_csq_keys = [
        "AF",
        "1000GAF",
        "1000GP3_AF",
    ];

    let mut thousandg_freqs: Vec<f64> = Vec::new();
    let mut gnomad_freqs: Vec<f64> = Vec::new();

    for (key, value) in entry {
        // Frequency fields end with AF or POPMAX
        if !(key.ends_with("AF") || key.ends_with("POPMAX")) {
            continue;
        }

        if value.is_empty()
            || value == "."
            || value.chars().all(char::is_alphabetic)
        {
            continue;
        }

        let Ok(freq) = value.parse::<f64>() else {
            continue;
        };

        if thousand_genomes_csq_keys.contains(&key.as_str()) {
            transcript.insert(
                "thousand_g_maf",
                Bson::Double(freq),
            );
            continue;
        }

        if gnomad_csq_keys.contains(&key.as_str()) {
            transcript.insert(
                "gnomad_maf",
                Bson::Double(freq),
            );
            continue;
        }

        if key == "EXAC_MAX_AF" {
            transcript.insert(
                "exac_max",
                Bson::Double(freq),
            );
            transcript.insert(
                "exac_maf",
                Bson::Double(freq),
            );
            continue;
        }

        // Remaining population-specific frequencies
        if key.contains("GNOMAD") {
            gnomad_freqs.push(freq);
        } else {
            thousandg_freqs.push(freq);
        }
    }

    if let Some(max) = thousandg_freqs.iter().max_by(|a, b| a.total_cmp(b)) {
        transcript.insert(
            "thousandg_max",
            Bson::Double(*max),
        );
    }

    if let Some(max) = gnomad_freqs.iter().max_by(|a, b| a.total_cmp(b)) {
        transcript.insert(
            "gnomad_max",
            Bson::Double(*max),
        );
    }
}

/// Parse ClinVar annotations from a VEP transcript entry.
///
/// Extracts ClinVar variation identifiers, clinical significance, review
/// status, and clinical significance terms when available.
fn parse_clinvar_annotations(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let clinvar_id = entry
        .get("CLINVAR_CLNVID")
        .or_else(|| entry.get("CLINVAR"));

    if let Some(clinvar_id) = clinvar_id {
        transcript.insert(
            "clinvar_clnvid",
            Bson::String(clinvar_id.clone()),
        );

        if let Some(clnsig) = entry.get("CLINVAR_CLNSIG") {
            transcript.insert(
                "clinvar_clnsig",
                Bson::String(clnsig.to_lowercase()),
            );
        }

        if let Some(revstat) = entry.get("CLINVAR_CLNREVSTAT") {
            transcript.insert(
                "clinvar_revstat",
                Bson::String(revstat.to_lowercase()),
            );
        }
    }

    let clnsig = entry
        .get("CLIN_SIG")
        .or_else(|| entry.get("ClinVar_CLNSIG"));

    if let Some(clnsig) = clnsig.filter(|value| !value.is_empty()) {
        transcript.insert(
            "clnsig",
            Bson::Array(
                clnsig
                    .split('&')
                    .map(|value| Bson::String(value.to_string()))
                    .collect(),
            ),
        );
    }
}

/// Parse dbSNP identifiers from a VEP transcript entry.
///
/// Supports different VEP versions where dbSNP identifiers can be stored
/// in `EXISTING_VARIATION`, `RS_DBSNP150`, or `RS_DBSNP`. Only identifiers
/// starting with `rs` are retained.
fn parse_dbsnp(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let dbsnp_keys = [
        "EXISTING_VARIATION",
        "RS_DBSNP150",
        "RS_DBSNP",
    ];

    let mut dbsnp_ids: HashSet<String> = HashSet::new();

    for key in dbsnp_keys {
        if let Some(variant_ids) = entry.get(key) {
            for variant_id in variant_ids.split('&') {
                if variant_id.starts_with("rs") {
                    dbsnp_ids.insert(variant_id.to_string());
                }
            }
        }
    }

    transcript.insert(
        "dbsnp",
        Bson::Array(
            dbsnp_ids
                .into_iter()
                .map(Bson::String)
                .collect(),
        ),
    );
}

/// Parse COSMIC identifiers from a VEP transcript entry.
///
/// Extracts COSMIC identifiers from `EXISTING_VARIATION` (COSM/COSV
/// prefixes) and from the dedicated `COSMIC` field when available.
fn parse_cosmic(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let mut cosmic_ids: Vec<Bson> = Vec::new();

    if let Some(variant_ids) = entry.get("EXISTING_VARIATION") {
        for variant_id in variant_ids.split('&') {
            if variant_id.starts_with("COSM") || variant_id.starts_with("COSV") {
                cosmic_ids.push(Bson::String(variant_id.to_string()));
            }
        }
    }

    if let Some(cosmic_ids_entry) = entry.get("COSMIC") {
        for cosmic_id in cosmic_ids_entry.split('&') {
            cosmic_ids.push(Bson::String(cosmic_id.to_string()));
        }
    }

    transcript.insert(
        "cosmic",
        Bson::Array(cosmic_ids),
    );
}
