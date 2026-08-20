use mongodb::bson::{Bson, Document};

const BUILD_TRANSCRIPT_OPTIONAL_KEYS: &[&str] = &[
    "protein_id",
    "sift_prediction",
    "polyphen_prediction",
    "swiss_prot",
    "pfam_domain",
    "prosite_profile",
    "smart_domain",
    "biotype",
    "functional_annotations",
    "region_annotations",
    "exon",
    "intron",
    "strand",
    "coding_sequence_name",
    "protein_sequence_name",
    "superdups_fracmatch",
    "mane_select_transcript",
    "mane_plus_clinical_transcript",
];

/// Build a transcript annotation from parsed VCF information.
///
/// The transcript contains the transcript ID and HGNC ID together with
/// optional annotations present in the parsed transcript. The canonical
/// status is always included and defaults to `false`.
///
/// These annotations represent transcripts parsed from the VCF rather than
/// transcript definitions collected from Ensembl.
pub fn build_transcript(transcript: &Document) -> Document {
    let mut transcript_obj = Document::new();

    let transcript_id = transcript
        .get_str("transcript_id")
        .expect("Transcript must have a transcript_id");

    transcript_obj.insert("transcript_id", transcript_id);

    let hgnc_id = transcript
        .get("hgnc_id")
        .expect("Transcript must have an hgnc_id");

    transcript_obj.insert("hgnc_id", hgnc_id.clone());

    for key in BUILD_TRANSCRIPT_OPTIONAL_KEYS {
        if let Some(value) = transcript.get(*key) {
            if !matches!(value, Bson::Null) {
                transcript_obj.insert(*key, value.clone());
            }
        }
    }

    let is_canonical = transcript.get_bool("is_canonical").unwrap_or(false);

    transcript_obj.insert("is_canonical", is_canonical);

    transcript_obj
}
