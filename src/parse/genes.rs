use std::collections::HashMap;

use mongodb::bson::{Bson, Document};

use crate::models::consequence::SO_TERMS;

/// Parse transcript information and build gene annotations.
///
/// Transcripts are grouped by HGNC ID when available, otherwise by HGNC
/// symbol. Transcripts without a gene identifier are skipped.
///
/// For each gene, the transcript with the most severe consequence is selected
/// according to the Sequence Ontology consequence ranking.
///
/// Args:
///     transcripts: Parsed VEP transcript annotations.
///
/// Returns:
///     A list of gene documents containing transcripts and the most severe
///     annotation information.
pub fn parse_genes(transcripts: &[Document]) -> Vec<Document> {
    let mut genes_to_transcripts: HashMap<String, Vec<Document>> = HashMap::new();

    // Group transcripts by gene
    for transcript in transcripts {
        let hgnc_id = transcript
            .get_str("hgnc_id")
            .ok();

        let hgnc_symbol = transcript
            .get_str("hgnc_symbol")
            .ok();

        let gene_identifier = hgnc_id.or(hgnc_symbol);

        let Some(gene_identifier) = gene_identifier else {
            continue;
        };

        genes_to_transcripts
            .entry(gene_identifier.to_string())
            .or_default()
            .push(transcript.clone());
    }

    let mut genes = Vec::new();

    for (_, gene_transcripts) in genes_to_transcripts {
        let mut most_severe_rank = u32::MAX;
        let mut most_severe_consequence: Option<String> = None;
        let mut most_severe_transcript: Option<Document> = None;
        let mut most_severe_region: Option<String> = None;

        let mut most_severe_sift = None;
        let mut most_severe_polyphen = None;
        let mut most_severe_spliceai_score = None;
        let mut most_severe_spliceai_position = None;
        let mut spliceai_prediction = None;

        let mut hgvs_identifier = None;
        let mut exon = None;
        let mut canonical_transcript = None;

        let mut hgnc_id = None;
        let mut hgnc_symbol = None;

        for transcript in &gene_transcripts {
            hgnc_id = transcript.get("hgnc_id").cloned();
            hgnc_symbol = transcript.get("hgnc_symbol").cloned();

            if hgvs_identifier.is_none() {
                hgvs_identifier = transcript.get("coding_sequence_name").cloned();
            }

            if exon.is_none() {
                exon = transcript.get("exon").cloned();
            }

            if let Some(Bson::Array(consequences)) =
                transcript.get("functional_annotations")
            {
                for consequence in consequences {
                    let Some(consequence) = consequence.as_str() else {
                        continue;
                    };

                    let Some(so_term) = SO_TERMS.get(consequence) else {
                        continue;
                    };

                    if so_term.rank > most_severe_rank {
                        continue;
                    }

                    most_severe_rank = so_term.rank;
                    most_severe_consequence =
                        Some(consequence.to_string());

                    most_severe_transcript =
                        Some(transcript.clone());

                    most_severe_region =
                        Some(so_term.region.to_string());

                    most_severe_sift =
                        transcript.get("sift_prediction").cloned();

                    most_severe_polyphen =
                        transcript.get("polyphen_prediction").cloned();

                    most_severe_spliceai_score =
                        transcript.get("spliceai_delta_score").cloned();

                    most_severe_spliceai_position =
                        transcript.get("spliceai_delta_position").cloned();

                    spliceai_prediction =
                        transcript.get("spliceai_prediction").cloned();
                }
            }

            if transcript
                .get_bool("is_canonical")
                .unwrap_or(false)
            {
                canonical_transcript =
                    transcript.get("transcript_id").cloned();

                if transcript.get("coding_sequence_name").is_some() {
                    hgvs_identifier =
                        transcript.get("coding_sequence_name").cloned();

                    exon =
                        transcript.get("exon").cloned();
                }
            }
        }

        let mut gene = Document::new();

        gene.insert(
            "transcripts",
            Bson::Array(
                gene_transcripts
                    .into_iter()
                    .map(Bson::Document)
                    .collect(),
            ),
        );

        gene.insert(
            "most_severe_transcript",
            most_severe_transcript
                .map(Bson::Document)
                .unwrap_or(Bson::Null),
        );

        gene.insert(
            "most_severe_consequence",
            most_severe_consequence
                .map(Bson::String)
                .unwrap_or(Bson::Null),
        );

        gene.insert(
            "most_severe_sift",
            most_severe_sift.unwrap_or(Bson::Null),
        );

        gene.insert(
            "most_severe_polyphen",
            most_severe_polyphen.unwrap_or(Bson::Null),
        );

        gene.insert(
            "most_severe_spliceai_score",
            most_severe_spliceai_score.unwrap_or(Bson::Null),
        );

        gene.insert(
            "most_severe_spliceai_position",
            most_severe_spliceai_position.unwrap_or(Bson::Null),
        );

        gene.insert(
            "spliceai_prediction",
            spliceai_prediction.unwrap_or(Bson::Null),
        );

        gene.insert(
            "hgnc_id",
            hgnc_id.unwrap_or(Bson::Null),
        );

        gene.insert(
            "hgnc_symbol",
            hgnc_symbol.unwrap_or(Bson::Null),
        );

        gene.insert(
            "region_annotation",
            most_severe_region
                .map(Bson::String)
                .unwrap_or(Bson::Null),
        );

        gene.insert(
            "hgvs_identifier",
            hgvs_identifier.unwrap_or(Bson::Null),
        );

        gene.insert(
            "canonical_transcript",
            canonical_transcript.unwrap_or(Bson::Null),
        );

        gene.insert(
            "exon",
            exon.unwrap_or(Bson::Null),
        );

        genes.push(gene);
    }

    genes
}

