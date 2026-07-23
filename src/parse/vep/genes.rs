use std::collections::HashMap;
use std::collections::HashSet;
use mongodb::bson::{doc, Bson, Document};

use crate::models::consequence::SO_TERMS;

/// Parse gene information from transcript annotations.
///
/// Transcripts are grouped by gene using the HGNC identifier when available,
/// falling back to the HGNC symbol if the identifier is missing.
///
/// For each gene, all associated transcripts are stored and the transcript
/// with the most severe consequence is selected according to the SO_TERMS
/// consequence ranking. Gene-level annotations such as the most severe
/// consequence, region, SIFT/PolyPhen predictions, SpliceAI information,
/// canonical transcript and HGVS identifier are extracted from the selected
/// transcript.
///
/// Transcripts without a valid gene identifier are skipped.
///
/// # Arguments
///
/// * `transcripts` - A slice of parsed VEP transcript BSON documents.
///
/// # Returns
///
/// A vector of BSON documents, where each document represents a gene and
/// contains its transcripts and gene-level annotations.
pub fn parse_genes(transcripts: &[Document]) -> Vec<Document> {
    let mut genes_to_transcripts: HashMap<String, Vec<Document>> = HashMap::new();

    // Group transcripts by gene
    for (idx, transcript) in transcripts.iter().enumerate() {
        let hgnc_id = match transcript.get("hgnc_id") {
            Some(Bson::String(value)) if !value.is_empty() => Some(value.clone()),
            Some(Bson::Int32(value)) => Some(value.to_string()),
            Some(Bson::Int64(value)) => Some(value.to_string()),
            _ => None,
        };

        let hgnc_symbol = transcript
            .get_str("hgnc_symbol")
            .ok()
            .filter(|value| !value.is_empty())
            .map(str::to_string);

        let gene_identifier = hgnc_id.or(hgnc_symbol);

        let Some(gene_identifier) = gene_identifier else {
            println!("Skipping transcript {}: no gene identifier", idx);
            continue;
        };

        genes_to_transcripts
            .entry(gene_identifier)
            .or_default()
            .push(transcript.clone());
    }

    let mut genes = Vec::new();

    for (_gene_id, gene_transcripts) in genes_to_transcripts {
        let mut most_severe_rank = u32::MAX;

        let mut most_severe_consequence: Option<Bson> = None;
        let mut most_severe_region: Option<Bson> = None;

        let mut most_severe_sift: Option<Bson> = None;
        let mut most_severe_polyphen: Option<Bson> = None;

        let mut most_severe_spliceai_score: Option<Bson> = None;
        let mut most_severe_spliceai_position: Option<Bson> = None;
        let mut spliceai_prediction: Option<Bson> = None;

        let mut hgvs_identifier: Option<Bson> = None;
        let mut exon: Option<Bson> = None;
        let mut canonical_transcript: Option<Bson> = None;

        let mut hgnc_id: Option<Bson> = None;
        let mut hgnc_symbol: Option<Bson> = None;

        for transcript in &gene_transcripts {
            if hgnc_id.is_none() {
                hgnc_id = transcript.get("hgnc_id").cloned();
            }

            if hgnc_symbol.is_none() {
                hgnc_symbol = transcript.get("hgnc_symbol").cloned();
            }

            if hgvs_identifier.is_none() {
                hgvs_identifier = transcript.get("coding_sequence_name").cloned();
            }

            if exon.is_none() {
                exon = transcript.get("exon").cloned();
            }

            if let Ok(consequences) = transcript.get_array("functional_annotations") {
                for consequence in consequences {
                    let Some(consequence) = consequence.as_str() else {
                        continue;
                    };

                    let Some(so_term) = SO_TERMS.get(consequence) else {
                        println!("Unknown consequence: {}", consequence);
                        continue;
                    };

                    if so_term.rank > most_severe_rank {
                        continue;
                    }

                    most_severe_rank = so_term.rank;

                    most_severe_consequence =
                        Some(Bson::String(consequence.to_string()));

                    most_severe_region =
                        Some(Bson::String(so_term.region.to_string()));

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

                    exon = transcript.get("exon").cloned();
                }
            }
        }

        genes.push(doc! {
            "transcripts": Bson::Array(
                gene_transcripts
                    .into_iter()
                    .map(Bson::Document)
                    .collect()
            ),
            "most_severe_consequence": most_severe_consequence,
            "most_severe_sift": most_severe_sift,
            "most_severe_polyphen": most_severe_polyphen,
            "most_severe_spliceai_score": most_severe_spliceai_score,
            "most_severe_spliceai_position": most_severe_spliceai_position,
            "spliceai_prediction": spliceai_prediction,
            "hgnc_id": hgnc_id,
            "hgnc_symbol": hgnc_symbol,
            "region_annotation": most_severe_region,
            "hgvs_identifier": hgvs_identifier,
            "canonical_transcript": canonical_transcript,
            "exon": exon,
        });
    }

    println!("parse_genes: returning {} genes", genes.len());

    genes
}


/// Collect HGNC identifiers from parsed genes and variant annotations.
///
/// HGNC identifiers are collected from the gene annotations. For STR variants,
/// Stranger can annotate the HGNC identifier directly in the INFO field
/// (`HGNCId`), which is added as an additional identifier. If no genes were
/// parsed, a minimal gene document containing only the HGNC identifier is
/// created.
///
/// The collected identifiers are stored in the variant as `hgnc_ids`.
pub fn set_hgnc_ids(variant: &mut Document) {
    let mut hgnc_ids: HashSet<String> = HashSet::new();

    // Collect HGNC IDs from parsed genes
    if let Some(Bson::Array(genes)) = variant.get("genes") {
        for gene in genes {
            let Some(gene) = gene.as_document() else {
                continue;
            };

            if let Some(Bson::Int64(hgnc_id)) = gene.get("hgnc_id") {
                hgnc_ids.insert(hgnc_id.to_string());
            }

            if let Some(Bson::Int32(hgnc_id)) = gene.get("hgnc_id") {
                hgnc_ids.insert(hgnc_id.to_string());
            }

            if let Some(Bson::String(hgnc_id)) = gene.get("hgnc_id") {
                hgnc_ids.insert(hgnc_id.clone());
            }
        }
    }

    // STR HGNC IDs are annotated by Stranger
    if let Some(Bson::String(str_hgnc_id)) = variant.get("HGNCId") {
        hgnc_ids.insert(str_hgnc_id.clone());

        let has_genes = matches!(
            variant.get("genes"),
            Some(Bson::Array(genes)) if !genes.is_empty()
        );

        if !has_genes {
            variant.insert(
                "genes",
                Bson::Array(vec![
                    Bson::Document(doc! {
                        "hgnc_id": str_hgnc_id,
                    })
                ]),
            );
        }
    }

    if !hgnc_ids.is_empty() {
        variant.insert(
            "hgnc_ids",
            Bson::Array(
                hgnc_ids
                    .into_iter()
                    .map(Bson::String)
                    .collect(),
            ),
        );
    }
}