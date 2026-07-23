use std::collections::HashMap;

use mongodb::bson::{Bson, Document};
use crate::models::gene::GeneAnnotation;

/// Group transcripts into genes.
///
/// Genes are identified using HGNC ID when available, otherwise HGNC symbol.
/// Transcript documents do not contain gene information; gene annotations are
/// provided separately.
///
/// Returns a list of gene documents containing their transcripts.
pub fn parse_genes(
    transcripts: Vec<Document>,
    annotations: Vec<GeneAnnotation>,
) -> Vec<Document> {
    let mut genes: HashMap<String, Document> = HashMap::new();

    for (transcript, annotation) in transcripts.into_iter().zip(annotations) {
        let identifier = annotation
            .hgnc_id
            .clone()
            .or(annotation.hgnc_symbol.clone());

        let Some(identifier) = identifier else {
            continue;
        };

        let gene = genes
            .entry(identifier.clone())
            .or_insert_with(|| {
                let mut doc = Document::new();

                doc.insert(
                    "hgnc_id",
                    annotation.hgnc_id.clone(),
                );

                doc.insert(
                    "hgnc_symbol",
                    annotation.hgnc_symbol.clone(),
                );

                doc.insert(
                    "transcripts",
                    Bson::Array(Vec::new()),
                );

                doc
            });

        if let Ok(transcripts) = gene.get_array_mut("transcripts") {
            transcripts.push(Bson::Document(transcript));
        }
    }

    genes.into_values().collect()
}

