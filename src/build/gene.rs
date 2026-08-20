use crate::build::transcript::build_transcript;
use mongodb::bson::{Bson, Document};
use std::collections::HashMap;

use crate::models::consequence::{CONSEQUENCE, FEATURE_TYPES, SO_TERMS};

/// Build and add complete gene annotations to a variant.
///
/// Only parsed genes containing an HGNC ID are processed. Each gene is
/// enriched with information from the HGNC database and its transcripts
/// are converted to complete transcript objects.
///
/// To avoid uploading excessive amounts of data, processing stops after
/// the 31st gene and the variant is marked as having missing data.
pub fn add_genes(
    variant: &mut Document,
    gene_list: &[Bson],
    hgncid_to_gene: &HashMap<i32, Document>,
) {
    let mut genes = Vec::new();

    for (index, gene) in gene_list.iter().enumerate() {
        let Some(gene) = gene.as_document() else {
            continue;
        };

        if !gene.contains_key("hgnc_id") {
            continue;
        }

        let gene_obj = build_gene(gene, hgncid_to_gene);
        genes.push(Bson::Document(gene_obj));

        if index > 30 {
            variant.insert("missing_data", true);
            break;
        }
    }

    if !genes.is_empty() {
        variant.insert("genes", Bson::Array(genes));
    }
}

/// Copy an optional field from a parsed gene document to the built gene
/// document.
///
/// The field is only added when the source document contains the requested
/// key. This is used for optional gene annotations such as SpliceAI, HGVS,
/// canonical transcript, and exon information.
fn insert_optional_gene_field(
    gene_obj: &mut Document,
    gene: &Document,
    source_key: &str,
    target_key: &str,
) {
    if let Some(value) = gene.get(source_key) {
        gene_obj.insert(target_key, value.clone());
    }
}

/// Build a complete gene annotation from parsed VCF information.
///
/// The gene is enriched with information from the HGNC database, including
/// symbol, Ensembl ID, description, inheritance models, and phenotypes.
/// Transcript annotations and variant-level gene annotations are also added
/// when available. If the gene is not found in the HGNC database, the
/// annotations available from the VCF are still returned.
pub fn build_gene(gene: &Document, hgncid_to_gene: &HashMap<i32, Document>) -> Document {
    let mut gene_obj = Document::new();

    let Some(hgnc_id) = gene.get_i32("hgnc_id").ok() else {
        return gene_obj;
    };

    gene_obj.insert("hgnc_id", hgnc_id);

    // Get gene information from the database, if available.
    if let Some(hgnc_gene) = hgncid_to_gene.get(&hgnc_id) {
        if let Ok(value) = hgnc_gene.get_str("hgnc_symbol") {
            gene_obj.insert("hgnc_symbol", value);
        }

        if let Ok(value) = hgnc_gene.get_str("ensembl_id") {
            gene_obj.insert("ensembl_id", value);
        }

        if let Ok(value) = hgnc_gene.get_str("description") {
            gene_obj.insert("description", value);
        }

        if let Ok(value) = hgnc_gene.get_array("inheritance_models")
            && !value.is_empty()
        {
            gene_obj.insert("inheritance", value.clone());
        }

        if let Ok(value) = hgnc_gene.get_array("phenotypes")
            && !value.is_empty()
        {
            gene_obj.insert("phenotypes", value.clone());
        }
    }

    // Build transcripts.
    let transcripts = gene
        .get_array("transcripts")
        .ok()
        .map(|transcripts| {
            transcripts
                .iter()
                .filter_map(Bson::as_document)
                .map(build_transcript)
                .map(Bson::Document)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();

    gene_obj.insert("transcripts", Bson::Array(transcripts));

    // Functional annotation.
    if let Ok(value) = gene.get_str("most_severe_consequence") {
        if SO_TERMS.contains_key(value) {
            gene_obj.insert("functional_annotation", value);
        } else {
            eprintln!("Invalid functional annotation {}", value);
        }
    }

    // Region annotation.
    if let Ok(value) = gene.get_str("region_annotation") {
        if FEATURE_TYPES.contains(&value) {
            gene_obj.insert("region_annotation", value);
        } else {
            eprintln!("Invalid region annotation {}", value);
        }
    }

    // SIFT prediction.
    if let Ok(value) = gene.get_str("most_severe_sift") {
        if CONSEQUENCE.contains(&value) {
            gene_obj.insert("sift_prediction", value);
        } else {
            eprintln!("Invalid sift prediction {}", value);
        }
    }

    // PolyPhen prediction.
    if let Ok(value) = gene.get_str("most_severe_polyphen") {
        if CONSEQUENCE.contains(&value) {
            gene_obj.insert("polyphen_prediction", value);
        } else {
            eprintln!("Invalid polyphen prediction {}", value);
        }
    }

    insert_optional_gene_field(
        &mut gene_obj,
        gene,
        "most_severe_spliceai_score",
        "spliceai_score",
    );
    insert_optional_gene_field(
        &mut gene_obj,
        gene,
        "most_severe_spliceai_position",
        "spliceai_position",
    );
    insert_optional_gene_field(
        &mut gene_obj,
        gene,
        "spliceai_prediction",
        "spliceai_prediction",
    );
    insert_optional_gene_field(&mut gene_obj, gene, "hgvs_identifier", "hgvs_identifier");
    insert_optional_gene_field(
        &mut gene_obj,
        gene,
        "canonical_transcript",
        "canonical_transcript",
    );
    insert_optional_gene_field(&mut gene_obj, gene, "exon", "exon");

    gene_obj
}
