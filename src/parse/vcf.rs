use rust_htslib::bcf::{Read, Reader};
use mongodb::bson::{doc, self};
use std::collections::HashMap;
use mongodb::bson::{Bson, Document};
use crate::parse::coordinates::parse_coordinates;
use crate::parse::alleles::parse_alleles;
use crate::parse::filters::parse_filters;
use crate::parse::compounds::parse_compounds;
use crate::parse::ids::parse_ids;
use crate::parse::rank_scores::parse_rank_scores;
use crate::parse::genetic_models::parse_genetic_models;
use crate::parse::info::{parse_info_int, parse_info_string, parse_custom_data};
use crate::parse::strs::set_str_info;
use crate::parse::meis::set_mei_info;
use crate::parse::fusions::set_fusion_info;
use crate::parse::genotypes::{parse_genotypes, validate_sample_mapping};
use crate::parse::mt_annotations::{set_mitomap_associated_diseases, set_hmtvar};
use crate::parse::vep::{parse_vep_header, parse_vep_transcripts};
use crate::parse::genes::parse_genes;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;
use crate::models::cytoband::Cytoband;
use crate::models::sample::SampleInfo;


/// Processes a VCF file and parses each record according to the variant category.
///
/// The function reads the VCF file at the provided path, iterates through all
/// variant records, and dispatches parsing based on the specified category
/// (for example, "snv", "sv", or "str").
///
/// # Arguments
///
/// * `path` - Path to the input VCF file.
/// * `category` - Variant category used to select the appropriate parser.
/// * `variant_type` - Variant type (clinical or research).
/// * `case_id` - _id of a case.
/// * `cytobands` - A list of parsed cytobands, reflecting the case genome build
/// * `sample_mapping`, a list of samples, with ID and expected position on the VCF
///
/// # Panics
///
/// Panics if the VCF file cannot be opened or if a record cannot be read.
pub fn process_vcf(path: &str, category: VariantCategory, variant_type: VariantType, case_id: &str, cytobands: &HashMap<String, Vec<Cytoband>>, sample_mapping: &HashMap<String, SampleInfo>) {

    let mut vcf = Reader::from_path(path)
        .expect("couldn't open input vcf");

    let header = vcf.header().clone();
    let vep_header = parse_vep_header(&header);

    if let Err(error) = validate_sample_mapping(vcf.header(), &sample_mapping) {
        eprintln!("Sample mapping validation failed: {}", error);
        return;
    }

    for result in vcf.records() {
        let record = result.unwrap();
        let case_id = case_id.to_string();
        let coordinates = parse_coordinates(&record, &header, cytobands, &category);
        let variant_type = variant_type.to_string();
        let (reference, alternative) = parse_alleles(&record, category);
        let ids = parse_ids(&coordinates.chromosome, &coordinates.position, &reference, &alternative, &case_id, &variant_type);
        let filters = parse_filters(&record, &header);
        let compound_info = record
            .info(b"Compounds")
            .string()
            .ok()
            .flatten()
            .and_then(|values| values.first().map(|value| {
                String::from_utf8_lossy(value).to_string()
            }));
        let compounds = parse_compounds(compound_info, &case_id, &variant_type);
        let compounds_bson = bson::to_bson(&compounds).expect("Failed to convert compounds to BSON");
        let (rank_score, norm_rank_score) = parse_rank_scores(&record, &case_id);
        let genetic_models = parse_genetic_models(&record, &case_id);
        let samples = parse_genotypes(&record, &sample_mapping, category);

        // This structure contains fields common to all variants categories
        let mut variant = doc! {
            "simple_id": ids.simple_id,
            "variant_id": ids.variant_id,
            "display_name": ids.display_name,
            "document_id": ids.document_id,
            "case_id": case_id,

            "compounds": compounds_bson,

            "rank_score": rank_score,
            "norm_rank_score": norm_rank_score,

            "type": variant_type,

            "chromosome": coordinates.chromosome,
            "end_chrom": coordinates.end_chrom,
            "position": coordinates.position as i64,
            "end": coordinates.end as i64,
            "length": coordinates.length,

            "category": category.to_string(),
            "sub_category": coordinates.sub_category,

            "reference": reference,
            "alternative": alternative,

            "cytoband_start": coordinates.cytoband_start,
            "cytoband_end": coordinates.cytoband_end,

            "filters": filters,
            "quality": record.qual(),

            "genetic_models": genetic_models,

            "samples": samples,
        };

        if coordinates.mate_id.is_some(){
            variant.insert("mate_id", coordinates.mate_id);
        }

        let azlength = parse_info_string(&record, b"AZLENGTH").and_then(|value| value.parse::<i32>().ok());
        if let Some(value) = azlength {
            variant.insert("azlength", value);
        }

        let azqual = parse_info_string(&record, b"AZQUAL").and_then(|value| value.parse::<f64>().ok());
        if let Some(value) = azqual {
            variant.insert("azqual", value);
        }

        if let Some(custom) = parse_custom_data(parse_info_string(&record, b"SCOUT_CUSTOM")) {
                variant.insert("custom", custom);
            }

        let id = record.id();
        let variant_id = String::from_utf8_lossy(&id);

        if variant_id.contains("rs") {
            variant.insert("dbsnp_id", variant_id.to_string());
        }

        set_mitomap_associated_diseases(&record, &mut variant);
        set_hmtvar(&record, &mut variant);

        match category {
            VariantCategory::Str => {
                set_str_info(&record, &mut variant);
            }

            VariantCategory::Mei => {
                set_mei_info(&record, &mut variant);
            }

            VariantCategory::Fusion => {
                set_fusion_info(&record, &mut variant);
                return; // Setting of genes and transcripts is handled specifically by set_fusion_info for this category
            }

            VariantCategory::Cancer | VariantCategory::CancerSv => {
                if let Some(value) = parse_info_int(&record, b"SOMATICSCORE") {
                    variant.insert(
                        "somatic_score",
                        bson::Bson::Int32(value),
                    );
                }
            }

            _ => {}
        }
        let (parsed_transcripts, gene_annotations) = parse_vep_transcripts(&record, &vep_header, &mut variant);
        let genes = parse_genes(parsed_transcripts, gene_annotations);
        variant.insert(
            "genes",
            Bson::Array(
                genes
                    .into_iter()
                    .map(Bson::Document)
                    .collect(),
            ),
        );

        println!("{:#?}\n\n", variant["genes"]); 
            
    }

}