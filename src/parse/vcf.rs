use rust_htslib::bcf::{Read, Reader};
use mongodb::bson::{doc, self};
use std::collections::HashMap;
use crate::parse::coordinates::parse_coordinates;
use crate::parse::alleles::parse_alleles;
use crate::parse::filters::parse_filters;
use crate::parse::compounds::parse_compounds;
use crate::parse::ids::parse_ids;
use crate::parse::rank_scores::parse_rank_scores;
use crate::parse::genetic_models::parse_genetic_models;
use crate::parse::info::parse_info_int;
use crate::parse::strs::set_str_info;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;
use crate::models::cytoband::Cytoband;


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
///
/// # Panics
///
/// Panics if the VCF file cannot be opened or if a record cannot be read.
pub fn process_vcf(path: &str, category: VariantCategory, variant_type: VariantType, case_id: &str, cytobands: &HashMap<String, Vec<Cytoband>>) {

    let mut vcf = Reader::from_path(path)
        .expect("couldn't open input vcf");

    let header = vcf.header().clone();

    for result in vcf.records() {
        let case_id = case_id.to_string();
        let variant_type = variant_type.to_string();
        let record = result.unwrap();
        let coordinates = parse_coordinates(&record, &header, cytobands);
        let (reference, alternative) = parse_alleles(&record, category);
        let ids = parse_ids(&coordinates.chromosome, &coordinates.position, &reference, &alternative, &case_id, &variant_type,);
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

            "mate_id": coordinates.mate_id,
            "cytoband_start": coordinates.cytoband_start,
            "cytoband_end": coordinates.cytoband_end,

            "filters": filters,
            "quality": record.qual(),

            "genetic_models": genetic_models,
        };

        match category {
            VariantCategory::Str => {
                set_str_info(&record, &mut variant);
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

        println!("{:#?}", variant);
            
    }

}