use crate::models::case::SampleConfig;
use crate::models::cytoband::Cytoband;
use crate::models::sample::SampleInfo;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;
use crate::parse::alleles::parse_alleles;
use crate::parse::callers::parse_callers;
use crate::parse::compounds::parse_compounds;
use crate::parse::conservations::parse_conservations;
use crate::parse::coordinates::parse_coordinates;
use crate::parse::filters::parse_filters;
use crate::parse::frequencies::{
    add_frequencies, parse_frequencies, parse_mei_frequencies, parse_sv_frequencies,
};
use crate::parse::fusions::set_fusion_info;
use crate::parse::genetic_models::parse_genetic_models;
use crate::parse::genotypes::{parse_genotypes, validate_sample_mapping};
use crate::parse::header::{
    parse_local_archive_header, parse_rank_results_header, parse_vep_header,
};
use crate::parse::ids::parse_ids;
use crate::parse::info::{parse_custom_data, parse_info_int, parse_info_string};
use crate::parse::loqusdb_frequencies::add_loqus_archive_frequencies;
use crate::parse::meis::set_mei_info;
use crate::parse::mt_annotations::{set_hmtvar, set_mitomap_associated_diseases};
use crate::parse::onco_clnsig::parse_clnsig_onc;
use crate::parse::rank_scores::{parse_rank_result, parse_rank_score_other, parse_rank_scores};
use crate::parse::severity::set_severity_predictions;
use crate::parse::strs::set_str_info;
use crate::parse::vep::clnsig::{build_clnsig, parse_clnsig};
use crate::parse::vep::genes::{parse_genes, set_hgnc_ids};
use crate::parse::vep::transcripts::parse_vep_transcripts;
use mongodb::bson::{self, Bson, Document, doc};
use rust_htslib::bcf::{Read, Reader};
use std::collections::{HashMap, HashSet};

/// Builds a mapping between configured samples and their positions in a VCF.
///
/// Sample identity is taken from the case configuration, while the sample
/// index is determined from the VCF header and may therefore differ between
/// VCF files.
///
/// Returns an error if a configured sample cannot be found in the VCF.
pub fn parse_sample_mapping(
    samples: &[SampleConfig],
    vcf_samples: &[String],
) -> Result<HashMap<String, SampleInfo>, Box<dyn std::error::Error>> {
    let mut sample_mapping = HashMap::new();

    for sample in samples {
        let sample_index = vcf_samples
            .iter()
            .position(|vcf_sample| vcf_sample == &sample.sample_id)
            .ok_or_else(|| format!("Sample {} not found in VCF", sample.sample_id))?;

        sample_mapping.insert(
            sample.sample_id.clone(),
            SampleInfo {
                display_name: sample
                    .sample_name
                    .clone()
                    .unwrap_or_else(|| sample.sample_id.clone()),
                vcf_index: sample_index,
            },
        );
    }

    Ok(sample_mapping)
}

/// Adds gene panel information to a parsed variant.
///
/// Collects all gene panels associated with the variant's HGNC IDs and adds
/// the unique panel names to the `panels` field when any are found.
pub fn link_gene_panels(variant: &mut Document, gene_to_panels: &HashMap<i32, HashSet<String>>) {
    let Some(hgnc_ids) = variant.get_array("hgnc_ids").ok() else {
        return;
    };

    let mut panel_names = HashSet::new();

    for hgnc_id in hgnc_ids {
        if let Some(hgnc_id) = hgnc_id.as_i32() {
            if let Some(gene_panels) = gene_to_panels.get(&hgnc_id) {
                panel_names.extend(gene_panels.iter().cloned());
            }
        }
    }

    if !panel_names.is_empty() {
        variant.insert("panels", panel_names.into_iter().collect::<Vec<_>>());
    }
}

/// Processes a VCF file and parses each record according to the variant category.
///
/// The function reads the VCF file at the provided path, determines the sample
/// order from the VCF header, and builds a sample mapping using the samples
/// configured for the case. The mapping is created separately for each VCF
/// because sample order may differ between VCF files.
///
/// # Arguments
///
/// * `path` - Path to the input VCF file.
/// * `category` - Variant category used to select the appropriate parser.
/// * `variant_type` - Variant type for the VCF.
/// * `case_id` - ID of the case.
/// * `cytobands` - Parsed cytobands corresponding to the case genome build.
/// * `samples` - Samples configured for the case.
///
/// # Panics
///
/// Panics if the VCF file cannot be opened, if the sample mapping cannot be
/// created, or if a record cannot be read.
pub fn process_vcf(
    path: &str,
    category: VariantCategory,
    variant_type: VariantType,
    case_id: &str,
    cytobands: &HashMap<String, Vec<Cytoband>>,
    samples: &[SampleConfig],
    gene_to_panels: &HashMap<i32, HashSet<String>>,
) {
    let mut vcf = Reader::from_path(path).expect("couldn't open input vcf");

    let header = vcf.header().clone();

    // Sample order can differ between VCF files, so build the mapping from
    // the header of the current VCF.
    let vcf_samples: Vec<String> = header
        .samples()
        .iter()
        .map(|sample| String::from_utf8_lossy(sample).to_string())
        .collect();

    let sample_mapping =
        parse_sample_mapping(samples, &vcf_samples).expect("Failed to build sample mapping");

    let vep_header = parse_vep_header(&header);
    let rank_results_header = parse_rank_results_header(&header);

    let local_archive_info = parse_local_archive_header(path);

    if let Err(error) = validate_sample_mapping(vcf.header(), &sample_mapping) {
        eprintln!("Sample mapping validation failed: {}", error);
        return;
    }

    let mut variant_count = 0;

    for result in vcf.records() {
        let record = result.unwrap();

        let coordinates = parse_coordinates(&record, &header, cytobands, &category);

        let variant_type = variant_type.to_string();

        let (reference, alternative) = parse_alleles(&record, category);

        let ids = parse_ids(
            &coordinates.chromosome,
            &coordinates.position,
            &reference,
            &alternative,
            case_id,
            &variant_type,
        );

        let filters = parse_filters(&record, &header);
        let callers = parse_callers(&record, category, &filters);

        let compound_info = record
            .info(b"Compounds")
            .string()
            .ok()
            .flatten()
            .and_then(|values| {
                values
                    .first()
                    .map(|value| String::from_utf8_lossy(value).to_string())
            });

        let compounds = parse_compounds(compound_info, case_id, &variant_type);

        let compounds_bson =
            bson::to_bson(&compounds).expect("Failed to convert compounds to BSON");

        let (rank_score, norm_rank_score) = parse_rank_scores(&record, case_id);

        let genetic_models = parse_genetic_models(&record, case_id);

        let samples = parse_genotypes(&record, &sample_mapping, category);

        // This structure contains fields common to all variant categories.
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

        if coordinates.mate_id.is_some() {
            variant.insert("mate_id", coordinates.mate_id);
        }

        let azlength =
            parse_info_string(&record, b"AZLENGTH").and_then(|value| value.parse::<i32>().ok());

        if let Some(value) = azlength {
            variant.insert("azlength", value);
        }

        let azqual =
            parse_info_string(&record, b"AZQUAL").and_then(|value| value.parse::<f64>().ok());

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

        let mut frequencies = bson::Document::new();

        match category {
            VariantCategory::Snv => {
                if let Some(rank_score_other) = parse_rank_score_other(&record) {
                    variant.insert("rank_score_other", rank_score_other);
                }
            }

            VariantCategory::Str => {
                set_str_info(&record, &mut variant);
            }

            VariantCategory::Mei => {
                set_mei_info(&record, &mut variant);
                frequencies.extend(parse_mei_frequencies(&record));
            }

            VariantCategory::Fusion => {
                set_fusion_info(&record, &mut variant);
                return;
            }

            VariantCategory::Cancer | VariantCategory::CancerSv => {
                if let Some(value) = parse_info_int(&record, b"SOMATICSCORE") {
                    variant.insert("somatic_score", bson::Bson::Int32(value));
                }

                if parse_info_string(&record, b"MSK_MVL").is_some() {
                    variant.insert("mvl_tag", true);
                }
            }

            _ => {}
        }

        if matches!(category, VariantCategory::Sv | VariantCategory::CancerSv) {
            frequencies.extend(parse_sv_frequencies(&record));
        }

        let parsed_transcripts = parse_vep_transcripts(&record, &vep_header, &mut variant);

        let genes = parse_genes(&parsed_transcripts);

        variant.insert(
            "genes",
            Bson::Array(genes.into_iter().map(Bson::Document).collect()),
        );

        set_hgnc_ids(&mut variant);

        let clnsig_predictions = parse_clnsig(&record, &parsed_transcripts);

        if !clnsig_predictions.is_empty() {
            variant.insert(
                "clnsig",
                Bson::Array(
                    clnsig_predictions
                        .into_iter()
                        .map(build_clnsig)
                        .map(Bson::Document)
                        .collect(),
                ),
            );
        }

        let clnsig_onc_predictions = parse_clnsig_onc(&record);

        if !clnsig_onc_predictions.is_empty() {
            variant.insert(
                "clnsig_onc",
                Bson::Array(
                    clnsig_onc_predictions
                        .into_iter()
                        .map(Bson::Document)
                        .collect(),
                ),
            );
        }

        frequencies.extend(parse_frequencies(&record, &parsed_transcripts));

        if !frequencies.is_empty() {
            add_frequencies(&mut variant, &frequencies);
        }

        add_loqus_archive_frequencies(&record, &mut variant, local_archive_info.as_ref());

        set_severity_predictions(&mut variant, &record, &parsed_transcripts);

        parse_conservations(&record, &parsed_transcripts, &mut variant);

        variant.extend(callers);

        if let Some(rank_score_results) = parse_rank_result(&record, &rank_results_header) {
            variant.insert("rank_score_results", Bson::Array(rank_score_results));
        }

        link_gene_panels(&mut variant, gene_to_panels);

        println!("{:#?}\n", variant);
        variant_count += 1;
    }

    println!("Parsed {} variants from {}", variant_count, path);
}
