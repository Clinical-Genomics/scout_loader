use rust_htslib::bcf::{Read, Reader};
use crate::parse::coordinates::parse_coordinates;
use crate::parse::alleles::parse_alleles;
use crate::parse::quality::parse_filters;
use crate::models::variant::Variant;
use crate::models::variant::VariantCategory;
use crate::models::variant::VariantType;


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
///
/// # Panics
///
/// Panics if the VCF file cannot be opened or if a record cannot be read.
pub fn process_vcf(path: &str, category: VariantCategory, variant_type: VariantType, case_id: &str) {

    let mut vcf = Reader::from_path(path)
        .expect("couldn't open input vcf");

    let header = vcf.header().clone();

    for result in vcf.records() {
        let record = result.unwrap();
        let (chromosome, position, end) = parse_coordinates(&record, &header);
        let (reference, alternative) = parse_alleles(&record, category);
        let filters = parse_filters(&record, &header);

        let variant = Variant {
            case_id: case_id.to_string(),
            r#type: variant_type.to_string(),
            chromosome: chromosome,
            position: position,
            end: end,
            reference: reference,
            alternative: alternative,
            filters: filters
        };
        println!("{:#?}", variant);
            
    }

}