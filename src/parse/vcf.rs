use rust_htslib::bcf::{Read, Reader};
use crate::parse::coordinates::parse_coordinates;
use crate::parse::alleles::parse_alleles;
use crate::models::variant::Variant;
use crate::models::variant::VariantCategory;


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
///
/// # Panics
///
/// Panics if the VCF file cannot be opened or if a record cannot be read.
pub fn process_vcf(path: &str, category: VariantCategory) {

    let mut vcf = Reader::from_path(path)
        .expect("couldn't open input vcf");

    let header = vcf.header().clone();

    for result in vcf.records() {
        let record = result.unwrap();
        let (chromosome, position, end) = parse_coordinates(&record, &header);
        let (reference, alternative) = parse_alleles(&record, category);

        let variant = Variant {
            chromosome: chromosome,
            position: position,
            end: end,
            reference: reference,
            alternative: alternative
        };
        println!("{:#?}", variant);
            
    }

}