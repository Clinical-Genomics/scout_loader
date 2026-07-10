use rust_htslib::bcf::{Read, Reader};
use crate::parse::coordinates::parse_coordinates;
use crate::models::variant::Variant;

pub fn process_vcf(path: &str, category: &str) {

    let mut vcf = Reader::from_path(path)
        .expect("couldn't open input vcf");

    let header = vcf.header().clone();

    for result in vcf.records() {
        let record = result.unwrap();
        let (chromosome, position, end) = parse_coordinates(&record, &header);

        let variant = Variant {
            chromosome: chromosome,
            position: position,
            end: end,
        };
        println!("{:#?}", variant);
            
    }

}