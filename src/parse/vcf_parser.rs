use rust_htslib::bcf::{Read, Reader};

pub fn process_vcf(path: &str) {

    let mut vcf = Reader::from_path(path)
        .expect("couldn't open input vcf");

    for result in vcf.records() {
        let record = result.unwrap();
         println!("{:?}", record);
    }

}