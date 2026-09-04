use rust_htslib::bcf::{Read, Reader};
use scout_loader::parse::genetic_models::parse_genetic_models;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_parse_genetic_models() {
    let vcf = b"##fileformat=VCFv4.2
##INFO=<ID=GeneticModels,Number=1,Type=String,Description=\"Genetic models\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tT\t.\t.\tGeneticModels=internal_id:XD_dn|XR_dn
";

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(vcf).unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();
    let mut record = reader.empty_record();

    reader
        .read(&mut record)
        .expect("Failed to read VCF")
        .expect("Failed to parse VCF record");

    let genetic_models = parse_genetic_models(&record);

    assert_eq!(
        genetic_models,
        vec!["XD_dn".to_string(), "XR_dn".to_string()]
    );
}
