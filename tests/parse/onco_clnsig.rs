use mongodb::bson::doc;
use rust_htslib::bcf::{Read, Reader};
use scout_loader::parse::onco_clnsig::parse_clnsig_onc;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_parse_clnsig_onc() {
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=2,length=243199373>
##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">
##INFO=<ID=ONC,Number=1,Type=String,Description=\"Oncogenicity classification\">
##INFO=<ID=ONCDN,Number=1,Type=String,Description=\"Oncogenicity description\">
##INFO=<ID=ONCREVSTAT,Number=1,Type=String,Description=\"Review status\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
2\t16082317\t.\tC\tT\t311\tPASS\tCLNVID=376460;ONC=Likely_oncogenic;ONCDN=Neoplasm;ONCREVSTAT=criteria_provided%2C_single_submitter
";

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(vcf.as_bytes()).unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    let result = parse_clnsig_onc(&record);

    assert_eq!(
        result,
        vec![doc! {
            "accession": 376460_i64,
            "value": "likely_oncogenic",
            "revstat": "criteria_provided,single_submitter",
            "dn": "Neoplasm",
        }]
    );
}
