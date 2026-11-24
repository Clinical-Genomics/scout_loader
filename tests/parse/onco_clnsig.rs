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

#[test]
fn test_parse_clnsig_onc_review_status_with_commas() {
    // GIVEN an oncogenicity review status with literal commas
    let record = crate::vcf_record::record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">",
            "##INFO=<ID=ONC,Number=.,Type=String,Description=\"Oncogenicity classification\">",
            "##INFO=<ID=ONCREVSTAT,Number=.,Type=String,Description=\"Oncogenicity review status\">",
        ],
        "1\t100\t.\tA\tG\t.\tPASS\tCLNVID=376460;ONC=Likely_oncogenic;ONCREVSTAT=criteria_provided,_single_submitter",
    );

    // WHEN parsing the oncogenicity classification
    let onc = parse_clnsig_onc(&record);

    // THEN the full review status is kept
    assert_eq!(onc.len(), 1);
    assert_eq!(
        onc[0].get_str("revstat"),
        Ok("criteria_provided,single_submitter")
    );
}

#[test]
fn test_parse_clnsig_onc_repeats_full_disease_list_for_each_classification() {
    // GIVEN multiple oncogenicity classifications and diseases
    let record = crate::vcf_record::record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">",
            "##INFO=<ID=ONC,Number=.,Type=String,Description=\"Oncogenicity classification\">",
            "##INFO=<ID=ONCDN,Number=1,Type=String,Description=\"Oncogenicity description\">",
        ],
        "1\t100\t.\tA\tG\t.\tPASS\tCLNVID=376460;ONC=Oncogenic/Likely_oncogenic;ONCDN=Neoplasm|Lymphoma/leukemia",
    );

    // WHEN parsing the oncogenicity classifications
    let onc = parse_clnsig_onc(&record);

    // THEN each classification has the complete disease list
    assert_eq!(
        onc,
        vec![
            doc! {
                "accession": 376460_i64,
                "value": "oncogenic",
                "revstat": "",
                "dn": "Neoplasm,Lymphoma/leukemia",
            },
            doc! {
                "accession": 376460_i64,
                "value": "likely_oncogenic",
                "revstat": "",
                "dn": "Neoplasm,Lymphoma/leukemia",
            },
        ]
    );
}
