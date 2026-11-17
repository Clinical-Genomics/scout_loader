use crate::vcf_record::record_from_vcf as parse_test_record;
use mongodb::bson::doc;
use scout_loader::parse::vep::clnsig::parse_clnsig;

#[test]
fn test_parse_clnsig_from_info() {
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=11,length=135086622>
##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description=\"Clinical significance\">
##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description=\"Review status\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
11\t108106541\t.\tTATCTC\tT\t205\tPASS\tCLNVID=185501;CLNREVSTAT=criteria_provided%2C_multiple_submitters%2C_no_conflicts;CLNSIG=Pathogenic/Likely_pathogenic
";

    let record = parse_test_record(vcf);
    let transcripts = vec![doc! {
        "clnsig": ["likely_pathogenic", "pathogenic"],
    }];

    let result = parse_clnsig(&record, &transcripts);

    assert_eq!(
        result,
        vec![
            doc! {
                "value": "pathogenic",
                "accession": 185501,
                "revstat": "criteria_provided,multiple_submitters,no_conflicts",
            },
            doc! {
                "value": "likely_pathogenic",
                "accession": 185501,
                "revstat": "criteria_provided,multiple_submitters,no_conflicts",
            },
        ]
    );
}

#[test]
fn test_parse_clnsig_from_csq_when_info_is_missing() {
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=11,length=135086622>
##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">
##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description=\"Review status\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
11\t108106541\t.\tTATCTC\tT\t205\tPASS\tCLNVID=185501;CLNREVSTAT=criteria_provided%2C_multiple_submitters%2C_no_conflicts
";

    let record = parse_test_record(vcf);
    let transcripts = vec![doc! {
        "clnsig": ["likely_pathogenic", "pathogenic"],
    }];

    let result = parse_clnsig(&record, &transcripts);

    assert_eq!(
        result,
        vec![
            doc! {
                "value": "likely_pathogenic",
                "accession": 185501,
                "revstat": "criteria_provided,multiple_submitters,no_conflicts",
            },
            doc! {
                "value": "pathogenic",
                "accession": 185501,
                "revstat": "criteria_provided,multiple_submitters,no_conflicts",
            },
        ]
    );
}
