use crate::vcf_record::{record, record_from_vcf as parse_test_record};
use mongodb::bson::doc;
use rust_htslib::bcf::Record;
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

fn clnsig_record(info: &str) -> Record {
    record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">",
            "##INFO=<ID=CLNACC,Number=.,Type=String,Description=\"Variant Accession and Versions\">",
            "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">",
            "##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description=\"Review status\">",
        ],
        &format!("1\t100\t.\tA\tG\t.\tPASS\t{info}"),
    )
}

#[test]
fn test_parse_clnsig_review_status_with_commas() {
    // GIVEN a ClinVar annotation as written in the ClinVar VCF, with literal commas
    let record = clnsig_record(
        "CLNVID=12345;CLNSIG=Pathogenic,_low_penetrance;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts",
    );

    // WHEN parsing the ClinVar significance
    let clnsig = parse_clnsig(&record, &[]);

    // THEN the full significance and review status are kept
    assert_eq!(
        clnsig,
        vec![doc! {
            "value": "pathogenic,low_penetrance",
            "accession": 12345,
            "revstat": "criteria_provided,multiple_submitters,no_conflicts",
        }]
    );
}

#[test]
fn test_parse_clnsig_old_format_multiple_accessions() {
    // GIVEN the old numeric ClinVar format with one entry per accession
    // (from scout's demo 643594.research.vcf.gz)
    let record =
        clnsig_record("CLNACC=RCV000195198.2,RCV000195199.1;CLNSIG=0,5;CLNREVSTAT=single,mult");

    // WHEN parsing the ClinVar significance
    let clnsig = parse_clnsig(&record, &[]);

    // THEN all accessions are parsed, not only the first
    assert_eq!(
        clnsig,
        vec![
            doc! { "value": 0, "accession": "RCV000195198.2", "revstat": "single" },
            doc! { "value": 5, "accession": "RCV000195199.1", "revstat": "mult" },
        ]
    );
}
