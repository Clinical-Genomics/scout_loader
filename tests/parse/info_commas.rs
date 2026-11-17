//! Regression tests for String INFO values containing commas.
//!
//! rust-htslib splits String INFO values on commas, and `parse_info_string`
//! used to return only the first part. cyvcf2 (used by scout) returns the
//! whole string, so e.g. a ClinVar review status of
//! `criteria_provided,_multiple_submitters,_no_conflicts` was stored as
//! `criteria_provided` by scout_loader.

use mongodb::bson::{Bson, Document, doc};
use rust_htslib::bcf::{Read, Reader, Record};
use scout_loader::parse::info::{parse_custom_data, parse_info_string};
use scout_loader::parse::onco_clnsig::parse_clnsig_onc;
use scout_loader::parse::rank_scores::parse_rank_score_other;
use scout_loader::parse::strs::set_str_info;
use scout_loader::parse::vep::clnsig::parse_clnsig;
use std::io::Write;
use tempfile::NamedTempFile;

const META: &str = "\
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##INFO=<ID=NUMBER_ONE,Number=1,Type=String,Description=\"Number=1 string\">
##INFO=<ID=NUMBER_DOT,Number=.,Type=String,Description=\"Number=. string\">
##INFO=<ID=CLNVID,Number=1,Type=String,Description=\"ClinVar Variation ID\">
##INFO=<ID=CLNACC,Number=.,Type=String,Description=\"Variant Accession and Versions\">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description=\"Review status\">
##INFO=<ID=ONC,Number=.,Type=String,Description=\"Oncogenicity classification\">
##INFO=<ID=ONCREVSTAT,Number=.,Type=String,Description=\"Oncogenicity review status\">
##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description=\"Custom annotations\">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description=\"Motifs that the tandem repeat is composed of\">
##INFO=<ID=PathologicStruc,Number=1,Type=String,Description=\"Struc components counting towards pathogenicity\">
##INFO=<ID=SourceDisplay,Number=1,Type=String,Description=\"Source for variant definition, display\">
##INFO=<ID=MivmirScore,Number=1,Type=Float,Description=\"Mivmir score\">
##INFO=<ID=MivmirExplanation,Number=1,Type=String,Description=\"Mivmir explanation\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

fn record(info: &str) -> Record {
    let mut file = NamedTempFile::new().unwrap();
    writeln!(file, "{META}1\t100\t.\tA\tG\t.\tPASS\t{info}").unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();
    reader.records().next().unwrap().unwrap()
}

#[test]
fn test_parse_info_string_keeps_commas() {
    // GIVEN String INFO fields with commas, declared with Number=1 and Number=.
    let record = record("NUMBER_ONE=a,b,c;NUMBER_DOT=x,,y");

    // THEN the full value is returned, including empty parts
    assert_eq!(
        parse_info_string(&record, b"NUMBER_ONE").as_deref(),
        Some("a,b,c")
    );
    assert_eq!(
        parse_info_string(&record, b"NUMBER_DOT").as_deref(),
        Some("x,,y")
    );

    // AND missing fields are still None
    assert_eq!(parse_info_string(&record, b"CLNSIG"), None);
}

#[test]
fn test_parse_clnsig_review_status_with_commas() {
    // GIVEN a ClinVar annotation as written in the ClinVar VCF, with literal commas
    let record = record(
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
    let record = record("CLNACC=RCV000195198.2,RCV000195199.1;CLNSIG=0,5;CLNREVSTAT=single,mult");

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

#[test]
fn test_parse_clnsig_onc_review_status_with_commas() {
    // GIVEN an oncogenicity review status with literal commas
    let record =
        record("CLNVID=376460;ONC=Likely_oncogenic;ONCREVSTAT=criteria_provided,_single_submitter");

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
fn test_parse_scout_custom_with_several_pairs() {
    // GIVEN a SCOUT_CUSTOM annotation with two key/value pairs
    // (from scout's customannotations_one.vcf.gz)
    let record = record("SCOUT_CUSTOM=key1|val1,key2|val2");

    // WHEN parsing it as the loader does
    let custom = parse_custom_data(parse_info_string(&record, b"SCOUT_CUSTOM"));

    // THEN both pairs are kept
    assert_eq!(
        custom,
        Some(Bson::Array(vec![
            Bson::Array(vec!["key1".into(), "val1".into()]),
            Bson::Array(vec!["key2".into(), "val2".into()]),
        ]))
    );
}

#[test]
fn test_set_str_info_with_commas() {
    // GIVEN TRGT/stranger STR annotations containing commas
    // (from scout's NIST.trgt.stranger.vcf.gz)
    let record = record(
        "MOTIFS=AAAAT,GAAAT;PathologicStruc=[0, 1, 2, 6];SourceDisplay=Paucar, Nilsson, Engvall et al medRxiv 2023-10-03",
    );
    let mut variant = Document::new();

    // WHEN setting the STR info
    set_str_info(&record, &mut variant);

    // THEN the full values are stored, as scout stores them
    assert_eq!(variant.get_str("str_motifs"), Ok("AAAAT,GAAAT"));
    assert_eq!(variant.get_str("str_pathologic_struc"), Ok("[0, 1, 2, 6]"));
    assert_eq!(
        variant
            .get_document("str_source")
            .unwrap()
            .get_str("display"),
        Ok("Paucar, Nilsson, Engvall et al medRxiv 2023-10-03")
    );
}

#[test]
fn test_parse_mivmir_explanation_with_several_entries() {
    // GIVEN a Mivmir explanation with several entries
    // (format from scout's test_parse_rank_score_mivmir_gicam)
    let record = record(
        "MivmirScore=0.0045;MivmirExplanation=[CSQ_CLINVAR_CLNSIG=0.48,CSQ_CLINVAR_CLNREVSTAT=0.24,]",
    );

    // WHEN parsing other rank scores
    let parsed = parse_rank_score_other(&record).unwrap();

    // THEN all entries are part of the description
    assert_eq!(
        parsed.get_document("Mivmir").unwrap().get_document("desc"),
        Ok(&doc! { "CSQ_CLINVAR_CLNSIG": 0.48, "CSQ_CLINVAR_CLNREVSTAT": 0.24 })
    );
}
