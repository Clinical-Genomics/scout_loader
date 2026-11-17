use std::collections::HashMap;
use std::io::Write;

use rust_htslib::bcf::{Read, Reader, Record};
use tempfile::NamedTempFile;

use scout_loader::models::sample::SampleInfo;
use scout_loader::models::variant::VariantCategory;
use scout_loader::parse::genotypes::{
    get_copy_number, get_phase_set, parse_format_entry, parse_genotypes,
};

/// Write a VCF string to a temporary file and read all of its records.
fn read_records(vcf: &str) -> Vec<Record> {
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(vcf.as_bytes()).unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();
    reader.records().map(|record| record.unwrap()).collect()
}

#[test]
fn test_parse_format_entry() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=ADFL,Number=1,Type=String,Description=\"Flanking depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\t.\tA\t<STR8>\t.\tPASS\t.\tADFL\t33/33
";

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(vcf.as_bytes()).unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();
    let mut record = reader.empty_record();
    reader.read(&mut record).unwrap().unwrap();

    assert_eq!(
        parse_format_entry(&record, 0, b"ADFL"),
        (Some(33), Some(33))
    );
}

#[test]
fn test_get_copy_number() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\t.\tA\t<CNV>\t.\tPASS\t.\tCN\t3
1\t200\t.\tA\t<CNV>\t.\tPASS\t.\tCN\t.
1\t300\t.\tA\t<CNV>\t.\tPASS\t.\tCN\t-1
";

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(vcf.as_bytes()).unwrap();

    let mut reader = Reader::from_path(file.path()).unwrap();

    let mut record = reader.empty_record();
    reader.read(&mut record).unwrap().unwrap();
    assert_eq!(get_copy_number(&record, 0), Some(3));

    let mut record = reader.empty_record();
    reader.read(&mut record).unwrap().unwrap();
    assert_eq!(get_copy_number(&record, 0), None);

    let mut record = reader.empty_record();
    reader.read(&mut record).unwrap().unwrap();
    assert_eq!(get_copy_number(&record, 0), None);
}

#[test]
fn test_get_phase_set() {
    // GIVEN records with a phase set, a missing phase set and a negative phase set
    let records = read_records(
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
1\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:12345\t1|0:100
1\t200\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0/1:.\t0/1:.
1\t300\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:-1\t0|1:7
",
    );

    // THEN the phase set is parsed per sample
    assert_eq!(get_phase_set(&records[0], 0), Some(12345));
    assert_eq!(get_phase_set(&records[0], 1), Some(100));

    // THEN missing values give no phase set
    assert_eq!(get_phase_set(&records[1], 0), None);
    assert_eq!(get_phase_set(&records[1], 1), None);

    // THEN negative values are rejected
    assert_eq!(get_phase_set(&records[2], 0), None);
    assert_eq!(get_phase_set(&records[2], 1), Some(7));
}

#[test]
fn test_get_phase_set_without_phase_set() {
    // GIVEN a record without a PS FORMAT field
    let records = read_records(
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1
",
    );

    // THEN no phase set is returned
    assert_eq!(get_phase_set(&records[0], 0), None);
}

#[test]
fn test_get_phase_set_string_type() {
    // GIVEN a PS FORMAT field declared as a String
    let records = read_records(
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set identifier\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:54321
1\t200\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:abc
",
    );

    // THEN numeric values are parsed and non-numeric values are ignored
    assert_eq!(get_phase_set(&records[0], 0), Some(54321));
    assert_eq!(get_phase_set(&records[1], 0), None);
}

#[test]
fn test_parse_genotypes_phase_set() {
    // GIVEN a phased record where only one sample has a phase set
    let records = read_records(
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPHASED\tUNPHASED
1\t100\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:95\t0/1:.
",
    );

    let sample_mapping = HashMap::from([
        (
            "PHASED".to_string(),
            SampleInfo {
                display_name: "phased".to_string(),
                vcf_index: 0,
                sample_type: None,
            },
        ),
        (
            "UNPHASED".to_string(),
            SampleInfo {
                display_name: "unphased".to_string(),
                vcf_index: 1,
                sample_type: None,
            },
        ),
    ]);

    // WHEN parsing the genotypes
    let (genotypes, _, _) = parse_genotypes(&records[0], &sample_mapping, VariantCategory::Snv);

    let genotype = |sample_id: &str| {
        genotypes
            .iter()
            .find(|gt| gt.get_str("sample_id").unwrap() == sample_id)
            .unwrap()
    };

    // THEN the phased sample has a phase set alongside its phased genotype call
    let phased = genotype("PHASED");
    assert_eq!(phased.get_str("genotype_call").unwrap(), "0|1");
    assert_eq!(phased.get_i32("phase_set").unwrap(), 95);

    // THEN the sample with a missing PS has no phase_set key
    let unphased = genotype("UNPHASED");
    assert_eq!(unphased.get_str("genotype_call").unwrap(), "0/1");
    assert!(!unphased.contains_key("phase_set"));
}
