use std::collections::HashMap;

use rust_htslib::bcf::Record;

use scout_loader::models::sample::SampleInfo;
use scout_loader::models::variant::VariantCategory;
use scout_loader::parse::genotypes::{
    get_copy_number, get_pathologic_struc, get_phase_set, parse_format_entry,
    parse_format_entry_trgt_mc, parse_genotypes, parse_pathologic_struc,
};

use crate::vcf_record::{record, records_from_vcf};

/// PathologicStruc as written by older stranger versions, e.g. `[0, 1, 2, 6]`.
const OLD_PATHOLOGIC_STRUC_HEADER: &str = "##INFO=<ID=PathologicStruc,Number=1,Type=String,Description=\"Struc components counting towards pathogenicity as list of motif numbers\">";
/// PathologicStruc as an Integer list, e.g. `0,1,2,6` (stranger#124).
const NEW_PATHOLOGIC_STRUC_HEADER: &str = "##INFO=<ID=PathologicStruc,Number=.,Type=Integer,Description=\"Struc components counting towards pathogenicity\">";
/// PathologicStruc as a comma-separated String list, e.g. `1,2`.
const STRING_LIST_PATHOLOGIC_STRUC_HEADER: &str = "##INFO=<ID=PathologicStruc,Number=.,Type=String,Description=\"Struc components counting towards pathogenicity\">";
const COPY_NUMBER_HEADER: &str =
    "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">";

#[test]
fn test_parse_format_entry() {
    // GIVEN a slash-separated String FORMAT entry
    let record = record(
        &["##FORMAT=<ID=ADFL,Number=1,Type=String,Description=\"Flanking depth\">"],
        "1\t100\t.\tA\t<STR8>\t.\tPASS\t.\tADFL\t33/33",
    );

    // WHEN parsing the entry for the sample
    // THEN both alleles are returned
    assert_eq!(
        parse_format_entry(&record, 0, b"ADFL"),
        (Some(33), Some(33))
    );
}

#[test]
fn test_get_copy_number() {
    // GIVEN a positive, a missing and a negative copy number
    let copy_number = |value: &str| {
        let record = record(
            &[COPY_NUMBER_HEADER],
            &format!("1\t100\t.\tA\t<CNV>\t.\tPASS\t.\tCN\t{value}"),
        );
        get_copy_number(&record, 0)
    };

    // THEN only the positive copy number is returned
    assert_eq!(copy_number("3"), Some(3));
    assert_eq!(copy_number("."), None);
    assert_eq!(copy_number("-1"), None);
}

#[test]
fn test_get_phase_set() {
    // GIVEN records with a phase set, a missing phase set and a negative phase set
    let records = records_from_vcf(
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
    let records = records_from_vcf(
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
    let records = records_from_vcf(
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
    let records = records_from_vcf(
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

/// Build a TRGT-like single sample record with the given PathologicStruc header and INFO.
fn trgt_record(pathologic_struc_header: &str, info: &str, gt: &str, mc: &str) -> Record {
    record(
        &[
            "##contig=<ID=1,length=249250621>",
            pathologic_struc_header,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=MC,Number=.,Type=String,Description=\"Motif counts\">",
        ],
        &format!("1\t100\t.\tA\t<STR>\t.\tPASS\t{info}\tGT:MC\t{gt}:{mc}"),
    )
}

#[test]
fn test_parse_pathologic_struc() {
    assert_eq!(parse_pathologic_struc("[0, 1, 2, 6]"), vec![0, 1, 2, 6]);
    assert_eq!(parse_pathologic_struc("[0]"), vec![0]);
    assert_eq!(parse_pathologic_struc("0,1,2,6"), vec![0, 1, 2, 6]);
    assert_eq!(parse_pathologic_struc("[]"), Vec::<usize>::new());
}

#[test]
fn test_get_pathologic_struc_old_string_format() {
    // GIVEN the Python list format written by older stranger versions
    let record = trgt_record(
        OLD_PATHOLOGIC_STRUC_HEADER,
        "PathologicStruc=[0, 1, 2, 6]",
        "0/1",
        "1_1_1_1_1,1_1_1_1_1",
    );

    // THEN all indexes are parsed, not only the part before the first comma
    assert_eq!(get_pathologic_struc(&record), Some(vec![0, 1, 2, 6]));

    let record = trgt_record(
        OLD_PATHOLOGIC_STRUC_HEADER,
        "PathologicStruc=[0]",
        "0/1",
        "1_1,1_1",
    );
    assert_eq!(get_pathologic_struc(&record), Some(vec![0]));
}

#[test]
fn test_get_pathologic_struc_new_integer_format() {
    // GIVEN the Number=.,Type=Integer format from stranger#124
    let record = trgt_record(
        NEW_PATHOLOGIC_STRUC_HEADER,
        "PathologicStruc=0,1,2,6",
        "0/1",
        "1_1_1_1_1,1_1_1_1_1",
    );
    assert_eq!(get_pathologic_struc(&record), Some(vec![0, 1, 2, 6]));

    // AND a single index, including index 0
    let record = trgt_record(
        NEW_PATHOLOGIC_STRUC_HEADER,
        "PathologicStruc=0",
        "0/1",
        "1_1,1_1",
    );
    assert_eq!(get_pathologic_struc(&record), Some(vec![0]));
}

#[test]
fn test_get_pathologic_struc_string_list_format() {
    // GIVEN a comma-separated list declared as String
    let record = trgt_record(
        STRING_LIST_PATHOLOGIC_STRUC_HEADER,
        "PathologicStruc=1,2",
        "0/1",
        "1_1_1,1_1_1",
    );
    assert_eq!(get_pathologic_struc(&record), Some(vec![1, 2]));
}

#[test]
fn test_get_pathologic_struc_missing() {
    let record = trgt_record(NEW_PATHOLOGIC_STRUC_HEADER, ".", "0/1", "1_1,1_1");
    assert_eq!(get_pathologic_struc(&record), None);
}

#[test]
fn test_parse_format_entry_trgt_mc_pathologic_struc() {
    // GIVEN motif counts for five motifs, where motifs 0, 1, 2 and 6 (if present) count
    // towards the pathologic size, in both the old and the new format
    for (header, info) in [
        (OLD_PATHOLOGIC_STRUC_HEADER, "PathologicStruc=[0, 1, 2, 6]"),
        (NEW_PATHOLOGIC_STRUC_HEADER, "PathologicStruc=0,1,2,6"),
    ] {
        let record = trgt_record(header, info, "0/1", "1_2_3_100_200,10_20_30_100_200");

        // THEN only the pathologic motifs are summed
        assert_eq!(
            parse_format_entry_trgt_mc(&record, 0),
            (Some(1 + 2 + 3), Some(10 + 20 + 30)),
            "{info}"
        );
    }
}

#[test]
fn test_parse_format_entry_trgt_mc_single_index_zero() {
    // GIVEN only motif 0 counts, in the new format with a single value 0
    let record = trgt_record(
        NEW_PATHOLOGIC_STRUC_HEADER,
        "PathologicStruc=0",
        "0/1",
        "5_100,7_100",
    );

    // THEN the other motifs are not counted
    assert_eq!(parse_format_entry_trgt_mc(&record, 0), (Some(5), Some(7)));
}

#[test]
fn test_parse_format_entry_trgt_mc_without_pathologic_struc() {
    // GIVEN no PathologicStruc
    let record = trgt_record(NEW_PATHOLOGIC_STRUC_HEADER, ".", "0/1", "5_100,7_100");

    // THEN all motifs are counted
    assert_eq!(
        parse_format_entry_trgt_mc(&record, 0),
        (Some(105), Some(107))
    );
}
