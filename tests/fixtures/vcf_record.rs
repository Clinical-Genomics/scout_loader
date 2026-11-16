//! Helpers for building `rust_htslib` VCF records from inline VCF text.
//!
//! These play the same role as the `cyvcf2_variant` / `mock_variant` fixtures
//! in the scout test suite: a test describes a variant as a VCF line and gets
//! back a real `Record` (with its header attached via `record.header()`).

#![allow(dead_code)]

use rust_htslib::bcf::{Read, Reader, Record};
use std::io::Write;
use tempfile::NamedTempFile;

/// Parse the first record of a complete VCF text.
pub fn record_from_vcf(vcf: &str) -> Record {
    records_from_vcf(vcf)
        .into_iter()
        .next()
        .expect("VCF contains no records")
}

/// Parse all records of a complete VCF text.
pub fn records_from_vcf(vcf: &str) -> Vec<Record> {
    let mut file = NamedTempFile::new().expect("failed to create temp VCF");
    file.write_all(vcf.as_bytes())
        .expect("failed to write temp VCF");

    let mut reader = Reader::from_path(file.path()).expect("failed to open temp VCF");
    reader
        .records()
        .map(|record| record.expect("failed to parse VCF record"))
        .collect()
}

/// Build a record from `##` meta lines and a single tab-separated data line.
///
/// `##fileformat` and the `#CHROM` column header are added automatically. If
/// the data line has more than eight columns, a `FORMAT` column and one sample
/// column per extra field (named `SAMPLE1`, `SAMPLE2`, ...) are declared.
pub fn record(meta_lines: &[&str], data_line: &str) -> Record {
    let columns = data_line.split('\t').count();

    let mut header = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if columns > 8 {
        header.push_str("\tFORMAT");
        for sample in 1..columns - 8 {
            header.push_str(&format!("\tSAMPLE{sample}"));
        }
    }

    let mut vcf = String::from("##fileformat=VCFv4.2\n");
    for line in meta_lines {
        vcf.push_str(line);
        vcf.push('\n');
    }
    vcf.push_str(&header);
    vcf.push('\n');
    vcf.push_str(data_line);
    vcf.push('\n');

    record_from_vcf(&vcf)
}
