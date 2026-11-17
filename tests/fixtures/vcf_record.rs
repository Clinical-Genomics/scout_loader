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

/// Build a VEP-annotated record from a CSQ format string and a CSQ value.
///
/// The CSQ header is written the way VEP does it (`... Format: A|B|C`) and the
/// field indices are derived with the loader's own `parse_vep_header`, so the
/// returned indices match what `parse_vcf` would use for the same file.
///
/// `extra_info` is appended to the INFO column (e.g. `"COSMIC=COSV1"`); any
/// keys used there must be declared in `extra_meta`.
pub fn vep_record(
    csq_format: &str,
    csq: &str,
    extra_meta: &[&str],
    extra_info: &str,
) -> (Record, std::collections::HashMap<String, usize>) {
    let csq_header = format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {csq_format}\">"
    );

    let mut meta = vec!["##contig=<ID=1,length=249250621>", csq_header.as_str()];
    meta.extend_from_slice(extra_meta);

    let info = if extra_info.is_empty() {
        format!("CSQ={csq}")
    } else {
        format!("CSQ={csq};{extra_info}")
    };

    let record = record(&meta, &format!("1\t38592898\t.\tT\tC\t.\tPASS\t{info}"));

    let field_indices = scout_loader::parse::header::parse_vep_header(record.header())
        .into_iter()
        .enumerate()
        .map(|(index, field)| (field, index))
        .collect();

    (record, field_indices)
}
