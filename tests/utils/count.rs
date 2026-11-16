use flate2::Compression;
use flate2::write::GzEncoder;
use scout_loader::utils::count::count_vcf_records;
use std::io::Write;
use tempfile::{Builder, NamedTempFile};

const VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\tPASS\t.
1\t200\t.\tC\tT\t.\tPASS\t.
1\t300\t.\tG\tA\t.\tPASS\t.
";

#[test]
fn test_count_vcf_records_plain() {
    // GIVEN an uncompressed VCF with 2 header lines and 3 records
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(VCF.as_bytes()).unwrap();

    // WHEN counting
    let count = count_vcf_records(file.path().to_str().unwrap()).unwrap();

    // THEN all lines, including the header, are counted
    assert_eq!(count, 5);
}

#[test]
fn test_count_vcf_records_gzipped() {
    // GIVEN the same VCF, gzipped
    let file = Builder::new().suffix(".vcf.gz").tempfile().unwrap();
    let mut encoder = GzEncoder::new(file.reopen().unwrap(), Compression::default());
    encoder.write_all(VCF.as_bytes()).unwrap();
    encoder.finish().unwrap();

    // WHEN counting
    let count = count_vcf_records(file.path().to_str().unwrap()).unwrap();

    // THEN the decompressed lines are counted
    assert_eq!(count, 5);
}

#[test]
fn test_count_vcf_records_missing_file() {
    // GIVEN paths that do not exist
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("missing.vcf");
    let gzipped = dir.path().join("missing.vcf.gz");

    // THEN counting fails for both
    assert!(count_vcf_records(plain.to_str().unwrap()).is_err());
    assert!(count_vcf_records(gzipped.to_str().unwrap()).is_err());
}

#[test]
fn test_count_vcf_records_invalid_gzip() {
    // GIVEN a .gz file that is not gzip-compressed
    let mut file = Builder::new().suffix(".vcf.gz").tempfile().unwrap();
    file.write_all(VCF.as_bytes()).unwrap();

    // THEN counting fails
    assert!(count_vcf_records(file.path().to_str().unwrap()).is_err());
}
