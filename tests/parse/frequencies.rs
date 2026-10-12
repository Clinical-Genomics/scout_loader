use std::io::Write;

use rust_htslib::bcf::{Read, Reader};
use tempfile::NamedTempFile;

use scout_loader::parse::frequencies::{GNOMAD_INFO_MAX_KEYS, parse_frequencies};

#[test]
fn test_parse_gnomad_popmax() {
    for key in GNOMAD_INFO_MAX_KEYS {
        let vcf = format!(
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID={key},Number=1,Type=Float,Description=\"gnomAD population maximum allele frequency\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t100\t.\tA\tT\t.\t.\t{key}=0.05\n"
        );

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(vcf.as_bytes()).unwrap();

        let mut reader = Reader::from_path(file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();

        let frequencies = parse_frequencies(&record, &[]);

        assert!(
            (frequencies.get_f64("gnomad_max").unwrap() - 0.05).abs() < 1e-6,
            "Failed to parse INFO field: {key}",
        );
    }
}
