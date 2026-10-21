use std::io::Write;

use rust_htslib::bcf::{Read, Reader};
use tempfile::NamedTempFile;

use scout_loader::parse::genotypes::parse_format_entry;

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
