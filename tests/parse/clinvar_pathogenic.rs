use rust_htslib::bcf::{self, Read};
use scout_loader::parse::clinvar_pathogenic::is_pathogenic_transcripts;
use std::collections::HashMap;
use std::io::Write;

#[test]
fn test_is_pathogenic_transcripts() {
    let mut vcf = tempfile::NamedTempFile::new().unwrap();

    writeln!(
        vcf,
        r#"##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|CLNSIG|CLINVAR_CLNSIG">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	17428382	.	C	T	54	PASS	CSQ=T|Benign|Conflicting_classifications_of_pathogenicity"#
    )
    .unwrap();

    let mut reader = bcf::Reader::from_path(vcf.path()).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    let field_indices =
        HashMap::from([("CLNSIG".to_string(), 1), ("CLINVAR_CLNSIG".to_string(), 2)]);

    assert!(is_pathogenic_transcripts(&record, &field_indices));
}
