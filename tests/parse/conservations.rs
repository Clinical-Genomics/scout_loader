//! Conservation parsing tests, ported from scout's
//! `tests/parse/test_parse_conservation.py`.

use crate::vcf_record::record;
use mongodb::bson::{Bson, Document, doc};
use scout_loader::parse::conservations::parse_conservations;

const META: &[&str] = &[
    "##contig=<ID=1,length=249250621>",
    "##INFO=<ID=dbNSFP_GERP___RS,Number=.,Type=Float,Description=\"GERP++ RS\">",
    "##INFO=<ID=dbNSFP_phastCons100way_vertebrate,Number=.,Type=Float,Description=\"phastCons\">",
    "##INFO=<ID=dbNSFP_phyloP100way_vertebrate,Number=.,Type=Float,Description=\"phyloP\">",
];

fn conservations(info: &str, transcripts: &[Document]) -> Document {
    let record = record(META, &format!("1\t880086\t.\tT\tC\t.\tPASS\t{info}"));
    let mut variant = Document::new();

    parse_conservations(&record, transcripts, &mut variant);

    variant
}

fn strings(values: &[&str]) -> Bson {
    Bson::Array(values.iter().map(|value| Bson::from(*value)).collect())
}

#[test]
fn test_parse_conservation() {
    // GIVEN a variant with GERP information in INFO
    let variant = conservations("dbNSFP_GERP___RS=3.7", &[]);

    // THEN it is conserved
    assert_eq!(
        variant.get("gerp_conservation"),
        Some(&strings(&["Conserved (3.7)"]))
    );
}

#[test]
fn test_parse_conservation_multiple_terms() {
    // GIVEN a variant with multiple GERP annotations
    let variant = conservations("dbNSFP_GERP___RS=3.7,-0.34", &[]);

    // THEN all terms are returned
    assert_eq!(
        variant.get("gerp_conservation"),
        Some(&strings(&["Conserved (3.7)", "NotConserved (-0.34)"]))
    );
}

#[test]
fn test_parse_conservations() {
    // GIVEN a variant with all conservation annotations in INFO
    let variant = conservations(
        "dbNSFP_GERP___RS=4.6,0;dbNSFP_phastCons100way_vertebrate=0.8;dbNSFP_phyloP100way_vertebrate=2.4",
        &[],
    );

    // THEN all terms are returned, compared against each threshold
    // (gerp >= 2.0, phast >= 0.8, phylop >= 2.5)
    assert_eq!(
        variant,
        doc! {
            "gerp_conservation": ["Conserved (4.6)", "NotConserved (0.0)"],
            "phast_conservation": ["Conserved (0.8)"],
            "phylop_conservation": ["NotConserved (2.4)"],
        }
    );
}

#[test]
fn test_parse_conservation_csq() {
    // GIVEN a transcript with multiple conservation annotations, and nothing in INFO
    let transcript = doc! { "gerp": "0&4.6", "phast": "0.8", "phylop": "2.4" };

    let variant = conservations(".", &[transcript]);

    // THEN the transcript values are used
    assert_eq!(
        variant,
        doc! {
            "gerp_conservation": ["NotConserved (0.0)", "Conserved (4.6)"],
            "phast_conservation": ["Conserved (0.8)"],
            "phylop_conservation": ["NotConserved (2.4)"],
        }
    );
}

#[test]
fn test_parse_conservation_info_takes_precedence() {
    // GIVEN GERP in both INFO and the transcript, and phast only in the transcript
    let transcript = doc! { "gerp": "0", "phast": "0.9&not_a_number" };

    let variant = conservations("dbNSFP_GERP___RS=5", &[transcript]);

    // THEN INFO wins for GERP, the transcript is used for phast (ignoring
    // invalid values), and phylop is empty
    assert_eq!(
        variant,
        doc! {
            "gerp_conservation": ["Conserved (5.0)"],
            "phast_conservation": ["Conserved (0.9)"],
            "phylop_conservation": [],
        }
    );
}

#[test]
fn test_parse_conservation_only_first_transcript_is_used() {
    // GIVEN two transcripts where only the second has conservation info
    let transcripts = [doc! {}, doc! { "gerp": "4.6" }];

    let variant = conservations(".", &transcripts);

    // THEN nothing is parsed
    assert_eq!(variant.get("gerp_conservation"), Some(&Bson::Array(vec![])));
}


#[test]
fn test_parse_conservation_missing() {
    // GIVEN no conservation information at all
    let variant = conservations(".", &[]);

    // THEN empty lists are set for all keys
    assert_eq!(
        variant,
        doc! {
            "gerp_conservation": [],
            "phast_conservation": [],
            "phylop_conservation": [],
        }
    );
}
