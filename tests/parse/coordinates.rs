//! Coordinate parsing tests, ported from scout's
//! `tests/parse/test_parse_coordinates.py` and extended with scout_loader
//! specifics (chr-prefix normalisation, MEIs, cytobands, MATEID).

use crate::vcf_record::record;
use rust_htslib::bcf::Record;
use scout_loader::models::cytoband::Cytoband;
use scout_loader::models::variant::{Coordinates, VariantCategory};
use scout_loader::parse::coordinates::{
    get_cytoband_coordinates, get_end_chrom, normalize_chromosome, parse_coordinates, sv_end,
    sv_length,
};
use std::collections::HashMap;

const SV_META: &[&str] = &[
    "##contig=<ID=1,length=249250621>",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV\">",
    "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"Mate breakend ID\">",
];

fn coordinates(record: &Record, category: VariantCategory) -> Coordinates {
    parse_coordinates(record, record.header(), &HashMap::new(), &category)
}

fn sv_record(ref_allele: &str, alt: &str, info: &str) -> Record {
    record(
        SV_META,
        &format!("1\t724779\t.\t{ref_allele}\t{alt}\t.\tPASS\t{info}"),
    )
}

fn test_cytobands() -> HashMap<String, Vec<Cytoband>> {
    HashMap::from([(
        "1".to_string(),
        vec![
            Cytoband {
                start: 0,
                end: 2_300_000,
                name: "p36.33".to_string(),
            },
            Cytoband {
                start: 2_300_000,
                end: 5_400_000,
                name: "p36.32".to_string(),
            },
        ],
    )])
}

#[test]
fn test_parse_coordinates_snv() {
    // GIVEN a simple SNV
    let record = record(
        &["##contig=<ID=1,length=249250621>"],
        "1\t880086\t.\tT\tC\t.\tPASS\t.",
    );

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Snv);

    // THEN they are correctly parsed
    assert_eq!(coordinates.chromosome, "1");
    assert_eq!(coordinates.position, 880086);
    assert_eq!(coordinates.end, 880086);
    assert_eq!(coordinates.end_chrom, "1");
    assert_eq!(coordinates.length, 1);
    assert_eq!(coordinates.sub_category, "snv");
    assert_eq!(coordinates.mate_id, None);
    assert_eq!(coordinates.cytoband_start, None);
    assert_eq!(coordinates.cytoband_end, None);
}

#[test]
fn test_parse_coordinates_indel() {
    // GIVEN an insertion
    let record = record(
        &["##contig=<ID=1,length=249250621>"],
        "1\t880086\t.\tA\tACCC\t.\tPASS\t.",
    );

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Snv);

    // THEN it is an indel with length |len(ALT) - len(REF)|
    assert_eq!(coordinates.position, 880086);
    assert_eq!(coordinates.end, 880086);
    assert_eq!(coordinates.sub_category, "indel");
    assert_eq!(coordinates.length, 3);
}

#[test]
fn test_parse_coordinates_deletion_end_from_ref() {
    // GIVEN a deletion spanning several reference bases
    let record = record(
        &["##contig=<ID=1,length=249250621>"],
        "1\t880086\t.\tACGT\tA\t.\tPASS\t.",
    );

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Snv);

    // THEN the end covers the reference allele
    assert_eq!(coordinates.end, 880089);
    assert_eq!(coordinates.sub_category, "indel");
    assert_eq!(coordinates.length, 3);
}

#[test]
fn test_parse_coordinates_strips_chr_prefix() {
    // GIVEN a variant on a chr-prefixed contig
    let record = record(
        &["##contig=<ID=chr1,length=248956422>"],
        "chr1\t1000\t.\tA\tG\t.\tPASS\t.",
    );

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Snv);

    // THEN the prefix is removed from both chromosome and end chromosome
    assert_eq!(coordinates.chromosome, "1");
    assert_eq!(coordinates.end_chrom, "1");
}

#[test]
fn test_parse_coordinates_translocation() {
    // GIVEN a BND to another contig
    let record = sv_record("N", "N[hs37d5:12060532[", "SVTYPE=BND");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN end and end_chrom come from the ALT field and the length is "infinite"
    assert_eq!(coordinates.position, 724779);
    assert_eq!(coordinates.end, 12060532);
    assert_eq!(coordinates.end_chrom, "hs37d5");
    assert_eq!(coordinates.length, 1_000_000_000);
    assert_eq!(coordinates.sub_category, "bnd");
}

#[test]
fn test_parse_coordinates_translocation_same_chrom() {
    // GIVEN a BND to the same chromosome
    let record = sv_record("N", "N[1:12060532[", "SVTYPE=BND");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN the length is the distance between the breakends
    assert_eq!(coordinates.position, 724779);
    assert_eq!(coordinates.end, 12060532);
    assert_eq!(coordinates.end_chrom, "1");
    assert_eq!(coordinates.length, 12060532 - 724779);
    assert_eq!(coordinates.sub_category, "bnd");
}

#[test]
fn test_parse_coordinates_translocation_2() {
    // GIVEN a BND to a contig name containing a dot
    let record = sv_record("N", "N[GL000232.1:25141[", "SVTYPE=BND");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN the full contig name is kept
    assert_eq!(coordinates.end, 25141);
    assert_eq!(coordinates.end_chrom, "GL000232.1");
    assert_eq!(coordinates.length, 1_000_000_000);
    assert_eq!(coordinates.sub_category, "bnd");
}

#[test]
fn test_parse_coordinates_translocation_chr_prefixed_mate() {
    // GIVEN a BND whose mate contig is chr-prefixed
    let record = sv_record("N", "]chr2:321682]N", "SVTYPE=BND;MATEID=bnd_2");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN the mate chromosome is normalised and the mate id is kept
    assert_eq!(coordinates.end_chrom, "2");
    assert_eq!(coordinates.end, 321682);
    assert_eq!(coordinates.mate_id.as_deref(), Some("bnd_2"));
}

#[test]
fn test_parse_coordinates_sgl_bnd() {
    // GIVEN a single-breakend BND with the SGL svtype alias
    let record = sv_record("G", ".G", "SVTYPE=SGL");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN it is a bnd with END == POS and unknown length
    assert_eq!(coordinates.position, 724779);
    assert_eq!(coordinates.end, 724779);
    assert_eq!(coordinates.end_chrom, "1");
    assert_eq!(coordinates.length, -1);
    assert_eq!(coordinates.sub_category, "bnd");
}

#[test]
fn test_parse_coordinates_sv_type_from_alt() {
    // GIVEN a symbolic deletion without SVTYPE
    let record = sv_record("N", "<DEL>", "END=725711");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN the sub category is inferred from ALT and END is used
    assert_eq!(coordinates.sub_category, "del");
    assert_eq!(coordinates.end, 725711);
    assert_eq!(coordinates.length, 725711 - 724779);
}

#[test]
fn test_parse_coordinates_sgl_bnd_without_svtype() {
    // GIVEN a single-breakend without SVTYPE
    let record = sv_record("G", "GTCA.", ".");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN it is inferred to be a bnd
    assert_eq!(coordinates.sub_category, "bnd");
    assert_eq!(coordinates.end, 724779);
}

#[test]
fn test_parse_coordinates_cnvnator_deletion() {
    // GIVEN a CNVnator-style deletion with a negative SVLEN
    let record = sv_record("N", "<DEL>", "SVTYPE=DEL;END=734778;SVLEN=-10000");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN the length is the absolute SVLEN
    assert_eq!(coordinates.sub_category, "del");
    assert_eq!(coordinates.end, 734778);
    assert_eq!(coordinates.length, 10000);
}

#[test]
fn test_parse_coordinates_dup_end_from_svlen() {
    // GIVEN a duplication with SVLEN but no END
    let record = sv_record("N", "<DUP>", "SVTYPE=DUP;SVLEN=500");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN END is derived from POS + SVLEN
    assert_eq!(coordinates.sub_category, "dup");
    assert_eq!(coordinates.end, 724779 + 500);
    assert_eq!(coordinates.length, 500);
}

#[test]
fn test_parse_coordinates_imprecise_insertion_without_length() {
    // GIVEN an imprecise insertion with END == POS and no SVLEN
    let record = sv_record("N", "<INS>", "SVTYPE=INS;END=724779");

    // WHEN parsing the coordinates
    let coordinates = coordinates(&record, VariantCategory::Sv);

    // THEN the length is unknown
    assert_eq!(coordinates.sub_category, "ins");
    assert_eq!(coordinates.length, -1);
}

#[test]
fn test_parse_coordinates_mei() {
    // GIVEN a mobile element insertion
    let record = record(
        &["##contig=<ID=1,length=249250621>"],
        "1\t724779\t.\tA\tACTGCTG\t.\tPASS\t.",
    );

    // WHEN parsing the coordinates as an MEI
    let coordinates = coordinates(&record, VariantCategory::Mei);

    // THEN the sub category is mei and the length is the inserted length
    assert_eq!(coordinates.sub_category, "mei");
    assert_eq!(coordinates.length, 6);
}

#[test]
fn test_parse_coordinates_mei_same_length_alleles() {
    // GIVEN an MEI with equally long REF and ALT
    let record = record(
        &["##contig=<ID=1,length=249250621>"],
        "1\t724779\t.\tA\tG\t.\tPASS\t.",
    );

    // WHEN parsing the coordinates as an MEI
    let coordinates = coordinates(&record, VariantCategory::Mei);

    // THEN the length falls back to the ALT length
    assert_eq!(coordinates.sub_category, "mei");
    assert_eq!(coordinates.length, 1);
}

#[test]
fn test_parse_coordinates_cytobands() {
    // GIVEN a deletion spanning two cytobands
    let record = sv_record("N", "<DEL>", "SVTYPE=DEL;END=3000000");

    // WHEN parsing the coordinates with cytoband definitions
    let coordinates = parse_coordinates(
        &record,
        record.header(),
        &test_cytobands(),
        &VariantCategory::Sv,
    );

    // THEN start and end cytobands are annotated
    assert_eq!(coordinates.cytoband_start.as_deref(), Some("p36.33"));
    assert_eq!(coordinates.cytoband_end.as_deref(), Some("p36.32"));
}

#[test]
fn test_get_cytoband_coordinates() {
    let cytobands = test_cytobands();

    // Positions inside the bands, with and without chr prefix
    assert_eq!(
        get_cytoband_coordinates(&cytobands, "1", 1).as_deref(),
        Some("p36.33")
    );
    assert_eq!(
        get_cytoband_coordinates(&cytobands, "chr1", 2_300_000).as_deref(),
        Some("p36.32")
    );

    // Band ends are exclusive; outside the bands or chromosomes returns None
    assert_eq!(get_cytoband_coordinates(&cytobands, "1", 5_400_000), None);
    assert_eq!(get_cytoband_coordinates(&cytobands, "2", 100), None);
    assert_eq!(get_cytoband_coordinates(&cytobands, "1", -1), None);
}

// Ported from scout's `tests/parse/test_parse_coordinates.py`.
#[test]
fn test_get_sv_length_small_ins() {
    // GIVEN an insertion with the whole sequence in ALT (POS == END)
    assert_eq!(sv_length(144343218, 144343218, "1", "1", Some(296)), 296);
}

#[test]
fn test_get_sv_length_large_ins_no_length() {
    // GIVEN an imprecise insertion without SVLEN
    assert_eq!(sv_length(133920667, 133920667, "1", "1", None), -1);
}

#[test]
fn test_get_sv_length_translocation() {
    // GIVEN a translocation between different chromosomes
    assert_eq!(sv_length(726044, 0, "1", "10", None), i32::MAX);
}

#[test]
fn test_get_sv_length_translocation_same_chrom() {
    assert_eq!(
        sv_length(726044, 16457990, "1", "1", None),
        16457990 - 726044
    );
}

#[test]
fn test_get_sv_length_cnvnator_del() {
    // GIVEN a CNVnator deletion with negative SVLEN
    assert_eq!(sv_length(1, 10000, "1", "1", Some(-10000)), 10000);
}

#[test]
fn test_get_sv_length_del_no_length() {
    assert_eq!(sv_length(869314, 870246, "1", "1", None), 870246 - 869314);
}

#[test]
fn test_get_sv_length_end_zero() {
    assert_eq!(sv_length(869314, 0, "1", "1", None), -1);
}

#[test]
fn test_get_end_tiddit_translocation() {
    // GIVEN a TIDDIT translocation, END comes from the ALT field
    assert_eq!(sv_end(724779, "N[hs37d5:12060532[", None, None), 12060532);
}

#[test]
fn test_get_end_bnd_alt_overrides_svend() {
    assert_eq!(sv_end(724779, "]1:5000]N", Some(1), None), 5000);
}

#[test]
fn test_get_end_deletion() {
    assert_eq!(sv_end(869314, "<DEL>", Some(870246), None), 870246);
}

#[test]
fn test_get_end_sgl_bnd() {
    // GIVEN a single-end BND, END is set to POS
    assert_eq!(sv_end(123455, "GTCA.", None, None), 123455);
}

#[test]
fn test_get_end_from_svlen() {
    assert_eq!(sv_end(1000, "<DUP>", None, Some(250)), 1250);
}

#[test]
fn test_get_end_fallback_to_pos() {
    assert_eq!(sv_end(1000, "<INS>", None, None), 1000);
}

#[test]
fn test_get_end_chrom() {
    assert_eq!(get_end_chrom("N[hs37d5:12060532[", "1"), "hs37d5");
    assert_eq!(get_end_chrom("]CHR2:100]N", "1"), "2");
    assert_eq!(get_end_chrom("<DEL>", "1"), "1");
    assert_eq!(get_end_chrom("N:", "1"), "1");
}

#[test]
fn test_normalize_chromosome() {
    assert_eq!(normalize_chromosome("chr1"), "1");
    assert_eq!(normalize_chromosome("Chr1"), "1");
    assert_eq!(normalize_chromosome("CHRX"), "X");
    assert_eq!(normalize_chromosome("MT"), "MT");
    assert_eq!(normalize_chromosome("ch"), "ch");
}
