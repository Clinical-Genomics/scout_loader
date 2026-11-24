//! Caller parsing tests, ported from scout's `tests/parse/test_parse_callers.py`.
//!
//! Unlike scout, which returns every relevant caller with `None` for callers
//! that did not make a call, scout_loader drops those callers from the
//! document. "Caller is None" in scout therefore maps to "key is absent" here.

use crate::vcf_record::record;
use mongodb::bson::{Document, doc};
use scout_loader::models::variant::VariantCategory;
use scout_loader::parse::callers::parse_callers;
use scout_loader::parse::filters::parse_filters;

const META: &[&str] = &[
    "##contig=<ID=1,length=249250621>",
    "##FILTER=<ID=MQRankSum,Description=\"MQRankSum filter\">",
    "##FILTER=<ID=MinSomaticScore,Description=\"Somatic score filter\">",
    "##FILTER=<ID=LowQual,Description=\"Low quality\">",
    "##INFO=<ID=FOUND_IN,Number=1,Type=String,Description=\"Program that called the variant\">",
    "##INFO=<ID=svdb_origin,Number=1,Type=String,Description=\"Callers merged by SVDB\">",
    "##INFO=<ID=set,Number=1,Type=String,Description=\"Source VCF for the merged record\">",
];

/// Parse callers for a variant with the given FILTER and INFO columns.
fn callers(category: VariantCategory, filter: &str, info: &str) -> Document {
    let record = record(META, &format!("1\t880086\t.\tT\tC\t.\t{filter}\t{info}"));
    let filters = parse_filters(&record, record.header());

    parse_callers(&record, category, &filters)
}

fn all_with_status(callers: &[&str], status: &str) -> Document {
    callers
        .iter()
        .map(|caller| (caller.to_string(), status.into()))
        .collect()
}

#[test]
fn test_parse_callers() {
    // GIVEN information that gatk and freebayes have passed
    let callers = callers(VariantCategory::Snv, "PASS", "set=gatk-freebayes");

    // THEN only gatk and freebayes are passed
    assert_eq!(callers, doc! { "freebayes": "Pass", "gatk": "Pass" });
}

#[test]
fn test_parse_callers_only_one() {
    // GIVEN information that only gatk called the variant
    let callers = callers(VariantCategory::Snv, "PASS", "set=gatk");

    // THEN only gatk is passed
    assert_eq!(callers, doc! { "gatk": "Pass" });
}

#[test]
fn test_parse_callers_complex() {
    // GIVEN a set where samtools filtered the variant
    let callers = callers(
        VariantCategory::Snv,
        "PASS",
        "set=gatk-filterInsamtools-freebayes",
    );

    // THEN samtools is filtered and the others pass
    assert_eq!(
        callers,
        doc! { "freebayes": "Pass", "gatk": "Pass", "samtools": "Filtered" }
    );
}

#[test]
fn test_parse_callers_complex_filter() {
    // GIVEN a FILTERed variant with a single filterIn caller
    let callers = callers(VariantCategory::Snv, "MQRankSum", "set=filterInsamtools");

    // THEN that caller gets the specific filter status
    assert_eq!(callers, doc! { "samtools": "Filtered - MQRankSum" });
}

#[test]
fn test_parse_cancer_callers_filterin_tnscope_umi() {
    // GIVEN a set where tnscope_umi filtered the variant
    let callers = callers(VariantCategory::Cancer, "PASS", "set=filterIntnscope_umi");

    // THEN only tnscope_umi is filtered
    assert_eq!(callers, doc! { "tnscope_umi": "Filtered" });
}

#[test]
fn test_parse_callers_intersection() {
    // GIVEN information that all callers agree on Pass
    let callers = callers(VariantCategory::Snv, "PASS", "set=Intersection");

    // THEN all callers pass
    assert_eq!(callers.get_str("gatk"), Ok("Pass"));
    assert_eq!(callers.get_str("freebayes"), Ok("Pass"));
    assert_eq!(callers.get_str("samtools"), Ok("Pass"));
}

#[test]
fn test_parse_callers_intersection_svdb_info() {
    // GIVEN both svdb_origin and set annotations
    let callers = callers(
        VariantCategory::Snv,
        "PASS",
        "svdb_origin=gatk|deepvariant;set=Intersection",
    );

    // THEN svdb_origin takes precedence and only its callers are passed
    assert_eq!(callers, doc! { "deepvariant": "Pass", "gatk": "Pass" });
}

#[test]
fn test_parse_callers_filtered_svdb_info() {
    // GIVEN svdb_origin information that manta made a filtered call
    let callers = callers(
        VariantCategory::CancerSv,
        "MinSomaticScore",
        "svdb_origin=manta",
    );

    // THEN the manta filter status is set
    assert_eq!(callers, doc! { "manta": "Filtered - MinSomaticScore" });
}

#[test]
fn test_parse_callers_filtered_svdb_info_multiple_callers() {
    // GIVEN a filtered variant merged from several callers by SVDB
    let callers = callers(
        VariantCategory::Sv,
        "LowQual",
        "svdb_origin=manta|tiddit|unknown_caller",
    );

    // THEN we cannot tell which caller filtered it, so all are "Filtered",
    // and callers not relevant to the category are ignored
    assert_eq!(callers, doc! { "manta": "Filtered", "tiddit": "Filtered" });
}

#[test]
fn test_parse_callers_filtered_found_in_info() {
    // GIVEN FOUND_IN information that manta made a filtered call
    let callers = callers(
        VariantCategory::CancerSv,
        "MinSomaticScore",
        "FOUND_IN=manta",
    );

    // THEN the manta filter status is set
    assert_eq!(callers, doc! { "manta": "Filtered - MinSomaticScore" });
}

#[test]
fn test_parse_callers_found_in_multiple_filters() {
    // GIVEN a single caller and several FILTER values
    let callers = callers(
        VariantCategory::Snv,
        "MQRankSum;LowQual",
        "FOUND_IN=deepvariant",
    );

    // THEN all filters are listed in the caller status
    assert_eq!(
        callers,
        doc! { "deepvariant": "Filtered - MQRankSum - LowQual" }
    );
}

#[test]
fn test_parse_callers_found_in_multiple_callers() {
    // GIVEN FOUND_IN with several callers, one with a pipe-separated suffix
    let callers = callers(
        VariantCategory::Snv,
        "PASS",
        "FOUND_IN=deepvariant|extra,gatk",
    );

    // THEN the caller id before the pipe is used and both pass
    assert_eq!(callers, doc! { "deepvariant": "Pass", "gatk": "Pass" });
}

#[test]
fn test_parse_callers_found_in_multiple_callers_filtered() {
    // GIVEN a filtered variant found by several callers
    let callers = callers(
        VariantCategory::Snv,
        "MQRankSum",
        "FOUND_IN=deepvariant,gatk",
    );

    // THEN all are marked as generically "Filtered"
    assert_eq!(
        callers,
        doc! { "deepvariant": "Filtered", "gatk": "Filtered" }
    );
}

#[test]
fn test_parse_callers_found_in_takes_precedence() {
    // GIVEN FOUND_IN, svdb_origin and set annotations
    let callers = callers(
        VariantCategory::Snv,
        "PASS",
        "FOUND_IN=deepvariant;svdb_origin=gatk;set=Intersection",
    );

    // THEN only FOUND_IN is used
    assert_eq!(callers, doc! { "deepvariant": "Pass" });
}

#[test]
fn test_parse_callers_filtered_all() {
    // GIVEN information that all callers agree on filtered
    let callers = callers(VariantCategory::Snv, "PASS", "set=FilteredInAll");

    // THEN all callers are filtered
    assert_eq!(callers.get_str("gatk"), Ok("Filtered"));
    assert_eq!(callers.get_str("freebayes"), Ok("Filtered"));
    assert_eq!(callers.get_str("samtools"), Ok("Filtered"));
}

#[test]
fn test_parse_sv_callers_filtered_all() {
    // GIVEN information that all SV callers agree on filtered
    let callers = callers(VariantCategory::Sv, "PASS", "set=FilteredInAll");

    // THEN all SV callers are filtered
    assert_eq!(callers.get_str("cnvnator"), Ok("Filtered"));
    assert_eq!(callers.get_str("delly"), Ok("Filtered"));
    assert_eq!(callers.get_str("tiddit"), Ok("Filtered"));
    assert_eq!(callers.get_str("manta"), Ok("Filtered"));
}

#[test]
fn test_parse_sv_callers_intersection() {
    // GIVEN information that all SV callers agree on Pass
    let callers = callers(VariantCategory::Sv, "PASS", "set=Intersection");

    // THEN all SV callers pass
    assert_eq!(callers.get_str("cnvnator"), Ok("Pass"));
    assert_eq!(callers.get_str("delly"), Ok("Pass"));
    assert_eq!(callers.get_str("tiddit"), Ok("Pass"));
    assert_eq!(callers.get_str("manta"), Ok("Pass"));
}

#[test]
fn test_parse_sv_callers_filterin_tiddit() {
    // GIVEN a set where manta passed and tiddit filtered
    let callers = callers(VariantCategory::Sv, "PASS", "set=manta-filterIntiddit");

    // THEN only manta and tiddit have a status
    assert_eq!(callers, doc! { "manta": "Pass", "tiddit": "Filtered" });
}

#[test]
fn test_parse_callers_found_in_deepvariant() {
    // GIVEN a VEP-annotated SNV found by deepvariant (scout: one_vep104_annotated_variant)
    let callers = callers(VariantCategory::Snv, "PASS", "FOUND_IN=deepvariant");

    // THEN deepvariant is found, and the GATK fallback is not applied
    assert_eq!(callers, doc! { "deepvariant": "Pass" });
}

#[test]
fn test_parse_callers_gatk_snv_fallback() {
    // GIVEN an SNV without any caller annotation
    let passed = callers(VariantCategory::Snv, "PASS", ".");
    let filtered = callers(VariantCategory::Snv, "MQRankSum;LowQual", ".");

    // THEN it is assumed to be a GATK call (older MIP versions)
    assert_eq!(passed, doc! { "gatk": "Pass" });
    assert_eq!(filtered, doc! { "gatk": "Filtered - MQRankSum - LowQual" });
}

#[test]
fn test_parse_callers_no_fallback_for_other_categories() {
    // GIVEN an SV without any caller annotation
    let callers = callers(VariantCategory::Sv, "PASS", ".");

    // THEN no callers are set
    assert_eq!(callers, Document::new());
}

#[test]
fn test_parse_callers_all() {
    // GIVEN all callers called a variant, for every category
    let cases: &[(VariantCategory, &[&str])] = &[
        (
            VariantCategory::Snv,
            &[
                "bcftools",
                "deepvariant",
                "freebayes",
                "gatk",
                "mutect2",
                "samtools",
                "mitorsaw",
            ],
        ),
        (
            VariantCategory::Sv,
            &[
                "cnvnator",
                "cnvpytor",
                "delly",
                "gatk",
                "gcnvcaller",
                "hificnv",
                "manta",
                "mitosalt",
                "mitorsaw",
                "severus",
                "sniffles",
                "tiddit",
            ],
        ),
        (
            VariantCategory::Cancer,
            &[
                "freebayes",
                "gatk",
                "mutect",
                "pindel",
                "tnscope",
                "tnscope_umi",
                "vardict",
            ],
        ),
        (
            VariantCategory::CancerSv,
            &[
                "ascat", "cnvkit", "dellysv", "dellycnv", "gatk", "igh_dux4", "manta", "tiddit",
            ],
        ),
        (VariantCategory::Mei, &["retroseq"]),
        (VariantCategory::Str, &["expansionhunter", "trgt"]),
        (
            VariantCategory::Fusion,
            &["arriba", "fusioncatcher", "starfusion"],
        ),
    ];

    for (category, expected) in cases {
        // WHEN parsing callers for an Intersection set
        let parsed = callers(*category, "PASS", "set=Intersection");

        // THEN exactly the callers of that category are passed
        assert_eq!(
            parsed,
            all_with_status(expected, "Pass"),
            "category {category:?}"
        );
    }
}
