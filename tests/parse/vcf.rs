use mongodb::bson::doc;
use scout_loader::models::case::CaseConfig;
use scout_loader::models::variant::VariantCategory;
use scout_loader::parse::vcf::should_load_variant;
use scout_loader::parser::select_vcfs;
use scout_loader::utils::hash::generate_md5_key;
use std::collections::HashSet;

fn test_variant() -> mongodb::bson::Document {
    doc! {
        "chromosome": "1",
        "position": 100_i64,
        "reference": "A",
        "alternative": "G",
        "simple_id": "1_100_A_G",
    }
}

fn test_case_config() -> CaseConfig {
    CaseConfig {
        owner: "test_institute".to_string(),
        family: "case_123".to_string(),
        human_genome_build: "37".to_string(),
        gene_panels: None,
        samples: vec![],
        rank_score_threshold: Some(-100),
        vcf_snv: Some("tests/fixtures/minimal_snv.vcf".into()),
        vcf_sv: Some("tests/fixtures/minimal_sv.vcf".into()),
        vcf_str: None,
        vcf_mei: None,
        vcf_cancer: None,
        vcf_cancer_sv: None,
        vcf_fusion: None,
        vcf_snv_research: Some("tests/fixtures/minimal_snv_research.vcf".into()),
        vcf_sv_research: Some("tests/fixtures/minimal_sv_research.vcf".into()),
        vcf_str_research: None,
        vcf_mei_research: None,
        vcf_cancer_research: None,
        vcf_cancer_sv_research: None,
        vcf_fusion_research: None,
        custom_images: None,
    }
}

#[test]
fn test_should_load_variant() {
    let threshold = 5;
    let managed_variant_ids = HashSet::new();
    let causative_variant_ids = HashSet::new();

    let cases = [
        // No rank score.
        (test_variant(), VariantCategory::Snv, true),
        // Rank score above threshold.
        (
            {
                let mut variant = test_variant();
                variant.insert("rank_score", 6);
                variant
            },
            VariantCategory::Snv,
            true,
        ),
        // Rank score equal to threshold.
        (
            {
                let mut variant = test_variant();
                variant.insert("rank_score", 5);
                variant
            },
            VariantCategory::Snv,
            false,
        ),
        // Rank score below threshold.
        (
            {
                let mut variant = test_variant();
                variant.insert("rank_score", 4);
                variant
            },
            VariantCategory::Snv,
            false,
        ),
        // Mitochondrial variant.
        (
            {
                let mut variant = test_variant();
                variant.insert("chromosome", "MT");
                variant.insert("rank_score", 1);
                variant
            },
            VariantCategory::Snv,
            true,
        ),
        // STR variant.
        (
            {
                let mut variant = test_variant();
                variant.insert("rank_score", 1);
                variant
            },
            VariantCategory::Str,
            true,
        ),
        // Pathogenic variant.
        (
            {
                let mut variant = test_variant();
                variant.insert("rank_score", 1);
                variant.insert(
                    "clnsig",
                    vec![doc! {
                        "value": "pathogenic",
                    }],
                );
                variant
            },
            VariantCategory::Snv,
            true,
        ),
    ];

    for (variant, category, expected) in cases {
        assert_eq!(
            should_load_variant(
                &variant,
                category,
                threshold,
                &managed_variant_ids,
                &causative_variant_ids,
            ),
            expected,
        );
    }
}

#[test]
fn should_load_managed_variant() {
    let variant = test_variant();

    let managed_variant_id = generate_md5_key(&[
        "1".to_string(),
        "100".to_string(),
        "A".to_string(),
        "G".to_string(),
        "clinical".to_string(),
    ]);

    let managed_variant_ids = HashSet::from([managed_variant_id]);
    let causative_variant_ids = HashSet::new();

    assert!(should_load_variant(
        &variant,
        VariantCategory::Snv,
        5,
        &managed_variant_ids,
        &causative_variant_ids,
    ));
}

#[test]
fn should_load_causative_variant() {
    let variant = test_variant();
    let managed_variant_ids = HashSet::new();

    let causative_variant_ids = HashSet::from(["1_100_A_G_clinical".to_string()]);

    assert!(should_load_variant(
        &variant,
        VariantCategory::Snv,
        5,
        &managed_variant_ids,
        &causative_variant_ids,
    ));
}

#[test]
fn test_selects_all_clinical_vcfs_by_default() {
    let config = test_case_config();

    let selected = select_vcfs(&config, None, false);

    assert_eq!(selected.len(), 2);
    assert_eq!(selected[0].1, VariantCategory::Snv);
    assert_eq!(selected[1].1, VariantCategory::Sv);
}

#[test]
fn test_selects_all_research_vcfs() {
    let config = test_case_config();

    let selected = select_vcfs(&config, None, true);

    assert_eq!(selected.len(), 2);
    assert_eq!(selected[0].1, VariantCategory::Snv);
    assert_eq!(selected[1].1, VariantCategory::Sv);
}

#[test]
fn test_selects_requested_clinical_categories() {
    let config = test_case_config();
    let categories = vec!["snv".to_string()];

    let selected = select_vcfs(&config, Some(&categories), false);

    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].1, VariantCategory::Snv);
}

#[test]
fn test_selects_requested_research_categories() {
    let config = test_case_config();
    let categories = vec!["sv".to_string()];

    let selected = select_vcfs(&config, Some(&categories), true);

    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].1, VariantCategory::Sv);
}
