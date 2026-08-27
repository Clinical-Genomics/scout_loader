use mongodb::bson::doc;

use scout_loader::models::variant::VariantCategory;
use scout_loader::parse::vcf::should_load_variant;
use scout_loader::utils::hash::generate_md5_key;
use std::collections::HashSet;

#[test]
fn test_should_load_variant() {
    let threshold = 5;
    let managed_variant_ids = HashSet::new();

    let cases = [
        // No rank score.
        (
            doc! {
                "chromosome": "1",
            },
            VariantCategory::Snv,
            true,
        ),
        // Rank score above threshold.
        (
            doc! {
                "rank_score": 6,
                "chromosome": "1",
            },
            VariantCategory::Snv,
            true,
        ),
        // Rank score equal to threshold.
        (
            doc! {
                "rank_score": 5,
                "chromosome": "1",
            },
            VariantCategory::Snv,
            false,
        ),
        // Rank score below threshold.
        (
            doc! {
                "rank_score": 4,
                "chromosome": "1",
            },
            VariantCategory::Snv,
            false,
        ),
        // Mitochondrial variant.
        (
            doc! {
                "rank_score": 1,
                "chromosome": "MT",
            },
            VariantCategory::Snv,
            true,
        ),
        // STR variant.
        (
            doc! {
                "rank_score": 1,
                "chromosome": "1",
            },
            VariantCategory::Str,
            true,
        ),
        // Pathogenic variant.
        (
            doc! {
                "rank_score": 1,
                "chromosome": "1",
                "clnsig": [
                    {
                        "value": "pathogenic",
                    }
                ],
            },
            VariantCategory::Snv,
            true,
        ),
    ];

    for (variant, category, expected) in cases {
        assert_eq!(
            should_load_variant(&variant, category, threshold, &managed_variant_ids),
            expected,
        );
    }
}

#[test]
fn should_load_managed_variant() {
    let variant = doc! {
        "chromosome": "5",
        "position": 112043220,
        "reference": "A",
        "alternative": "C",
        "rank_score": 0,
    };

    let managed_variant_id = generate_md5_key(&[
        "5".to_string(),
        "112043220".to_string(),
        "A".to_string(),
        "C".to_string(),
        "clinical".to_string(),
    ]);

    let managed_variant_ids = HashSet::from([managed_variant_id]);

    assert!(should_load_variant(
        &variant,
        VariantCategory::Snv,
        5,
        &managed_variant_ids,
    ));
}
