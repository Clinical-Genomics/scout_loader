use mongodb::bson::doc;

use scout_loader::models::variant::VariantCategory;
use scout_loader::parse::vcf::should_load_variant;

#[test]
fn test_should_load_variant() {
    let threshold = 5;

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
        assert_eq!(should_load_variant(&variant, category, threshold), expected);
    }
}
