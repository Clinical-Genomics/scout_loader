use std::collections::HashMap;
use std::sync::LazyLock;

pub struct SoTerm {
    pub rank: u32,
    pub region: &'static str,
}

pub static SO_TERMS: LazyLock<HashMap<&'static str, SoTerm>> =
    LazyLock::new(|| {
        HashMap::from([
            (
                "transcript_ablation",
                SoTerm { rank: 1, region: "exonic" },
            ),
            (
                "splice_donor_variant",
                SoTerm { rank: 2, region: "splicing" },
            ),
            (
                "splice_acceptor_variant",
                SoTerm { rank: 3, region: "splicing" },
            ),
            (
                "stop_gained",
                SoTerm { rank: 4, region: "exonic" },
            ),
            (
                "frameshift_variant",
                SoTerm { rank: 5, region: "exonic" },
            ),
            (
                "stop_lost",
                SoTerm { rank: 6, region: "exonic" },
            ),
            (
                "start_lost",
                SoTerm { rank: 7, region: "exonic" },
            ),
            (
                "missense_variant",
                SoTerm { rank: 11, region: "exonic" },
            ),
            (
                "synonymous_variant",
                SoTerm { rank: 20, region: "exonic" },
            ),
            (
                "intron_variant",
                SoTerm { rank: 31, region: "intronic" },
            ),
            (
                "upstream_gene_variant",
                SoTerm { rank: 34, region: "upstream" },
            ),
            (
                "downstream_gene_variant",
                SoTerm { rank: 35, region: "downstream" },
            ),
            (
                "intergenic_variant",
                SoTerm { rank: 43, region: "intergenic_variant" },
            ),
            (
                "sequence_variant",
                SoTerm { rank: 44, region: "genomic_feature" },
            ),
        ])
    });
