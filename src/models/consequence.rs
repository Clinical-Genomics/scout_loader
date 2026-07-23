use std::collections::HashMap;
use std::sync::LazyLock;

pub struct SoTerm {
    pub rank: u32,
    pub region: &'static str,
}

pub static SO_TERMS: LazyLock<HashMap<&'static str, SoTerm>> =
    LazyLock::new(|| {
        HashMap::from([
            ("transcript_ablation", SoTerm { rank: 1, region: "exonic" }),
            ("splice_donor_variant", SoTerm { rank: 2, region: "splicing" }),
            ("splice_acceptor_variant", SoTerm { rank: 3, region: "splicing" }),
            ("stop_gained", SoTerm { rank: 4, region: "exonic" }),
            ("frameshift_variant", SoTerm { rank: 5, region: "exonic" }),
            ("stop_lost", SoTerm { rank: 6, region: "exonic" }),
            ("start_lost", SoTerm { rank: 7, region: "exonic" }),
            ("initiator_codon_variant", SoTerm { rank: 8, region: "exonic" }),
            ("inframe_insertion", SoTerm { rank: 9, region: "exonic" }),
            ("inframe_deletion", SoTerm { rank: 10, region: "exonic" }),
            ("missense_variant", SoTerm { rank: 11, region: "exonic" }),
            ("protein_altering_variant", SoTerm { rank: 12, region: "exonic" }),
            ("transcript_amplification", SoTerm { rank: 13, region: "exonic" }),
            ("regulatory_region_ablation", SoTerm { rank: 14, region: "regulatory_region" }),
            ("splice_region_variant", SoTerm { rank: 15, region: "splicing" }),
            ("splice_donor_5th_base_variant", SoTerm { rank: 16, region: "splicing" }),
            ("splice_donor_region_variant", SoTerm { rank: 17, region: "splicing" }),
            ("splice_polypyrimidine_tract_variant", SoTerm { rank: 18, region: "splicing" }),
            ("incomplete_terminal_codon_variant", SoTerm { rank: 19, region: "exonic" }),
            ("synonymous_variant", SoTerm { rank: 20, region: "exonic" }),
            ("start_retained_variant", SoTerm { rank: 21, region: "exonic" }),
            ("stop_retained_variant", SoTerm { rank: 22, region: "exonic" }),
            ("coding_sequence_variant", SoTerm { rank: 23, region: "exonic" }),
            ("mature_miRNA_variant", SoTerm { rank: 24, region: "ncRNA_exonic" }),
            ("5_prime_UTR_variant", SoTerm { rank: 25, region: "5UTR" }),
            ("3_prime_UTR_variant", SoTerm { rank: 26, region: "3UTR" }),
            ("non_coding_transcript_exon_variant", SoTerm { rank: 27, region: "ncRNA_exonic" }),
            ("non_coding_exon_variant", SoTerm { rank: 28, region: "ncRNA_exonic" }),
            ("non_coding_transcript_variant", SoTerm { rank: 29, region: "ncRNA_exonic" }),
            ("nc_transcript_variant", SoTerm { rank: 30, region: "ncRNA_exonic" }),
            ("intron_variant", SoTerm { rank: 31, region: "intronic" }),
            ("NMD_transcript_variant", SoTerm { rank: 32, region: "ncRNA" }),
            ("coding_transcript_variant", SoTerm { rank: 33, region: "exonic" }),
            ("upstream_gene_variant", SoTerm { rank: 34, region: "upstream" }),
            ("downstream_gene_variant", SoTerm { rank: 35, region: "downstream" }),
            ("TFBS_ablation", SoTerm { rank: 36, region: "TFBS" }),
            ("TFBS_amplification", SoTerm { rank: 37, region: "TFBS" }),
            ("TF_binding_site_variant", SoTerm { rank: 38, region: "TFBS" }),
            ("regulatory_region_amplification", SoTerm { rank: 39, region: "regulatory_region" }),
            ("regulatory_region_variant", SoTerm { rank: 40, region: "regulatory_region" }),
            ("feature_elongation", SoTerm { rank: 41, region: "genomic_feature" }),
            ("feature_truncation", SoTerm { rank: 42, region: "genomic_feature" }),
            ("intergenic_variant", SoTerm { rank: 43, region: "intergenic_variant" }),
            ("sequence_variant", SoTerm { rank: 44, region: "genomic_feature" }),
        ])
    });
