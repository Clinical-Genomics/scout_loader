use std::collections::HashMap;
use std::sync::LazyLock;

pub struct SoTerm {
    pub region: &'static str,
}

pub static SO_TERMS: LazyLock<HashMap<&'static str, SoTerm>> =
    LazyLock::new(|| {
        HashMap::from([
            ("transcript_ablation", SoTerm { region: "exonic" }),
            ("splice_donor_variant", SoTerm { region: "splicing" }),
            ("splice_acceptor_variant", SoTerm { region: "splicing" }),
            ("stop_gained", SoTerm { region: "exonic" }),
            ("frameshift_variant", SoTerm { region: "exonic" }),
            ("stop_lost", SoTerm { region: "exonic" }),
            ("start_lost", SoTerm { region: "exonic" }),
            ("initiator_codon_variant", SoTerm { region: "exonic" }),
            ("inframe_insertion", SoTerm { region: "exonic" }),
            ("inframe_deletion", SoTerm { region: "exonic" }),
            ("missense_variant", SoTerm { region: "exonic" }),
            ("protein_altering_variant", SoTerm { region: "exonic" }),
            ("transcript_amplification", SoTerm { region: "exonic" }),
            ("regulatory_region_ablation", SoTerm { region: "regulatory_region" }),
            ("splice_region_variant", SoTerm { region: "splicing" }),
            ("splice_donor_5th_base_variant", SoTerm { region: "splicing" }),
            ("splice_donor_region_variant", SoTerm { region: "splicing" }),
            ("splice_polypyrimidine_tract_variant", SoTerm { region: "splicing" }),
            ("incomplete_terminal_codon_variant", SoTerm { region: "exonic" }),
            ("synonymous_variant", SoTerm { region: "exonic" }),
            ("start_retained_variant", SoTerm { region: "exonic" }),
            ("stop_retained_variant", SoTerm { region: "exonic" }),
            ("coding_sequence_variant", SoTerm { region: "exonic" }),
            ("mature_miRNA_variant", SoTerm { region: "ncRNA_exonic" }),
            ("5_prime_UTR_variant", SoTerm { region: "5UTR" }),
            ("3_prime_UTR_variant", SoTerm { region: "3UTR" }),
            ("non_coding_transcript_exon_variant", SoTerm { region: "ncRNA_exonic" }),
            ("non_coding_exon_variant", SoTerm { region: "ncRNA_exonic" }),
            ("non_coding_transcript_variant", SoTerm { region: "ncRNA_exonic" }),
            ("nc_transcript_variant", SoTerm { region: "ncRNA_exonic" }),
            ("intron_variant", SoTerm { region: "intronic" }),
            ("NMD_transcript_variant", SoTerm { region: "ncRNA" }),
            ("coding_transcript_variant", SoTerm { region: "exonic" }),
            ("upstream_gene_variant", SoTerm { region: "upstream" }),
            ("downstream_gene_variant", SoTerm { region: "downstream" }),
            ("TFBS_ablation", SoTerm { region: "TFBS" }),
            ("TFBS_amplification", SoTerm { region: "TFBS" }),
            ("TF_binding_site_variant", SoTerm { region: "TFBS" }),
            ("regulatory_region_amplification", SoTerm { region: "regulatory_region" }),
            ("regulatory_region_variant", SoTerm { region: "regulatory_region" }),
            ("feature_elongation", SoTerm { region: "genomic_feature" }),
            ("feature_truncation", SoTerm { region: "genomic_feature" }),
            ("intergenic_variant", SoTerm { region: "intergenic_variant" }),
            ("sequence_variant", SoTerm { region: "genomic_feature" }),
        ])
    });
