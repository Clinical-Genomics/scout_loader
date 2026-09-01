use scout_loader::models::case::CaseConfig;

pub fn test_case_config() -> CaseConfig {
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
