use scout_loader::parser::select_vcfs;

#[path = "fixtures/case.rs"]
mod case;
use case::test_case_config;

#[test]
fn test_select_vcfs_returns_all_categories() {
    let mut config = test_case_config();

    config.vcf_snv_mt = Some("snv_mt.vcf".into());
    config.vcf_str = Some("str.vcf".into());
    config.vcf_mei = Some("mei.vcf".into());
    config.vcf_fusion = Some("fusion.vcf".into());
    config.vcf_cancer = Some("cancer.vcf".into());
    config.vcf_cancer_sv = Some("cancer_sv.vcf".into());

    let selected = select_vcfs(&config, None, false);

    let paths: Vec<_> = selected.iter().map(|(path, _)| *path).collect();

    assert_eq!(
        paths,
        vec![
            config.vcf_snv.as_ref().unwrap(),
            config.vcf_snv_mt.as_ref().unwrap(),
            config.vcf_sv.as_ref().unwrap(),
            config.vcf_mei.as_ref().unwrap(),
            config.vcf_str.as_ref().unwrap(),
            config.vcf_fusion.as_ref().unwrap(),
            config.vcf_cancer.as_ref().unwrap(),
            config.vcf_cancer_sv.as_ref().unwrap(),
        ]
    );
}
