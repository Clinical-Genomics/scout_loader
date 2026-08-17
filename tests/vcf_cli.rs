use assert_cmd::Command;
use predicates::prelude::*;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixture_path(name: &str) -> PathBuf {
    repo_root().join("tests").join("fixtures").join(name)
}

#[test]
fn cli_processes_minimal_snv_vcf() {
    let fixture = fixture_path("minimal_snv.vcf");

    let mut cmd = Command::cargo_bin("scout_loader").expect("binary should build");

    cmd.current_dir(repo_root())
        .arg("--snv")
        .arg(&fixture)
        .arg("--variant-type")
        .arg("clinical")
        .arg("--case-id")
        .arg("case_123")
        .arg("--genome-build")
        .arg("GRCh37")
        .arg("--samples")
        .arg("SAMPLE1:sample-one:0");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("display_name"))
        .stdout(predicate::str::contains("1_123456_A_T_clinical"))
        .stdout(predicate::str::contains("sample_id"))
        .stdout(predicate::str::contains("SAMPLE1"))
        .stdout(predicate::str::contains("genotype_call"))
        .stdout(predicate::str::contains("0/1"));
}

#[test]
fn cli_processes_minimal_vep_snv_vcf() {
    let fixture = fixture_path("minimal_vep_snv.vcf");

    let mut cmd = Command::cargo_bin("scout_loader").expect("binary should build");

    cmd.current_dir(repo_root())
        .arg("--snv")
        .arg(&fixture)
        .arg("--variant-type")
        .arg("clinical")
        .arg("--case-id")
        .arg("case_123")
        .arg("--genome-build")
        .arg("GRCh37")
        .arg("--samples")
        .arg("SAMPLE1:sample-one:0");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123460_G_A_clinical"))
        .stdout(predicate::str::contains("genes"))
        .stdout(predicate::str::contains("BRCA2"))
        .stdout(predicate::str::contains("hgnc_ids"))
        .stdout(predicate::str::contains("1101"))
        .stdout(predicate::str::contains("missense_variant"));
}

#[test]
fn cli_processes_minimal_sv_vcf() {
    let fixture = fixture_path("minimal_sv.vcf");

    let mut cmd = Command::cargo_bin("scout_loader").expect("binary should build");

    cmd.current_dir(repo_root())
        .arg("--sv")
        .arg(&fixture)
        .arg("--variant-type")
        .arg("clinical")
        .arg("--case-id")
        .arg("case_123")
        .arg("--genome-build")
        .arg("GRCh37")
        .arg("--samples")
        .arg("SAMPLE1:sample-one:0");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123500_N_<DEL>_clinical"))
        .stdout(predicate::str::contains("sub_category"))
        .stdout(predicate::str::contains("del"))
        .stdout(predicate::str::contains("end"))
        .stdout(predicate::str::contains("123650"))
        .stdout(predicate::str::contains("sample_id"))
        .stdout(predicate::str::contains("0/1"));
}

#[test]
fn cli_processes_multiple_vcfs() {
    let snv_fixture = fixture_path("minimal_snv.vcf");
    let sv_fixture = fixture_path("minimal_sv.vcf");

    let mut cmd = Command::cargo_bin("scout_loader").expect("binary should build");

    cmd.current_dir(repo_root())
        .arg("--snv")
        .arg(&snv_fixture)
        .arg("--sv")
        .arg(&sv_fixture)
        .arg("--variant-type")
        .arg("clinical")
        .arg("--case-id")
        .arg("case_123")
        .arg("--genome-build")
        .arg("GRCh37")
        .arg("--samples")
        .arg("SAMPLE1:sample-one:0");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123456_A_T_clinical"))
        .stdout(predicate::str::contains("1_123500_N_<DEL>_clinical"));
}
