mod fixtures;

use assert_cmd::Command;
use fixtures::{fixture_path, TestDatabase};
use predicates::prelude::*;

fn run_with_config(case_config: &str, db_config: &str) -> Command {
    let mut cmd =
        Command::cargo_bin("scout_loader").expect("binary should build");

    cmd.current_dir(env!("CARGO_MANIFEST_DIR"))
        .env("TEST_ENV", "1")
        .arg("--config")
        .arg(fixture_path(db_config))
        .arg("--case-config")
        .arg(fixture_path(case_config));

    cmd
}

#[tokio::test]
async fn cli_processes_minimal_snv_vcf() {
    let test_db = TestDatabase::new().await;

    let mut cmd = run_with_config("minimal_case.yaml", "test_config.toml");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123456_A_T_clinical"));

    assert_eq!(test_db.count_variants().await, 1);
}

#[tokio::test]
async fn cli_processes_minimal_vep_snv_vcf() {
    let test_db = TestDatabase::new().await;

    let mut cmd =
        run_with_config("minimal_vep_case.yaml", "test_config.toml");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123460_G_A_clinical"))
        .stdout(predicate::str::contains("genes"))
        .stdout(predicate::str::contains("hgnc_ids"))
        .stdout(predicate::str::contains("1101"))
        .stdout(predicate::str::contains("missense_variant"));

    assert_eq!(test_db.count_variants().await, 1);
}

#[tokio::test]
async fn cli_processes_minimal_sv_vcf() {
    let test_db = TestDatabase::new().await;

    let mut cmd = run_with_config("minimal_case.yaml", "test_config.toml");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123500_N_<DEL>_clinical"))
        .stdout(predicate::str::contains("sub_category"))
        .stdout(predicate::str::contains("del"))
        .stdout(predicate::str::contains("end"))
        .stdout(predicate::str::contains("123650"))
        .stdout(predicate::str::contains("sample_id"))
        .stdout(predicate::str::contains("0/1"));

    assert_eq!(test_db.count_variants().await, 1);
}

#[tokio::test]
async fn cli_processes_multiple_vcfs() {
    let test_db = TestDatabase::new().await;

    let mut cmd = run_with_config("minimal_case.yaml", "test_config.toml");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("1_123456_A_T_clinical"))
        .stdout(predicate::str::contains("1_123500_N_<DEL>_clinical"));

    assert_eq!(test_db.count_variants().await, 2);
}