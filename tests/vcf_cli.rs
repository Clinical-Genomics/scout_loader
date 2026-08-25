mod fixtures;

use assert_cmd::Command;
use fixtures::{TestDatabase, fixture_path};
use predicates::prelude::*;

#[tokio::test]
async fn cli_processes_minimal_case() {
    let test_db = TestDatabase::new("minimal_case").await;
    let config_path = test_db.config_path();

    Command::cargo_bin("scout_loader")
        .unwrap()
        .env("TEST_ENV", "1")
        .args([
            "--config",
            config_path.to_str().unwrap(),
            "--case-config",
            fixture_path("minimal_case.yaml").to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Snv: 1 variants added"))
        .stdout(predicate::str::contains("Sv: 1 variants added"))
        .stdout(predicate::str::contains(
            "Total variants added for case case_123: 2",
        ));

    assert_eq!(test_db.count_variants().await, 2);

    test_db.cleanup().await;
}

#[tokio::test]
async fn cli_processes_minimal_vep_case() {
    let test_db = TestDatabase::new("minimal_vep_case").await;
    let config_path = test_db.config_path();

    Command::cargo_bin("scout_loader")
        .unwrap()
        .env("TEST_ENV", "1")
        .args([
            "--config",
            config_path.to_str().unwrap(),
            "--case-config",
            fixture_path("minimal_vep_case.yaml").to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("Snv: 1 variants added"))
        .stdout(predicate::str::contains(
            "Total variants added for case case_123: 1",
        ));

    assert_eq!(test_db.count_variants().await, 1);

    test_db.cleanup().await;
}
