#[test]
fn test_version() {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_scout_loader"))
        .arg("--version")
        .output()
        .expect("failed to run scout_loader");

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        format!("scout_loader {}", env!("CARGO_PKG_VERSION"))
    );
}
