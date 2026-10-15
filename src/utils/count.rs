use std::error::Error;
use std::io;
use std::process::{Command, Stdio};

/// Counts lines in a VCF file using `wc -l`.
///
/// Supports compressed and uncompressed files. For compressed files,
/// the contents are decompressed before counting.
pub fn count_vcf_records(path: &str) -> Result<usize, Box<dyn Error>> {
    let output = if path.ends_with(".gz") {
        let mut decompressor = Command::new("gzip")
            .args(["-cd", path])
            .stdout(Stdio::piped())
            .spawn()?;

        let output = Command::new("wc")
            .arg("-l")
            .stdin(
                decompressor
                    .stdout
                    .take()
                    .ok_or_else(|| io::Error::other("Failed to capture decompressor output"))?,
            )
            .output()?;

        let status = decompressor.wait()?;
        if !status.success() || !output.status.success() {
            return Err(io::Error::other("Failed to count VCF lines").into());
        }

        output
    } else {
        Command::new("wc").arg("-l").arg(path).output()?
    };

    if !output.status.success() {
        return Err(io::Error::other("Failed to count VCF lines").into());
    }

    let count = String::from_utf8(output.stdout)?
        .split_whitespace()
        .next()
        .ok_or_else(|| io::Error::other("Failed to parse VCF line count"))?
        .parse::<usize>()?;
    Ok(count)
}
