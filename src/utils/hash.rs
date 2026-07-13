/// Generates an MD5 hash from a list of string arguments.
///
/// The arguments are concatenated with a single space before computing
/// the hash, matching the behavior of the original Scout implementation.
///
/// # Arguments
///
/// * `args` - A slice of string-like values to hash.
///
/// # Returns
///
/// A hexadecimal string containing the MD5 digest.
pub fn generate_md5_key<S: AsRef<str>>(args: &[S]) -> String {
    let joined = args
        .iter()
        .map(|s| s.as_ref())
        .collect::<Vec<_>>()
        .join(" ");

    format!("{:x}", md5::compute(joined))
}