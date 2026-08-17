use serde::Deserialize;
use std::fs;

/// MongoDB configuration.
#[derive(Debug, Deserialize)]
pub struct Config {
    pub mongo_uri: String,
    pub mongo_dbname: String,
}

impl Config {
    /// Reads the configuration from a TOML file.
    pub fn from_file(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let contents = fs::read_to_string(path)?;
        let config = toml::from_str(&contents)?;

        Ok(config)
    }
}
