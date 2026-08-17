use mongodb::sync::Database;

use crate::config::Config;

/// Handles loading cases and manages the MongoDB connection.
pub struct Loader {
    db: Database,
}

impl Loader {
    /// Creates a new loader from the given configuration file.
    ///
    /// Reads the MongoDB configuration and establishes a database connection.
    pub fn new(config_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let config = Config::from_file(config_path)?;

        let client = mongodb::sync::Client::with_uri_str(&config.mongo_uri)?;
        let db = client.database(&config.mongo_dbname);

        Ok(Self { db })
    }
}
