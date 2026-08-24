use mongodb::{Client, bson::doc};
use std::path::PathBuf;

const MONGO_URI: &str = "mongodb://127.0.0.1:27017";

pub fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

pub struct TestDatabase {
    client: Client,
    name: String,
}

impl TestDatabase {
    pub async fn new(name: &str) -> Self {
        let client = Client::with_uri_str(MONGO_URI)
            .await
            .expect("Failed to connect to MongoDB");

        let name = format!("scout_loader_test_{name}");

        client
            .database(&name)
            .drop()
            .await
            .expect("Failed to clean test database");

        Self { client, name }
    }

    pub async fn count_variants(&self) -> u64 {
        self.client
            .database(&self.name)
            .collection::<mongodb::bson::Document>("variant")
            .count_documents(doc! {})
            .await
            .expect("Failed to count variants")
    }

    pub fn config_path(&self) -> PathBuf {
        let config = format!(
            r#"
mongo_uri = "{MONGO_URI}"
mongo_dbname = "{}"
"#,
            self.name
        );

        let path = std::env::temp_dir().join(format!("{}.toml", self.name));

        std::fs::write(&path, config).expect("Failed to write test config");

        path
    }

    pub async fn cleanup(self) {
        self.client
            .database(&self.name)
            .drop()
            .await
            .expect("Failed to drop test database");
    }
}
