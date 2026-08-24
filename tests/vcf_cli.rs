use mongodb::{Client, Database, bson::doc};
use std::path::PathBuf;

pub struct TestDatabase {
    db: Database,
}

impl TestDatabase {
    pub async fn new() -> Self {
        let client = Client::with_uri_str("mongodb://localhost:27017")
            .await
            .expect("failed to connect to MongoDB");

        let db = client.database("scout_loader_test");

        db.drop().await.expect("failed to reset test database");

        Self { db }
    }

    pub async fn count_variants(&self) -> u64 {
        self.db
            .collection::<mongodb::bson::Document>("variant")
            .count_documents(doc! {})
            .await
            .expect("failed to count variants")
    }
}

pub fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}
