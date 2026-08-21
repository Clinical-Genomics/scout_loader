mod config {
    include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/config.rs"));
}

mod loader {
    include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/loader.rs"));
}

use loader::Loader;
use mongodb::Client;
use mongodb::bson::{Document, doc};
use std::collections::HashSet;

const TEST_CONFIG: &str = "tests/fixtures/test_config.toml";

/// Returns a MongoDB client connected to the test MongoDB instance.
async fn test_client() -> Option<Client> {
    let uri = std::env::var("MONGODB_URI").ok()?;

    Some(
        Client::with_uri_str(&uri)
            .await
            .expect("failed to connect to MongoDB"),
    )
}

/// Inserts test genes for the HGNC mapping tests.
async fn insert_test_genes(collection: &mongodb::Collection<Document>) {
    collection
        .insert_many([
            doc! {
                "hgnc_id": 25662,
                "hgnc_symbol": "AAGAB",
                "build": "37",
            },
            doc! {
                "hgnc_id": 1001,
                "hgnc_symbol": "GENE1",
                "build": "37",
            },
            doc! {
                "hgnc_id": 1002,
                "hgnc_symbol": "GENE2",
                "build": "38",
            },
        ])
        .await
        .expect("failed to insert test genes");
}

/// Inserts test gene panels for the gene-to-panel mapping tests.
async fn insert_test_gene_panels(collection: &mongodb::Collection<Document>) {
    collection
        .insert_many([
            doc! {
                "panel_name": "panel1",
                "genes": [
                    {
                        "hgnc_id": 25662,
                        "symbol": "AAGAB",
                    },
                    {
                        "hgnc_id": 1001,
                        "symbol": "GENE1",
                    },
                ],
            },
            doc! {
                "panel_name": "panel2",
                "genes": [
                    {
                        "hgnc_id": 25662,
                        "symbol": "AAGAB",
                    },
                    {
                        "hgnc_id": 1002,
                        "symbol": "GENE2",
                    },
                ],
            },
        ])
        .await
        .expect("failed to insert test gene panels");
}

/// Tests that HGNC IDs are mapped to genes for the requested genome build.
#[tokio::test]
async fn hgncid_to_gene() {
    let Some(client) = test_client().await else {
        eprintln!("Skipping MongoDB Loader test: MONGODB_URI is not set");
        return;
    };

    let db = client.database("scout_loader_test");
    let collection = db.collection::<Document>("hgnc_gene");

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean hgnc_gene collection");

    insert_test_genes(&collection).await;

    let loader = Loader::new(TEST_CONFIG)
        .await
        .expect("failed to create Loader");

    let genes = loader
        .hgncid_to_gene("37")
        .await
        .expect("failed to build HGNC mapping");

    assert_eq!(genes.len(), 2);

    assert_eq!(
        genes[&25662]
            .get_str("hgnc_symbol")
            .expect("hgnc_symbol should exist"),
        "AAGAB"
    );

    assert_eq!(
        genes[&1001]
            .get_str("hgnc_symbol")
            .expect("hgnc_symbol should exist"),
        "GENE1"
    );

    assert!(!genes.contains_key(&1002));

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean hgnc_gene collection");
}

/// Tests that genes are mapped to all gene panels containing them.
#[tokio::test]
async fn gene_to_panels() {
    let Some(client) = test_client().await else {
        eprintln!("Skipping MongoDB Loader test: MONGODB_URI is not set");
        return;
    };

    let db = client.database("scout_loader_test");
    let collection = db.collection::<Document>("gene_panel");

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean gene_panel collection");

    insert_test_gene_panels(&collection).await;

    let loader = Loader::new(TEST_CONFIG)
        .await
        .expect("failed to create Loader");

    let panel_ids = vec!["panel1".to_string(), "panel2".to_string()];

    let gene_to_panels = loader
        .gene_to_panels(&panel_ids)
        .await
        .expect("failed to build gene-to-panel mapping");

    assert_eq!(
        gene_to_panels[&25662],
        HashSet::from(["panel1".to_string(), "panel2".to_string(),])
    );

    assert_eq!(gene_to_panels[&1001], HashSet::from(["panel1".to_string()]));

    assert_eq!(gene_to_panels[&1002], HashSet::from(["panel2".to_string()]));

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean gene_panel collection");
}

/// Tests that institute_exists returns true for an existing institute
/// and false for an institute that does not exist.
#[tokio::test]
async fn institute_exists() {
    let Some(client) = test_client().await else {
        eprintln!("Skipping MongoDB Loader test: MONGODB_URI is not set");
        return;
    };

    let db = client.database("scout_loader_test");
    let collection = db.collection::<Document>("institute");

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean institute collection");

    collection
        .insert_one(doc! {
            "_id": "test_institute",
        })
        .await
        .expect("failed to insert test institute");

    let loader = Loader::new(TEST_CONFIG)
        .await
        .expect("failed to create Loader");

    assert!(
        loader
            .institute_exists("test_institute")
            .await
            .expect("failed to check existing institute")
    );

    assert!(
        !loader
            .institute_exists("does_not_exist")
            .await
            .expect("failed to check nonexistent institute")
    );

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean institute collection");
}
