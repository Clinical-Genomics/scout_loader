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

#[tokio::test]
async fn hgncid_to_gene() {
    let Some(uri) = std::env::var("MONGODB_URI").ok() else {
        eprintln!("Skipping MongoDB Loader test: MONGODB_URI is not set");
        return;
    };

    let client = Client::with_uri_str(&uri)
        .await
        .expect("failed to connect to MongoDB");

    let config = config::Config::from_file("config.toml").expect("failed to load config");

    let db = client.database(&config.mongo_dbname);
    let collection = db.collection::<Document>("hgnc_gene");

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean hgnc_gene collection");

    insert_test_genes(&collection).await;

    let loader = Loader::new("config.toml")
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

#[tokio::test]
async fn gene_to_panels() {
    let Some(uri) = std::env::var("MONGODB_URI").ok() else {
        eprintln!("Skipping MongoDB Loader test: MONGODB_URI is not set");
        return;
    };

    let client = Client::with_uri_str(&uri)
        .await
        .expect("failed to connect to MongoDB");

    let config = config::Config::from_file("config.toml").expect("failed to load config");

    let db = client.database(&config.mongo_dbname);
    let collection = db.collection::<Document>("gene_panel");

    collection
        .delete_many(doc! {})
        .await
        .expect("failed to clean gene_panel collection");

    insert_test_gene_panels(&collection).await;

    let loader = Loader::new("config.toml")
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
