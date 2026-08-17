use mongodb::Client;
use mongodb::bson::{Document, doc};
use std::time::{SystemTime, UNIX_EPOCH};

mod models {
    pub mod variant {
        #[derive(Debug)]
        pub struct VariantIds {
            pub simple_id: String,
            pub variant_id: String,
            pub display_name: String,
            pub document_id: String,
        }
    }
}

mod utils {
    pub mod hash {
        pub fn generate_md5_key(parts: &[&str]) -> String {
            format!("test-md5:{}", parts.join("|"))
        }
    }
}

mod parse {
    pub mod ids {
        include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/parse/ids.rs"));
    }
}

use parse::ids::parse_ids;

fn unique_collection_name() -> String {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock should be after UNIX_EPOCH")
        .as_nanos();

    format!("scout_loader_smoke_{nanos}")
}

#[tokio::test]
async fn mongodb_smoke_roundtrip() {
    let Some(uri) = std::env::var("MONGODB_URI").ok() else {
        eprintln!("Skipping MongoDB smoke test: MONGODB_URI is not set");
        return;
    };

    let client = Client::with_uri_str(&uri)
        .await
        .expect("failed to connect to MongoDB");

    let db = client.database("scout_loader_test");
    let collection_name = unique_collection_name();
    let collection = db.collection::<Document>(&collection_name);

    let expected_id = "smoke_doc";
    let expected_value = 42;

    collection
        .insert_one(doc! {
            "_id": expected_id,
            "value": expected_value
        })
        .await
        .expect("failed to insert smoke document");

    let fetched = collection
        .find_one(doc! { "_id": expected_id })
        .await
        .expect("failed to query smoke document")
        .expect("expected smoke document to exist");

    assert_eq!(fetched.get_i32("value").ok(), Some(expected_value));

    collection
        .drop()
        .await
        .expect("failed to drop smoke test collection");
}

#[tokio::test]
async fn mongodb_parsed_variant_roundtrip() {
    let Some(uri) = std::env::var("MONGODB_URI").ok() else {
        eprintln!("Skipping MongoDB parsed-variant test: MONGODB_URI is not set");
        return;
    };

    let client = Client::with_uri_str(&uri)
        .await
        .expect("failed to connect to MongoDB");

    let db = client.database("scout_loader_test");
    let collection_name = unique_collection_name();
    let collection = db.collection::<Document>(&collection_name);

    let chrom = "1";
    let pos = 123_456_u64;
    let reference = "A";
    let alternative = "T";
    let case_id = "case_123";
    let variant_type = "clinical";

    let ids = parse_ids(chrom, &pos, reference, alternative, case_id, variant_type);

    let variant = doc! {
        "simple_id": ids.simple_id,
        "variant_id": ids.variant_id,
        "display_name": ids.display_name,
        "document_id": ids.document_id,
        "case_id": case_id,
        "chromosome": chrom,
        "position": pos as i64,
        "reference": reference,
        "alternative": alternative,
        "type": variant_type,
        "category": "snv",
        "length": 1_i64,
    };

    collection
        .insert_one(variant.clone())
        .await
        .expect("failed to insert parsed variant document");

    let fetched = collection
        .find_one(doc! { "document_id": &variant["document_id"] })
        .await
        .expect("failed to query parsed variant document")
        .expect("expected parsed variant document to exist");

    assert_eq!(
        fetched.get_str("document_id").ok(),
        variant.get_str("document_id").ok()
    );
    assert_eq!(
        fetched.get_str("simple_id").ok(),
        variant.get_str("simple_id").ok()
    );
    assert_eq!(
        fetched.get_str("display_name").ok(),
        variant.get_str("display_name").ok()
    );
    assert_eq!(fetched.get_str("case_id").ok(), Some(case_id));
    assert_eq!(fetched.get_str("chromosome").ok(), Some(chrom));
    assert_eq!(fetched.get_i64("position").ok(), Some(pos as i64));

    collection
        .drop()
        .await
        .expect("failed to drop parsed variant test collection");
}
