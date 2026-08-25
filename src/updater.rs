use futures::TryStreamExt;
use mongodb::{
    bson::{Document, doc},
    options::FindOptions,
};

use crate::{
    loader::Loader,
    models::variant::{VariantCategory, VariantType},
};

const UPDATE_BATCH_SIZE: usize = 5000;

/// Updates the manual rank for all variants in a case.
///
/// Variants are sorted by `rank_score` in descending order and assigned
/// a `variant_rank` starting at 1. Updates are sent to MongoDB in batches
/// of `UPDATE_BATCH_SIZE` to reduce the number of database round trips.
///
/// Only variants matching the supplied case, variant type, and category
/// are updated.
///
/// # Arguments
///
/// * `loader` - Loader used to access the MongoDB database.
/// * `case_id` - ID of the case whose variants should be ranked.
/// * `variant_type` - Variant type of the variants to update.
/// * `category` - Variant category of the variants to update.
///
/// # Errors
///
/// Returns an error if the variants cannot be queried or their ranks
/// cannot be updated.
pub async fn update_variant_rank(
    loader: &Loader,
    case_id: &str,
    variant_type: VariantType,
    category: VariantCategory,
) -> Result<(), Box<dyn std::error::Error>> {
    let collection = loader.variant_collection();

    let filter = doc! {
        "case_id": case_id,
        "variant_type": variant_type.to_string(),
        "category": category.to_string(),
    };

    let options = FindOptions::builder()
        .sort(doc! { "rank_score": -1 })
        .build();

    let mut cursor = collection.find(filter).with_options(options).await?;

    let mut updates = Vec::with_capacity(UPDATE_BATCH_SIZE);
    let mut variant_rank = 1i32;

    while let Some(variant) = cursor.try_next().await? {
        let variant_id = variant.get("_id").ok_or("Variant is missing _id")?.clone();

        updates.push(doc! {
            "q": {
                "_id": variant_id,
            },
            "u": {
                "$set": {
                    "variant_rank": variant_rank,
                },
            },
            "multi": false,
        });

        variant_rank += 1;

        if updates.len() >= UPDATE_BATCH_SIZE {
            execute_updates(loader, updates).await?;
            updates = Vec::with_capacity(UPDATE_BATCH_SIZE);
        }
    }

    if !updates.is_empty() {
        execute_updates(loader, updates).await?;
    }

    Ok(())
}

/// Executes a batch of variant rank updates using MongoDB's `update`
/// database command.
///
/// This works with MongoDB 6, 7, and 8 and avoids the MongoDB 8.0-only
/// client-level `bulkWrite` command.
async fn execute_updates(
    loader: &Loader,
    updates: Vec<Document>,
) -> Result<(), Box<dyn std::error::Error>> {
    let command = doc! {
        "update": "variant",
        "updates": updates,
        "ordered": false,
    };

    loader.database().run_command(command).await?;

    Ok(())
}
