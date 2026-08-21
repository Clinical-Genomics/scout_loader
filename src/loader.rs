use futures::TryStreamExt;
use mongodb::Database;
use mongodb::bson::{Document, doc};
use std::collections::{HashMap, HashSet};

use crate::config::Config;

/// Handles loading cases and manages the MongoDB connection.
pub struct Loader {
    db: Database,
}

impl Loader {
    /// Creates a new loader from the given configuration file.
    ///
    /// Reads the MongoDB configuration and establishes a database connection.
    pub async fn new(config_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let config = Config::from_file(config_path)?;

        let client = mongodb::Client::with_uri_str(&config.mongo_uri).await?;
        let db = client.database(&config.mongo_dbname);

        Ok(Self { db })
    }

    /// Build a mapping from HGNC ID to the corresponding gene document.
    ///
    /// Fetches genes from MongoDB for the specified genome build and uses the
    /// HGNC ID as the key.
    pub async fn hgncid_to_gene(
        &self,
        build: &str,
    ) -> Result<HashMap<i32, Document>, Box<dyn std::error::Error>> {
        let collection = self.db.collection::<Document>("hgnc_gene");

        let filter = doc! {
            "build": build
        };

        let mut genes = collection.find(filter).await?;

        let mut hgnc_dict = HashMap::new();

        while let Some(gene) = genes.try_next().await? {
            let hgnc_id = gene.get_i32("hgnc_id")?;
            hgnc_dict.insert(hgnc_id, gene);
        }

        Ok(hgnc_dict)
    }

    /// Check whether an institute exists in the database.
    ///
    /// # Arguments
    ///
    /// * `institute_id` - Identifier of the institute to look up.
    ///
    /// # Returns
    ///
    /// `Ok(true)` if the institute exists, `Ok(false)` if it does not,
    /// or a MongoDB error if the database query fails.
    pub async fn institute_exists(&self, institute_id: &str) -> mongodb::error::Result<bool> {
        let collection = self.db.collection::<Document>("institute");

        Ok(collection
            .find_one(doc! { "_id": institute_id })
            .await?
            .is_some())
    }

    /// Build a mapping of HGNC IDs to the gene panels containing each gene.
    ///
    /// Fetches the requested gene panels from MongoDB and collects the panel
    /// names for each HGNC ID.
    ///
    /// # Arguments
    ///
    /// * `panel_ids` - IDs of the gene panels to retrieve.
    ///
    /// # Returns
    ///
    /// A mapping from HGNC ID to the set of panel names containing that gene.
    pub async fn gene_to_panels(
        &self,
        panel_ids: &[String],
    ) -> Result<HashMap<i32, HashSet<String>>, Box<dyn std::error::Error>> {
        let collection = self.db.collection::<Document>("gene_panel");

        let filter = doc! {
            "panel_name": {
                "$in": panel_ids
            }
        };

        let mut panels = collection.find(filter).await?;

        let mut gene_dict: HashMap<i32, HashSet<String>> = HashMap::new();

        while let Some(panel) = panels.try_next().await? {
            let panel_name = panel.get_str("panel_name")?;
            let genes = panel.get_array("genes")?;

            for gene in genes {
                let gene = gene.as_document().ok_or("Invalid gene document")?;

                let hgnc_id = gene.get_i32("hgnc_id")?;

                gene_dict
                    .entry(hgnc_id)
                    .or_default()
                    .insert(panel_name.to_string());
            }
        }

        Ok(gene_dict)
    }

    pub async fn count_variants(&self) -> Result<u64, mongodb::error::Error> {
        let collection = self.db.collection::<Document>("variant");

        collection.count_documents(doc! {}).await
    }

    /// Loads a batch of variants into the database.
    ///
    /// Variants are inserted in bulk to reduce the number of database
    /// round trips. If the bulk insertion fails because some variants
    /// already exist, each variant is inserted individually.
    pub async fn load_variant_bulk(
        &self,
        variants: Vec<Document>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if variants.is_empty() {
            return Ok(());
        }

        let collection = self.db.collection::<Document>("variant");

        match collection.insert_many(variants.clone()).await {
            Ok(_) => Ok(()),
            Err(mongodb::error::Error { .. }) => {
                // If we need the same fallback behaviour as the old Python
                // implementation, insert the variants individually here.
                for variant in variants {
                    collection.insert_one(variant).await?;
                }

                Ok(())
            }
        }
    }
}
