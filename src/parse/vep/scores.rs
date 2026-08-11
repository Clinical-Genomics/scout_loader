use crate::HashMap;
use mongodb::bson::{Bson, Document};

/// Parse transcript-level CADD Phred score.
pub fn parse_cadd(transcript: &mut Document, entry: &HashMap<String, String>) {
    println!("CADD_PHRED entry: {:?}", entry.get("CADD_PHRED"));
    if let Some(cadd_phred) = entry.get("CADD_PHRED")
        && let Ok(value) = cadd_phred.parse::<f64>()
    {
        transcript.insert("cadd", Bson::Double(value));
    }
}
