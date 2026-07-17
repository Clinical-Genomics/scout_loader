use mongodb::bson::{doc, Document};
use rust_htslib::bcf::Record;

use crate::parse::info::parse_info_string;

/// Add mobile element insertion (MEI) annotations to a MongoDB variant document.
///
/// The MEINFO INFO field is expected to have the format:
///
/// `NAME,START,END,POLARITY`
///
/// If the field is present and correctly formatted, a nested `mei` document
/// is added to the variant:
///
/// ```text
/// "mei": {
///     "name": "...",
///     "polarity": "..."
/// }
/// ```
///
/// Missing or malformed MEINFO fields are ignored.
pub fn set_mei_info(record: &Record, variant: &mut Document) {
    let Some(mei_info) = parse_info_string(record, b"MEINFO") else {
        return;
    };

    let fields: Vec<&str> = mei_info.split(',').collect();

    if fields.len() != 4 {
        return;
    }

    variant.insert(
        "mei",
        doc! {
            "name": fields[0],
            "polarity": fields[3],
        },
    );
}
