use mongodb::bson::{doc, Document};
use rust_htslib::bcf::Record;

use crate::parse::info::parse_info_string_array;

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
    let Some(mei_info) = parse_info_string_array(record, b"MEINFO") else {
        return;
    };

    if mei_info.len() != 4 {
        return;
    }

    variant.insert("mei_name", mei_info[0].clone());
    variant.insert("mei_polarity", mei_info[3].clone());
}
