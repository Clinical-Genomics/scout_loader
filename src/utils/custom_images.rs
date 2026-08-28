use std::path::Path;

use crate::models::case::CustomImage;

/// Gets the custom images matching an STR variant.
///
/// Images with a `str_repid` matching the variant are included directly.
/// Images using the `{REPID}` placeholder are resolved using the variant's
/// `str_repid` and included only if the resulting image file exists on disk.
///
/// # Arguments
///
/// * `images` - Configured custom STR variant images.
/// * `str_repid` - The STR replicon ID of the variant.
///
/// # Returns
///
/// A vector containing the custom images applicable to the STR variant.
pub fn set_custom_images(images: &[CustomImage], str_repid: &str) -> Vec<CustomImage> {
    let mut custom_images = Vec::new();

    for image in images {
        if image.str_repid == str_repid {
            custom_images.push(image.clone());
        } else if image.str_repid == "{REPID}" {
            let path = image.path.replace("{REPID}", str_repid);

            if !Path::new(&path).exists() {
                continue;
            }

            let mut new_image = image.clone();
            new_image.str_repid = str_repid.to_string();
            new_image.path = path;
            new_image.title = image.title.replace("{REPID}", str_repid);
            new_image.description = image.description.replace("{REPID}", str_repid);

            custom_images.push(new_image);
        }
    }

    custom_images
}
