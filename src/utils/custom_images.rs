use glob::glob;

use crate::models::case::CustomImage;

/// Expands custom STR images using REPID and TRID wildcards.
///
/// Images using `{REPID}` or `{TRID}` are expanded for each matching file.
/// The wildcard is replaced in the path, title, and description.
///
/// Images without wildcards are included unchanged.
pub fn set_custom_images(images: &[CustomImage]) -> Vec<CustomImage> {
    let mut custom_images = Vec::new();

    for image in images {
        let (identifier, is_repid) = match (image.str_repid.as_deref(), image.str_trid.as_deref()) {
            (Some("{REPID}"), _) => ("{REPID}", true),
            (_, Some("{TRID}")) => ("{TRID}", false),
            _ => {
                custom_images.push(image.clone());
                continue;
            }
        };

        let Some((prefix, suffix)) = image.path.split_once(identifier) else {
            continue;
        };

        let pattern = format!("{prefix}*{suffix}");

        let Ok(paths) = glob(&pattern) else {
            continue;
        };

        for path in paths.flatten() {
            let path = path.to_string_lossy().to_string();

            let Some(matched_id) = path
                .strip_prefix(prefix)
                .and_then(|value| value.strip_suffix(suffix))
                .map(str::to_string)
            else {
                continue;
            };

            let mut new_image = image.clone();

            if is_repid {
                new_image.str_repid = Some(matched_id.clone());
            } else {
                new_image.str_trid = Some(matched_id.clone());
            }

            new_image.path = path;
            new_image.title = image.title.replace(identifier, &matched_id);
            new_image.description = image.description.replace(identifier, &matched_id);

            custom_images.push(new_image);
        }
    }

    custom_images
}
