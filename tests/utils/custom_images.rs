use std::fs::File;

use tempfile::tempdir;

use scout_loader::models::case::CustomImage;
use scout_loader::utils::custom_images::set_custom_images;

#[test]
fn test_set_custom_images_matching_repid() {
    let images = vec![
        CustomImage {
            title: "A png image".to_string(),
            str_repid: "AFF2".to_string(),
            description: "A very good description".to_string(),
            width: Some(500),
            height: Some(100),
            format: None,
            path: "image.png".to_string(),
        },
        CustomImage {
            title: "An SVG image".to_string(),
            str_repid: "DMD".to_string(),
            description: "Another very good description".to_string(),
            width: None,
            height: None,
            format: None,
            path: "image.svg".to_string(),
        },
    ];

    let result = set_custom_images(&images, "AFF2");

    assert_eq!(result.len(), 1);
    assert_eq!(result[0].str_repid, "AFF2");
    assert_eq!(result[0].title, "A png image");
}

#[test]
fn test_set_custom_images_replaces_repid_when_file_exists() {
    let temp_dir = tempdir().unwrap();
    let image_path = temp_dir.path().join("640x480_AFF2.svg");

    File::create(&image_path).unwrap();

    let template_path = temp_dir.path().join("640x480_{REPID}.svg");

    let images = vec![CustomImage {
        title: "Image for replicon {REPID}".to_string(),
        str_repid: "{REPID}".to_string(),
        description: "Usage of substitution".to_string(),
        width: None,
        height: None,
        format: None,
        path: template_path.to_string_lossy().to_string(),
    }];

    let result = set_custom_images(&images, "AFF2");

    assert_eq!(result.len(), 1);
    assert_eq!(result[0].str_repid, "AFF2");
    assert_eq!(result[0].title, "Image for replicon AFF2");
    assert_eq!(result[0].description, "Usage of substitution");
    assert_eq!(result[0].path, image_path.to_string_lossy());
}

#[test]
fn test_set_custom_images_does_not_add_missing_file() {
    let temp_dir = tempdir().unwrap();

    let images = vec![CustomImage {
        title: "Image for replicon {REPID}".to_string(),
        str_repid: "{REPID}".to_string(),
        description: "Usage of substitution".to_string(),
        width: None,
        height: None,
        format: None,
        path: temp_dir
            .path()
            .join("640x480_{REPID}.svg")
            .to_string_lossy()
            .to_string(),
    }];

    let result = set_custom_images(&images, "AFF2");

    assert!(result.is_empty());
}
