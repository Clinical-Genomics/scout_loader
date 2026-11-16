use std::fs::File;

use tempfile::tempdir;

use scout_loader::models::case::CustomImage;
use scout_loader::utils::custom_images::set_custom_images;

#[test]
fn test_set_custom_images_without_wildcards() {
    let images = vec![
        CustomImage {
            title: "A png image".to_string(),
            str_repid: Some("AFF2".to_string()),
            str_trid: None,
            description: "A very good description".to_string(),
            width: Some(500),
            height: Some(100),
            format: None,
            path: "image.png".to_string(),
        },
        CustomImage {
            title: "An SVG image".to_string(),
            str_repid: Some("DMD".to_string()),
            str_trid: None,
            description: "Another very good description".to_string(),
            width: None,
            height: None,
            format: None,
            path: "image.svg".to_string(),
        },
    ];

    let result = set_custom_images(&images);

    assert_eq!(result.len(), 2);
    assert_eq!(result[0].str_repid.as_deref(), Some("AFF2"));
    assert_eq!(result[0].title, "A png image");
    assert_eq!(result[1].str_repid.as_deref(), Some("DMD"));
}

#[test]
fn test_set_custom_images_replaces_repid_when_file_exists() {
    let temp_dir = tempdir().unwrap();
    let image_path = temp_dir.path().join("640x480_AFF2.svg");

    File::create(&image_path).unwrap();

    let template_path = temp_dir.path().join("640x480_{REPID}.svg");

    let images = vec![CustomImage {
        title: "Image for replicon {REPID}".to_string(),
        str_repid: Some("{REPID}".to_string()),
        str_trid: None,
        description: "Usage of {REPID}".to_string(),
        width: None,
        height: None,
        format: None,
        path: template_path.to_string_lossy().to_string(),
    }];

    let result = set_custom_images(&images);

    assert_eq!(result.len(), 1);
    assert_eq!(result[0].str_repid.as_deref(), Some("AFF2"));
    assert_eq!(result[0].title, "Image for replicon AFF2");
    assert_eq!(result[0].description, "Usage of AFF2");
    assert_eq!(result[0].path, image_path.to_string_lossy());
}

#[test]
fn test_set_custom_images_replaces_trid_when_file_exists() {
    let temp_dir = tempdir().unwrap();
    let image_path = temp_dir.path().join("640x480_TR123.svg");

    File::create(&image_path).unwrap();

    let template_path = temp_dir.path().join("640x480_{TRID}.svg");

    let images = vec![CustomImage {
        title: "Image for repeat {TRID}".to_string(),
        str_repid: None,
        str_trid: Some("{TRID}".to_string()),
        description: "Usage of {TRID}".to_string(),
        width: None,
        height: None,
        format: None,
        path: template_path.to_string_lossy().to_string(),
    }];

    let result = set_custom_images(&images);

    assert_eq!(result.len(), 1);
    assert_eq!(result[0].str_repid, None);
    assert_eq!(result[0].str_trid.as_deref(), Some("TR123"));
    assert_eq!(result[0].title, "Image for repeat TR123");
    assert_eq!(result[0].description, "Usage of TR123");
    assert_eq!(result[0].path, image_path.to_string_lossy());
}

#[test]
fn test_set_custom_images_does_not_add_missing_repid_file() {
    let temp_dir = tempdir().unwrap();

    let images = vec![CustomImage {
        title: "Image for replicon {REPID}".to_string(),
        str_repid: Some("{REPID}".to_string()),
        str_trid: None,
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

    let result = set_custom_images(&images);

    assert!(result.is_empty());
}

#[test]
fn test_set_custom_images_does_not_add_missing_trid_file() {
    let temp_dir = tempdir().unwrap();

    let images = vec![CustomImage {
        title: "Image for repeat {TRID}".to_string(),
        str_repid: None,
        str_trid: Some("{TRID}".to_string()),
        description: "Usage of substitution".to_string(),
        width: None,
        height: None,
        format: None,
        path: temp_dir
            .path()
            .join("640x480_{TRID}.svg")
            .to_string_lossy()
            .to_string(),
    }];

    let result = set_custom_images(&images);

    assert!(result.is_empty());
}
