use scout_loader::parse::frequencies::{
    parse_frequencies, GNOMAD_INFO_MAX_KEYS,
};

#[test]
fn test_parse_gnomad_popmax() {
    for key in GNOMAD_INFO_MAX_KEYS {
        let info = format!("{key}=0.05");
        let record = create_test_record(&info);

        let frequencies = parse_frequencies(&record, &[]);

        assert_eq!(
            frequencies.get_f64("gnomad_max").unwrap(),
            0.05,
            "Failed to parse INFO field: {key}",
        );
    }
}