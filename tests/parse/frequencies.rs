use mongodb::bson::{Bson, Document, doc};
use rust_htslib::bcf::Record;

use scout_loader::parse::frequencies::{
    GNOMAD_INFO_KEYS, GNOMAD_INFO_MAX_KEYS, add_frequencies, parse_frequencies,
    parse_mei_frequencies, parse_sv_frequencies,
};

use crate::vcf_record::record;

/// Build a record with the given INFO fields, each declared with `Number=1` and its `Type`.
fn info_record(fields: &[(&str, &str, &str)]) -> Record {
    let meta: Vec<String> = fields
        .iter()
        .map(|(key, kind, _)| {
            format!("##INFO=<ID={key},Number=1,Type={kind},Description=\"{key}\">")
        })
        .collect();
    let meta: Vec<&str> = meta.iter().map(String::as_str).collect();

    let info = fields
        .iter()
        .map(|(key, _, value)| format!("{key}={value}"))
        .collect::<Vec<_>>()
        .join(";");

    record(&meta, &format!("1\t100\t.\tA\tT\t.\t.\t{info}"))
}

/// Assert that `doc[key]` is a double close to `expected`.
fn assert_f64(doc: &Document, key: &str, expected: f64) {
    let value = doc
        .get_f64(key)
        .unwrap_or_else(|_| panic!("{key} missing or not a double in {doc}"));
    assert!(
        (value - expected).abs() < 1e-6,
        "{key}: expected {expected}, got {value}"
    );
}

#[test]
fn test_missing_frequency_uses_transcript_fallback() {
    // GIVEN a missing frequency in either a Float or String INFO field
    for field_type in ["Float", "String"] {
        for value in [".", "NaN", "0", "-1"] {
            let header = format!(
                "##INFO=<ID=GNOMADAF,Number=1,Type={field_type},Description=\"Frequency\">"
            );
            let record = record(
                &["##contig=<ID=1,length=249250621>", &header],
                &format!("1\t100\t.\tA\tG\t.\tPASS\tGNOMADAF={value}"),
            );

            // WHEN parsing with a transcript frequency
            let frequencies = parse_frequencies(&record, &[doc! { "gnomad_maf": 0.02 }]);

            // THEN placeholders do not block the transcript fallback
            assert_eq!(
                frequencies,
                doc! { "gnomad": 0.02 },
                "{field_type}: {value}"
            );
            assert!(parse_frequencies(&record, &[]).is_empty());
        }
    }
}

#[test]
fn test_missing_frequency_uses_next_info_key() {
    // GIVEN a missing primary frequency and a valid alternative INFO key
    let record = record(
        &[
            "##contig=<ID=1,length=249250621>",
            "##INFO=<ID=GNOMADAF,Number=1,Type=Float,Description=\"Frequency\">",
            "##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description=\"Alternative frequency\">",
        ],
        "1\t100\t.\tA\tG\t.\tPASS\tGNOMADAF=.;GNOMAD_AF=0.25",
    );

    // WHEN parsing frequencies THEN the next INFO key is used
    assert_eq!(parse_frequencies(&record, &[]), doc! { "gnomad": 0.25 });
}

#[test]
fn test_parse_gnomad_popmax() {
    for key in GNOMAD_INFO_MAX_KEYS {
        // GIVEN a record with a gnomAD popmax frequency under one of the known keys
        let record = info_record(&[(key, "Float", "0.05")]);

        // WHEN parsing frequencies
        let frequencies = parse_frequencies(&record, &[]);

        // THEN the frequency is stored as gnomad_max
        assert_f64(&frequencies, "gnomad_max", 0.05);
    }
}

#[test]
fn test_parse_gnomad() {
    for key in GNOMAD_INFO_KEYS {
        // GIVEN a record with a gnomAD frequency under one of the known keys,
        // including the SV (gnomad_svAF) and echtvar (gnomad_af) variants
        let record = info_record(&[(key, "Float", "0.117036")]);

        // WHEN parsing frequencies
        let frequencies = parse_frequencies(&record, &[]);

        // THEN the frequency is stored as gnomad
        assert_f64(&frequencies, "gnomad", 0.117036);
    }
}

#[test]
fn test_parse_frequencies() {
    // GIVEN a record with 1000G, ExAC, SweGen and CoLoRSdb frequencies
    let record = info_record(&[
        ("1000GAF", "Float", "0.01"),
        ("EXACAF", "Float", "0.001"),
        ("SWEGENAF", "Float", "0.02"),
        ("colorsdb_af", "Float", "0.4982"),
    ]);

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &[]);

    // THEN each frequency is stored under its own key
    assert_f64(&frequencies, "thousand_g", 0.01);
    assert_f64(&frequencies, "exac", 0.001);
    assert_f64(&frequencies, "swegen", 0.02);
    assert_f64(&frequencies, "colorsdb_af", 0.4982);
}

#[test]
fn test_parse_frequencies_string_type() {
    // GIVEN a frequency declared as a String INFO field
    let record = info_record(&[("GNOMADAF", "String", "0.003")]);

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &[]);

    // THEN the string is parsed as a frequency
    assert_f64(&frequencies, "gnomad", 0.003);
}

#[test]
fn test_parse_frequencies_placeholders() {
    // GIVEN frequencies with the placeholder values 0 and -1
    let record = info_record(&[
        ("GNOMADAF", "Float", "0"),
        ("EXACAF", "Float", "-1"),
        ("SWEGENAF", "String", "-1"),
    ]);

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &[]);

    // THEN none of them are stored
    assert!(
        frequencies.is_empty(),
        "unexpected frequencies: {frequencies}"
    );
}

#[test]
fn test_parse_frequencies_first_key_wins() {
    // GIVEN two gnomAD frequency keys
    let record = info_record(&[("GNOMADAF", "Float", "0.1"), ("GNOMAD_AF", "Float", "0.2")]);

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &[]);

    // THEN the first key in the key list is used
    assert_f64(&frequencies, "gnomad", 0.1);
}

#[test]
fn test_parse_frequencies_gnomad_counts() {
    // GIVEN gnomAD genome and exome allele counts and homozygote counts
    let record = info_record(&[
        ("GNOMADAC", "Integer", "12"),
        ("GNOMAD_nhomalt", "Integer", "3"),
        ("GNOMAD_nhomalt_XX", "Integer", "2"),
        ("GNOMAD_nhomalt_XY", "Integer", "1"),
        ("GNOMAD_WES_AC", "Integer", "40"),
        ("GNOMAD_WES_nhomalt", "Integer", "5"),
    ]);

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &[]);

    // THEN the counts are stored as 32-bit integers
    assert_eq!(frequencies.get("gnomad_obs"), Some(&Bson::Int32(12)));
    assert_eq!(frequencies.get("gnomad_nhomalt"), Some(&Bson::Int32(3)));
    assert_eq!(frequencies.get("gnomad_nhomalt_XX"), Some(&Bson::Int32(2)));
    assert_eq!(frequencies.get("gnomad_nhomalt_XY"), Some(&Bson::Int32(1)));
    assert_eq!(frequencies.get("gnomad_wes_obs"), Some(&Bson::Int32(40)));
    assert_eq!(frequencies.get("gnomad_wes_nhomalt"), Some(&Bson::Int32(5)));
}

#[test]
fn test_parse_frequencies_from_transcripts() {
    // GIVEN a record without INFO frequencies and VEP transcripts with frequencies
    let record = info_record(&[("DP", "Integer", "30")]);
    let transcripts = vec![
        doc! { "gnomad_maf": 0.02, "exac_maf": Bson::Null },
        doc! { "thousand_g_maf": 0.03, "gnomad_max": 0.04 },
    ];

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &transcripts);

    // THEN the transcript frequencies are used, skipping null values
    assert_f64(&frequencies, "gnomad", 0.02);
    assert_f64(&frequencies, "thousand_g", 0.03);
    assert_f64(&frequencies, "gnomad_max", 0.04);
    assert!(!frequencies.contains_key("exac"));
}

#[test]
fn test_parse_frequencies_info_takes_precedence_over_transcripts() {
    // GIVEN a record with an INFO frequency and a transcript with another frequency
    let record = info_record(&[("GNOMADAF", "Float", "0.1")]);
    let transcripts = vec![doc! { "thousand_g_maf": 0.03 }];

    // WHEN parsing frequencies
    let frequencies = parse_frequencies(&record, &transcripts);

    // THEN only the INFO frequency is used
    assert_f64(&frequencies, "gnomad", 0.1);
    assert!(!frequencies.contains_key("thousand_g"));
}

#[test]
fn test_parse_sv_frequencies_clingen_benign() {
    // GIVEN ClinGen CGH benign frequency and occurrence
    let record = info_record(&[
        ("clingen_cgh_benignAF", "Float", "0.01"),
        ("clingen_cgh_benign", "Integer", "3"),
    ]);

    // WHEN parsing SV frequencies
    let frequencies = parse_sv_frequencies(&record);

    // THEN the first key, the frequency, is used
    assert_f64(&frequencies, "clingen_cgh_benign", 0.01);
}

#[test]
fn test_parse_sv_frequencies_clingen_pathogenic() {
    // GIVEN ClinGen CGH pathogenic frequency, occurrence and OCC
    let record = info_record(&[
        ("clingen_cgh_pathogenicAF", "Float", "0.01"),
        ("clingen_cgh_pathogenic", "Integer", "3"),
        ("clingen_cgh_pathogenicOCC", "Integer", "1"),
    ]);

    // WHEN parsing SV frequencies
    let frequencies = parse_sv_frequencies(&record);

    // THEN the first key, the frequency, is used
    assert_f64(&frequencies, "clingen_cgh_pathogenic", 0.01);
}

#[test]
fn test_parse_sv_frequencies_clingen_mip() {
    // GIVEN Clinical Genomics MIP frequency and occurrence
    let record = info_record(&[
        ("clinical_genomics_mipAF", "Float", "0.01"),
        ("clinical_genomics_mipOCC", "Integer", "2"),
    ]);

    // WHEN parsing SV frequencies
    let frequencies = parse_sv_frequencies(&record);

    // THEN the frequency is stored as clingen_mip
    assert_f64(&frequencies, "clingen_mip", 0.01);
}

#[test]
fn test_parse_sv_frequencies_occurrence() {
    // GIVEN a ClinGen NGI occurrence count only, and a zero decipher count
    let record = info_record(&[
        ("clingen_ngi", "Integer", "1"),
        ("decipher", "Integer", "0"),
    ]);

    // WHEN parsing SV frequencies
    let frequencies = parse_sv_frequencies(&record);

    // THEN the occurrence is stored as an integer
    assert_eq!(frequencies.get("clingen_ngi"), Some(&Bson::Int32(1)));

    // THEN zero occurrences are ignored
    assert!(!frequencies.contains_key("decipher"));
}

#[test]
fn test_parse_mei_frequencies() {
    // GIVEN SweGen ALU and L1 MEI frequencies, and a SVA occurrence count
    let record = info_record(&[
        ("swegen_alu_FRQ", "Float", "0.2"),
        ("swegen_alu_OCC", "Integer", "50"),
        ("swegen_l1_FRQ", "Float", "0.4"),
        ("swegen_sva_OCC", "Integer", "3"),
    ]);

    // WHEN parsing MEI frequencies
    let frequencies = parse_mei_frequencies(&record);

    // THEN the frequency is used where present, otherwise the occurrence
    assert_f64(&frequencies, "swegen_alu", 0.2);
    assert_f64(&frequencies, "swegen_l1", 0.4);
    assert_eq!(frequencies.get("swegen_sva"), Some(&Bson::Int32(3)));
    assert!(!frequencies.contains_key("swegen_herv"));

    // THEN the maximum value is stored as swegen_mei_max
    assert_f64(&frequencies, "swegen_mei_max", 3.0);
}

#[test]
fn test_parse_mei_frequencies_missing() {
    // GIVEN a record without MEI frequencies
    let record = info_record(&[("swegen_alu_FRQ", "Float", ".")]);

    // WHEN parsing MEI frequencies
    let frequencies = parse_mei_frequencies(&record);

    // THEN no frequencies, and no maximum, are returned
    assert!(frequencies.is_empty());
}

#[test]
fn test_add_frequencies() {
    // GIVEN parsed frequencies, including SV and MEI frequencies
    let frequencies = doc! {
        "gnomad": 0.1,
        "gnomad_max": 0.2,
        "thousand_g": 0.3,
        "gnomad_nhomalt_XX": 2,
        "clingen_ngi": 4,
        "swegen_alu": 0.5,
    };
    let mut variant = Document::new();

    // WHEN adding them to the variant
    add_frequencies(&mut variant, &frequencies);

    // THEN they are renamed to the Scout variant keys
    assert_eq!(
        variant,
        doc! {
            "gnomad_frequency": 0.1,
            "max_gnomad_frequency": 0.2,
            "gnomad_nhomalt_xx_frequency": 2,
            "thousand_genomes_frequency": 0.3,
            "clingen_ngi": 4,
            "swegen_alu": 0.5,
        }
    );
}
