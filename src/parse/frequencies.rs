use crate::parse::info::parse_info_int;
use mongodb::bson::{Bson, Document};
use rust_htslib::bcf::Record;

pub const CG_KEYS: &[&str] = &["clinical_genomics_mipAF", "clinical_genomics_mipOCC"];

pub const CLINGEN_BENIGN_KEYS: &[&str] = &[
    "clingen_cgh_benignAF",
    "clingen_cgh_benign",
    "clingen_cgh_benignOCC",
];

pub const CLINGEN_NGI_KEYS: &[&str] = &["clingen_ngi", "clingen_ngiAF", "clingen_ngiOCC"];

pub const CLINGEN_PATHOGENIC_KEYS: &[&str] = &[
    "clingen_cgh_pathogenicAF",
    "clingen_cgh_pathogenic",
    "clingen_cgh_pathogenicOCC",
];

pub const DECIPHER_KEYS: &[&str] = &["decipherAF", "decipher"];

pub const EXAC_KEYS: &[&str] = &["EXACAF"];

pub const EXAC_MAX_KEYS: &[&str] = &["EXACMAXAF", "EXACMAX_AF", "ExAC_MAX_AF"];

pub const GNOMAD_INFO_KEYS: &[&str] = &[
    "GNOMADAF",
    "GNOMAD_AF",
    "gnomADg_AF",
    "gnomad_svAF",
    "gnomad_af",
];

pub const GNOMAD_INFO_MAX_KEYS: &[&str] = &[
    "gnomADg_AF_POPMAX",
    "GNOMADAF_popmax",
    "GNOMADAF_POPMAX",
    "GNOMADAF_MAX",
    "gnomad_popmax_af",
    "GNOMADAF_grpmax",
    "MAX_AF",
];

pub const GNOMAD_NHOMALT: &[&str] = &["GNOMAD_nhomalt"];

pub const GNOMAD_NHOMALT_XX: &[&str] = &["GNOMAD_nhomalt_XX"];

pub const GNOMAD_NHOMALT_XY: &[&str] = &["GNOMAD_nhomalt_XY"];

pub const GNOMAD_OBS: &[&str] = &["GNOMADAC"];

pub const GNOMAD_WES: &[&str] = &["GNOMAD_WES_AF"];

pub const GNOMAD_WES_MAX: &[&str] = &["GNOMAD_WES_AF_grpmax"];

pub const GNOMAD_WES_NHOMALT: &[&str] = &["GNOMAD_WES_nhomalt"];

pub const GNOMAD_WES_NHOMALT_XX: &[&str] = &["GNOMAD_WES_nhomalt_XX"];

pub const GNOMAD_WES_NHOMALT_XY: &[&str] = &["GNOMAD_WES_nhomalt_XY"];

pub const GNOMAD_WES_OBS: &[&str] = &["GNOMAD_WES_AC"];

pub const SWEGEN_KEYS: &[&str] = &["SWEGENAF", "swegenAF", "swegen"];

pub const SWEGEN_ALU_KEYS: &[&str] = &["swegen_alu_FRQ", "swegen_alu_OCC"];

pub const SWEGEN_HERV_KEYS: &[&str] = &["swegen_herv_FRQ", "swegen_herv_OCC"];

pub const SWEGEN_L1_KEYS: &[&str] = &["swegen_l1_FRQ", "swegen_l1_OCC"];

pub const SWEGEN_SVA_KEYS: &[&str] = &["swegen_sva_FRQ", "swegen_sva_OCC"];

pub const THOUSAND_GENOMES_KEYS: &[&str] = &["1000GAF"];

pub const THOUSAND_GENOMES_MAX_KEYS: &[&str] = &["1000G_MAX_AF"];

/// Parse a frequency value from a VCF INFO field.
///
/// Returns None if the field is missing, contains a placeholder value
/// (., 0, -1), or cannot be parsed as a frequency.
///
/// Supports both numeric INFO fields and string fields containing numeric
/// frequency values. Returns the first frequency value as f64.
fn parse_frequency(record: &Record, key: &[u8]) -> Option<f64> {
    record
        .info(key)
        .float()
        .ok()
        .flatten()
        .and_then(|values| values.first().copied())
        .filter(|v| *v != 0.0 && *v != -1.0)
        .map(|v| v as f64)
        .or_else(|| {
            record
                .info(key)
                .string()
                .ok()
                .flatten()
                .and_then(|values| {
                    values
                        .first()
                        .and_then(|value| std::str::from_utf8(value).ok())
                        .and_then(|value| value.parse::<f64>().ok())
                })
                .filter(|v| *v != 0.0 && *v != -1.0)
        })
}

/// Update frequency document from VCF INFO fields.
///
/// Searches the provided INFO keys in order. The first valid frequency
/// found is stored under `new_key`.
///
/// # Arguments
///
/// * `frequency` - Frequency BSON document to update.
/// * `record` - VCF record.
/// * `key_list` - INFO field keys to search.
/// * `new_key` - Database key where the frequency is stored.
pub fn update_frequency_from_vcf(
    frequency: &mut Document,
    record: &Record,
    key_list: &[&str],
    new_key: &str,
) {
    for key in key_list {
        if let Some(value) = parse_frequency(record, key.as_bytes()) {
            frequency.insert(new_key, Bson::Double(value));
            break;
        }
    }
}

/// Parse variant frequencies from VCF INFO fields or VEP transcripts.
///
/// Frequencies are first searched in INFO fields. If none are found, the
/// transcript-level VEP annotations are used.
///
/// Returns a BSON document containing available frequencies.
pub fn parse_frequencies(record: &Record, transcripts: &[Document]) -> Document {
    let mut frequencies = Document::new();

    let frequency_fields = [
        (EXAC_KEYS, "exac"),
        (EXAC_MAX_KEYS, "exac_max"),
        (GNOMAD_INFO_KEYS, "gnomad"),
        (SWEGEN_KEYS, "swegen"),
        (GNOMAD_INFO_MAX_KEYS, "gnomad_max"),
        (THOUSAND_GENOMES_KEYS, "thousand_g"),
        (THOUSAND_GENOMES_MAX_KEYS, "thousand_g_max"),
        (&["GNOMAD_MT_AF_HOM"][..], "gnomad_mt_homoplasmic"),
        (&["GNOMAD_MT_AF_HET"][..], "gnomad_mt_heteroplasmic"),
        (&["left_1000GAF"][..], "thousand_g_left"),
        (&["right_1000GAF"][..], "thousand_g_right"),
        (&["colorsdb_af"][..], "colorsdb_af"),
        (GNOMAD_WES, "gnomad_wes"),
        (GNOMAD_WES_MAX, "gnomad_wes_max"),
    ];

    for (keys, name) in frequency_fields {
        update_frequency_from_vcf(&mut frequencies, record, keys, name);
    }

    let integer_fields = [
        (GNOMAD_OBS, "gnomad_obs"),
        (GNOMAD_NHOMALT, "gnomad_nhomalt"),
        (GNOMAD_NHOMALT_XX, "gnomad_nhomalt_XX"),
        (GNOMAD_NHOMALT_XY, "gnomad_nhomalt_XY"),
        (GNOMAD_WES_OBS, "gnomad_wes_obs"),
        (GNOMAD_WES_NHOMALT, "gnomad_wes_nhomalt"),
        (GNOMAD_WES_NHOMALT_XX, "gnomad_wes_nhomalt_XX"),
        (GNOMAD_WES_NHOMALT_XY, "gnomad_wes_nhomalt_XY"),
    ];

    for (keys, name) in integer_fields {
        for key in keys {
            if let Some(value) = parse_info_int(record, key.as_bytes()) {
                frequencies.insert(name, Bson::Int32(value));
                break;
            }
        }
    }

    if frequencies.is_empty() {
        update_frequency_from_transcript(&mut frequencies, transcripts);
    }

    frequencies
}

/// Update frequencies from VEP transcript annotations.
///
/// Searches transcript-level annotations for population frequencies and adds
/// them to the frequency document.
pub fn update_frequency_from_transcript(frequencies: &mut Document, transcripts: &[Document]) {
    for transcript in transcripts {
        let frequency_fields = [
            ("exac_maf", "exac"),
            ("exac_max", "exac_max"),
            ("thousand_g_maf", "thousand_g"),
            ("thousandg_max", "thousand_g_max"),
            ("gnomad_maf", "gnomad"),
            ("gnomad_max", "gnomad_max"),
            ("gnomad_mt_homoplasmic", "gnomad_mt_homoplasmic"),
            ("gnomad_mt_heteroplasmic", "gnomad_mt_heteroplasmic"),
        ];

        for (transcript_key, frequency_key) in frequency_fields {
            if let Some(value) = transcript.get(transcript_key)
                && !matches!(value, Bson::Null)
            {
                frequencies.insert(frequency_key, value.clone());
            }
        }
    }
}

/// Add a frequency annotation to a variant if it exists.
///
/// Copies a value from the parsed frequency document to the target variant
/// using the provided source and target field names.
///
/// # Arguments
///
/// * `variant` - Mutable MongoDB document representing the parsed variant.
/// * `frequencies` - Document containing parsed frequency annotations.
/// * `source_key` - Key used in the frequency document.
/// * `target_key` - Field name to add to the variant document.
fn add_frequency_field(
    variant: &mut Document,
    frequencies: &Document,
    source_key: &str,
    target_key: &str,
) {
    if let Some(value) = frequencies.get(source_key) {
        variant.insert(target_key, value.clone());
    }
}

/// Add frequency annotations to a variant.
///
/// Frequencies are first collected from VCF INFO fields or transcript
/// annotations into an intermediate frequency document. This function maps
/// those internal keys to the field names expected in the Scout MongoDB
/// variant document.
///
/// Missing frequencies are ignored and are not added to the variant.
///
/// # Arguments
///
/// * `variant` - Mutable MongoDB document representing the parsed variant.
/// * `frequencies` - Document containing parsed frequency annotations.
pub fn add_frequencies(variant: &mut Document, frequencies: &Document) {
    let mappings = [
        ("colorsdb_af", "colorsdb_af"),
        ("exac", "exac_frequency"),
        ("exac_max", "max_exac_frequency"),
        ("gnomad", "gnomad_frequency"),
        ("gnomad_max", "max_gnomad_frequency"),
        (
            "gnomad_mt_heteroplasmic",
            "gnomad_mt_heteroplasmic_frequency",
        ),
        ("gnomad_mt_homoplasmic", "gnomad_mt_homoplasmic_frequency"),
        ("gnomad_nhomalt", "gnomad_nhomalt_frequency"),
        ("gnomad_nhomalt_XX", "gnomad_nhomalt_xx_frequency"),
        ("gnomad_nhomalt_XY", "gnomad_nhomalt_xy_frequency"),
        ("gnomad_obs", "gnomad_obs"),
        ("gnomad_wes", "gnomad_wes_frequency"),
        ("gnomad_wes_nhomalt", "gnomad_wes_nhomalt_frequency"),
        ("gnomad_wes_nhomalt_XX", "gnomad_wes_nhomalt_xx_frequency"),
        ("gnomad_wes_nhomalt_XY", "gnomad_wes_nhomalt_xy_frequency"),
        ("gnomad_wes_obs", "gnomad_wes_obs"),
        ("gnomad_wes_max", "max_gnomad_wes_frequency"),
        ("thousand_g", "thousand_genomes_frequency"),
        ("thousand_g_left", "thousand_genomes_frequency_left"),
        ("thousand_g_max", "max_thousand_genomes_frequency"),
        ("thousand_g_right", "thousand_genomes_frequency_right"),
    ];

    for (source, target) in mappings {
        add_frequency_field(variant, frequencies, source, target);
    }

    let sv_fields = [
        ("clingen_cgh_benign", "clingen_cgh_benign"),
        ("clingen_cgh_pathogenic", "clingen_cgh_pathogenic"),
        ("clingen_mip", "clingen_mip"),
        ("clingen_ngi", "clingen_ngi"),
        ("swegen", "swegen"),
        ("decipher", "decipher"),
    ];

    for (source, target) in sv_fields {
        add_frequency_field(variant, frequencies, source, target);
    }

    for key in [
        "swegen_alu",
        "swegen_herv",
        "swegen_l1",
        "swegen_sva",
        "swegen_mei_max",
    ] {
        add_frequency_field(variant, frequencies, key, key);
    }
}

/// Parse a structural-variant or MEI frequency/occurrence value from a VCF record.
///
/// Fields ending in `AF`, `FRQ`, or similar frequency identifiers are parsed
/// as floating-point values. Other fields are parsed as integer occurrence
/// counts. Zero, negative, missing, and placeholder values are ignored.
///
/// Returns the parsed value as BSON, or `None` if no valid value is found.
pub fn parse_sv_frequency(record: &Record, info_key: &str) -> Option<Bson> {
    let key = info_key.as_bytes();

    if info_key.to_uppercase().contains("AF") || info_key.to_uppercase().contains("FRQ") {
        let value = record.info(key).float().ok().flatten()?.first().copied()?;

        if value > 0.0 {
            Some(Bson::Double(value as f64))
        } else {
            None
        }
    } else {
        let value = record
            .info(key)
            .integer()
            .ok()
            .flatten()?
            .first()
            .copied()?;

        if value > 0 {
            Some(Bson::Int32(value))
        } else {
            None
        }
    }
}

/// Update an SV frequency document with the first available value from a
/// list of INFO keys.
///
/// The first key that contains a valid frequency or occurrence value is
/// stored under `new_key`. If none of the keys contain a valid value, the
/// document is left unchanged.
pub fn update_sv_frequency_from_vcf(
    frequencies: &mut Document,
    record: &Record,
    key_list: &[&str],
    new_key: &str,
) {
    for key in key_list {
        if let Some(value) = parse_sv_frequency(record, key) {
            frequencies.insert(new_key, value);
            break;
        }
    }
}

pub fn parse_sv_frequencies(record: &Record) -> Document {
    let mut frequencies = Document::new();

    update_sv_frequency_from_vcf(
        &mut frequencies,
        record,
        CLINGEN_BENIGN_KEYS,
        "clingen_cgh_benign",
    );
    update_sv_frequency_from_vcf(
        &mut frequencies,
        record,
        CLINGEN_PATHOGENIC_KEYS,
        "clingen_cgh_pathogenic",
    );
    update_sv_frequency_from_vcf(&mut frequencies, record, CLINGEN_NGI_KEYS, "clingen_ngi");
    update_sv_frequency_from_vcf(&mut frequencies, record, SWEGEN_KEYS, "swegen");
    update_sv_frequency_from_vcf(&mut frequencies, record, DECIPHER_KEYS, "decipher");
    update_sv_frequency_from_vcf(&mut frequencies, record, CG_KEYS, "clingen_mip");
    update_sv_frequency_from_vcf(&mut frequencies, record, &["colorsdb_af"], "colorsdb_af");

    frequencies
}

/// Parse MEI-specific frequency annotations from a VCF record.
///
/// Extracts SWEGEN frequencies for ALU, HERV, L1, and SVA insertions.
/// If any frequencies are found, also stores the maximum frequency as
/// `swegen_mei_max`.
///
/// Returns an empty document if no MEI frequencies are available.
pub fn parse_mei_frequencies(record: &Record) -> Document {
    let mut frequencies = Document::new();

    update_sv_frequency_from_vcf(&mut frequencies, record, SWEGEN_ALU_KEYS, "swegen_alu");

    update_sv_frequency_from_vcf(&mut frequencies, record, SWEGEN_HERV_KEYS, "swegen_herv");

    update_sv_frequency_from_vcf(&mut frequencies, record, SWEGEN_L1_KEYS, "swegen_l1");

    update_sv_frequency_from_vcf(&mut frequencies, record, SWEGEN_SVA_KEYS, "swegen_sva");

    let max_frequency = frequencies
        .values()
        .filter_map(|value| match value {
            Bson::Double(value) => Some(*value),
            Bson::Int32(value) => Some(*value as f64),
            Bson::Int64(value) => Some(*value as f64),
            _ => None,
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap());

    if let Some(max_frequency) = max_frequency {
        frequencies.insert("swegen_mei_max", max_frequency);
    }

    frequencies
}
