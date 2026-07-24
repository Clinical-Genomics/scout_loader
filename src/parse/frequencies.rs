use mongodb::bson::{Bson, Document};
use rust_htslib::bcf::Record;

pub const EXAC_KEYS: &[&str] = &[
    "EXACAF",
];

pub const EXAC_MAX_KEYS: &[&str] = &[
    "EXACMAXAF",
    "EXACMAX_AF",
];

pub const GNOMAD_INFO_KEYS: &[&str] = &[
    "GNOMADAF",
    "GNOMAD_AF",
    "gnomADg_AF",
    "gnomad_svAF",
    "gnomad_af",
];

pub const SWEGEN_KEYS: &[&str] = &[
    "SWEGENAF",
    "swegenAF",
    "swegen",
];

pub const GNOMAD_INFO_MAX_KEYS: &[&str] = &[
    "gnomADg_AF_POPMAX",
    "GNOMADAF_popmax",
    "GNOMADAF_POPMAX",
    "GNOMADAF_MAX",
    "gnomad_popmax_af",
];

pub const THOUSAND_GENOMES_KEYS: &[&str] = &[
    "1000GAF",
];

pub const THOUSAND_GENOMES_MAX_KEYS: &[&str] = &[
    "1000G_MAX_AF",
];

/// Parse a frequency value from a VCF INFO field.
///
/// Returns `None` if the field is missing or contains a placeholder value
/// (`.`, `0`, `-1`).
///
/// Returns the frequency as `f64` otherwise.
fn parse_frequency(record: &Record, key: &[u8]) -> Option<f64> {
    record
        .info(key)
        .float()
        .ok()
        .flatten()
        .and_then(|values| values.first().copied())
        .filter(|v| *v != 0.0 && *v != -1.0)
        .map(|v| v as f64)
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
            frequency.insert(
                new_key,
                Bson::Double(value),
            );
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
pub fn parse_frequencies(
    record: &Record,
    transcripts: &[Document],
) -> Document {
    let mut frequencies = Document::new();

    let frequency_fields = [
        (&EXAC_KEYS[..], "exac"),
        (&EXAC_MAX_KEYS[..], "exac_max"),
        (&GNOMAD_INFO_KEYS[..], "gnomad"),
        (&SWEGEN_KEYS[..], "swegen"),
        (&GNOMAD_INFO_MAX_KEYS[..], "gnomad_max"),
        (&THOUSAND_GENOMES_KEYS[..], "thousand_g"),
        (&THOUSAND_GENOMES_MAX_KEYS[..], "thousand_g_max"),
        (&["GNOMAD_MT_AF_HOM"][..], "gnomad_mt_homoplasmic"),
        (&["GNOMAD_MT_AF_HET"][..], "gnomad_mt_heteroplasmic"),
        (&["left_1000GAF"][..], "thousand_g_left"),
        (&["right_1000GAF"][..], "thousand_g_right"),
        (&["colorsdb_af"][..], "colorsdb_af"),
    ];

    for (keys, name) in frequency_fields {
        update_frequency_from_vcf(
            &mut frequencies,
            record,
            keys,
            name,
        );
    }

    if frequencies.is_empty() {
        update_frequency_from_transcript(
            &mut frequencies,
            transcripts,
        );
    }

    frequencies
}


/// Update frequencies from VEP transcript annotations.
///
/// Searches transcript-level annotations for population frequencies and adds
/// them to the frequency document.
pub fn update_frequency_from_transcript(
    frequencies: &mut Document,
    transcripts: &[Document],
) {
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
            if let Some(value) = transcript.get(transcript_key) {
                if !matches!(value, Bson::Null) {
                    frequencies.insert(
                        frequency_key,
                        value.clone(),
                    );
                }
            }
        }
    }
}