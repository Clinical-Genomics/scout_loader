use mongodb::bson::{Bson, Document};
use crate::HashMap;

/// Parse mitochondrial gnomAD allele frequencies from a VEP transcript.
///
/// Extracts homoplasmic and heteroplasmic mitochondrial allele frequencies
/// when available.
pub fn parse_mt_frequencies(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    if let Some(value) = entry
        .get("GNOMAD_MT_AF_HOM")
        .filter(|value| !value.is_empty())
    {
        if let Ok(value) = value.parse::<f64>() {
            transcript.insert(
                "gnomad_mt_homoplasmic",
                Bson::Double(value),
            );
        }
    }

    if let Some(value) = entry
        .get("GNOMAD_MT_AF_HET")
        .filter(|value| !value.is_empty())
    {
        if let Ok(value) = value.parse::<f64>() {
            transcript.insert(
                "gnomad_mt_heteroplasmic",
                Bson::Double(value),
            );
        }
    }
}


/// Parse variant population frequencies from VEP transcript annotations.
///
/// Supports VEP v90+ frequency fields, including 1000 Genomes, gnomAD,
/// and ExAC frequencies. Stores specific allele frequencies and calculates
/// maximum frequencies across population-specific annotations.
pub fn parse_variant_frequencies(
    transcript: &mut Document,
    entry: &HashMap<String, String>,
) {
    let gnomad_csq_keys = [
        "GNOMAD_AF",
        "GNOMADG_AF",
        "GNOMAD_EXOMES_AF",
        "GNOMAD_EXOMES_AF",
    ];

    let thousand_genomes_csq_keys = [
        "AF",
        "1000GAF",
        "1000GP3_AF",
    ];

    let mut thousandg_freqs: Vec<f64> = Vec::new();
    let mut gnomad_freqs: Vec<f64> = Vec::new();

    for (key, value) in entry {
        // Frequency fields end with AF or POPMAX
        if !(key.ends_with("AF") || key.ends_with("POPMAX")) {
            continue;
        }

        if value.is_empty()
            || value == "."
            || value.chars().all(char::is_alphabetic)
        {
            continue;
        }

        let Ok(freq) = value.parse::<f64>() else {
            continue;
        };

        if thousand_genomes_csq_keys.contains(&key.as_str()) {
            transcript.insert(
                "thousand_g_maf",
                Bson::Double(freq),
            );
            continue;
        }

        if gnomad_csq_keys.contains(&key.as_str()) {
            transcript.insert(
                "gnomad_maf",
                Bson::Double(freq),
            );
            continue;
        }

        if key == "EXAC_MAX_AF" {
            transcript.insert(
                "exac_max",
                Bson::Double(freq),
            );
            transcript.insert(
                "exac_maf",
                Bson::Double(freq),
            );
            continue;
        }

        // Remaining population-specific frequencies
        if key.contains("GNOMAD") {
            gnomad_freqs.push(freq);
        } else {
            thousandg_freqs.push(freq);
        }
    }

    if let Some(max) = thousandg_freqs.iter().max_by(|a, b| a.total_cmp(b)) {
        transcript.insert(
            "thousandg_max",
            Bson::Double(*max),
        );
    }

    if let Some(max) = gnomad_freqs.iter().max_by(|a, b| a.total_cmp(b)) {
        transcript.insert(
            "gnomad_max",
            Bson::Double(*max),
        );
    }
}