use mongodb::bson::{Document};
use rust_htslib::bcf::Record;
use crate::parse::info::{parse_info_int, parse_info_float};

/// Add LoqusDB archive observation frequencies, counts, and metadata to a variant.
///
/// Parses germline and cancer observation fields from VCF INFO annotations and
/// adds local archive metadata extracted from the VCF header.
///
/// Missing values and sentinel values (`-1`) are ignored.
///
/// Germline annotations:
/// - Obs / clinical_genomics_loqusObs / clin_obs -> `local_obs_old`
/// - Hom -> `local_obs_hom_old`
/// - clinical_genomics_loqusFrq / Frq -> `local_obs_old_freq`
///
/// Cancer annotations:
/// - Cancer_Germline
/// - Cancer_Somatic
/// - Cancer_Somatic_Panel
///
/// Each cancer annotation stores observation count, homozygous count,
/// and frequency when available.
///
/// Header metadata:
/// - Description -> `local_obs_old_desc`
/// - Date -> `local_obs_old_date`
/// - NrCases -> `local_obs_old_nr_cases`
pub fn add_loqus_archive_frequencies(
    record: &Record,
    variant: &mut Document,
    local_archive_info: Option<&Document>,
) {
    // Local archive metadata from VCF header
    if let Some(info) = local_archive_info {
        if let Ok(value) = info.get_str("Description") {
            variant.insert("local_obs_old_desc", value);
        }

        if let Ok(value) = info.get_str("Date") {
            variant.insert("local_obs_old_date", value);
        }

        if let Ok(value) = info.get_i32("NrCases") {
            variant.insert("local_obs_old_nr_cases", value);
        }
    }

    // Germline observations (SNVs and SVs)
    let obs = parse_info_int(record, b"Obs")
        .filter(|&v| v != -1)
        .or_else(|| {
            parse_info_int(record, b"clinical_genomics_loqusObs")
                .filter(|&v| v != -1)
        })
        .or_else(|| {
            parse_info_int(record, b"clin_obs")
                .filter(|&v| v != -1)
        });

    if let Some(value) = obs {
        variant.insert("local_obs_old", value);
    }

    if let Some(value) = parse_info_int(record, b"Hom").filter(|&v| v != -1) {
        variant.insert("local_obs_hom_old", value);
    }

    if let Some(value) = parse_info_float(record, b"clinical_genomics_loqusFrq")
        .or_else(|| parse_info_float(record, b"Frq"))
    {
        variant.insert("local_obs_old_freq", value);
    }

    // Cancer observations
    let cancer_sources = [
        ("Cancer_Germline", "cancer_germline"),
        ("Cancer_Somatic", "cancer_somatic"),
        ("Cancer_Somatic_Panel", "cancer_somatic_panel"),
    ];

    for (prefix, key) in cancer_sources {
        if let Some(value) = parse_info_int(
            record,
            format!("{prefix}_Obs").as_bytes(),
        )
        .filter(|&v| v != -1)
        {
            variant.insert(format!("local_obs_{key}_old"), value);
        }

        if let Some(value) = parse_info_int(
            record,
            format!("{prefix}_Hom").as_bytes(),
        )
        .filter(|&v| v != -1)
        {
            variant.insert(format!("local_obs_{key}_hom_old"), value);
        }

        if let Some(value) = parse_info_float(
            record,
            format!("{prefix}_Frq").as_bytes(),
        ) {
            variant.insert(format!("local_obs_{key}_old_freq"), value);
        }
    }
}