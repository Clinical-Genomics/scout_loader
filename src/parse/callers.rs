use crate::parse::info::parse_info_string;
use mongodb::bson::{Bson, Document};
use rust_htslib::bcf::Record;

use crate::models::variant::VariantCategory;

const PASS_STATUS: &str = "Pass";
const FILTERED_STATUS: &str = "Filtered";

const CALLERS: &[(&str, &[&str])] = &[
    (
        "snv",
        &[
            "bcftools",
            "deepvariant",
            "freebayes",
            "gatk",
            "mutect2",
            "samtools",
            "mitorsaw",
        ],
    ),
    (
        "cancer",
        &[
            "freebayes",
            "gatk",
            "mutect",
            "pindel",
            "tnscope",
            "tnscope_umi",
            "vardict",
        ],
    ),
    (
        "cancer_sv",
        &[
            "ascat", "cnvkit", "dellysv", "dellycnv", "gatk", "igh_dux4", "manta", "tiddit",
        ],
    ),
    ("mei", &["retroseq"]),
    (
        "sv",
        &[
            "cnvnator",
            "cnvpytor",
            "delly",
            "gatk",
            "gcnvcaller",
            "hificnv",
            "manta",
            "mitosalt",
            "mitorsaw",
            "severus",
            "sniffles",
            "tiddit",
        ],
    ),
    ("str", &["expansionhunter", "trgt"]),
    ("fusion", &["arriba", "fusioncatcher", "starfusion"]),
];

fn callers_for_category(category: VariantCategory) -> &'static [&'static str] {
    CALLERS
        .iter()
        .find(|(name, _)| *name == category.as_str())
        .map(|(_, callers)| *callers)
        .unwrap_or(&[])
}

/// Parse callers from a `FOUND_IN` INFO annotation.
///
/// Callers are comma-separated, with the caller ID before the first `|`.
/// If the variant is filtered and only one caller is present, that caller
/// receives the specific filter status. If multiple callers are present,
/// all callers are marked as `Filtered`.
fn get_callers_from_found_in(
    callers: &mut Document,
    found_in: &str,
    filter_status: Option<&str>,
) -> Document {
    let found_ins: Vec<&str> = found_in.split(',').collect();

    let call_status = caller_status(filter_status, found_ins.len());

    for found_in in found_ins {
        let called_by = found_in.split('|').next().unwrap_or("");

        if callers.contains_key(called_by) {
            callers.insert(called_by, Bson::String(call_status.clone()));
        }
    }

    callers.clone()
}

/// Parse callers from an `svdb_origin` INFO annotation.
///
/// Callers are pipe-separated. If the variant is filtered and only one
/// caller is present, that caller receives the specific filter status.
/// If multiple callers are present, all callers are marked as `Filtered`.
fn get_callers_from_svdb_origin(
    callers: &mut Document,
    svdb_origin: &str,
    filter_status: Option<&str>,
) -> Document {
    let svdb_callers: Vec<&str> = svdb_origin.split('|').collect();

    let call_status = caller_status(filter_status, svdb_callers.len());

    for called_by in svdb_callers {
        if callers.contains_key(called_by) {
            callers.insert(called_by, Bson::String(call_status.clone()));
        }
    }

    callers.clone()
}

/// Parse callers from a GATK `set` INFO annotation.
///
/// The `set` field contains dash-separated caller information. `Intersection`
/// marks all callers as passing, while `FilteredInAll` marks all callers as
/// filtered. Individual `filterIn` entries identify the callers that passed
/// or were filtered. Otherwise, a caller ID directly marks that caller as
/// passing.
fn get_callers_from_set(
    callers: &mut Document,
    info_set: &str,
    filter_status: Option<&str>,
) -> Document {
    let calls: Vec<&str> = info_set.split('-').collect();

    let mut call_status = caller_status(filter_status, calls.len());

    for call in calls {
        if call == "FilteredInAll" {
            for caller in callers.iter_mut() {
                *caller.1 = Bson::String(FILTERED_STATUS.to_string());
            }
            return callers.clone();
        }

        if call == "Intersection" {
            for caller in callers.iter_mut() {
                *caller.1 = Bson::String(PASS_STATUS.to_string());
            }
            return callers.clone();
        }

        if call.contains("filterIn") {
            if !call_status.starts_with("Filtered") {
                call_status = FILTERED_STATUS.to_string();
            }

            for (caller, status) in callers.iter_mut() {
                if call.contains(caller) {
                    *status = Bson::String(call_status.clone());
                }
            }
        } else if callers.contains_key(call) {
            callers.insert(call, Bson::String(PASS_STATUS.to_string()));
        }
    }

    callers.clone()
}

/// Determine the caller status from the variant filter and number of callers.
///
/// A passing variant gets `Pass`. For a filtered variant, a single caller
/// gets the specific filter status while multiple callers get `Filtered`.
fn caller_status(filter_status: Option<&str>, caller_count: usize) -> String {
    match filter_status {
        None => PASS_STATUS.to_string(),
        Some(status) if caller_count == 1 => {
            format!("Filtered - {}", status.replace(';', " - "))
        }
        Some(_) => FILTERED_STATUS.to_string(),
    }
}

/// Set GATK as the caller for an SNV when no caller information was found.
///
/// This is a fallback for older MIP versions where the GATK caller was not
/// explicitly recorded. The GATK status is `Pass` for an unfiltered variant,
/// or `Filtered - <status>` when the variant has a filter status.
fn get_callers_gatk_snv_fallback(callers: &mut Document, filter_status: Option<&str>) -> Document {
    let status = match filter_status {
        None => PASS_STATUS.to_string(),
        Some(status) => {
            format!("Filtered - {}", status.replace(';', " - "))
        }
    };

    if callers.contains_key("gatk") {
        callers.insert("gatk", Bson::String(status));
    }

    callers.clone()
}

/// Parse variant caller statuses from VCF INFO annotations.
///
/// Caller information is checked in the following order:
/// 1. `FOUND_IN`
/// 2. `svdb_origin`
/// 3. `set`
/// 4. GATK fallback for SNV variants
///
/// Only callers with an assigned status are included in the returned document.
pub fn parse_callers(record: &Record, category: VariantCategory, filters: &[String]) -> Document {
    let relevant_callers = callers_for_category(category);
    let mut callers = Document::new();

    for caller in relevant_callers {
        callers.insert(*caller, Bson::Null);
    }

    let filter_status = if filters.iter().any(|filter| filter == "PASS") {
        None
    } else {
        Some(filters.join(" - "))
    };

    if let Some(found_in) = parse_info_string(record, b"FOUND_IN") {
        get_callers_from_found_in(&mut callers, &found_in, filter_status.as_deref());
    } else if let Some(svdb_origin) = parse_info_string(record, b"svdb_origin") {
        get_callers_from_svdb_origin(&mut callers, &svdb_origin, filter_status.as_deref());
    } else if let Some(info_set) = parse_info_string(record, b"set") {
        get_callers_from_set(&mut callers, &info_set, filter_status.as_deref());
    } else if category == VariantCategory::Snv {
        get_callers_gatk_snv_fallback(&mut callers, filter_status.as_deref());
    }

    let null_callers: Vec<String> = callers
        .iter()
        .filter(|(_, value)| matches!(value, Bson::Null))
        .map(|(key, _)| key.to_string())
        .collect();

    for caller in null_callers {
        callers.remove(&caller);
    }

    callers
}
