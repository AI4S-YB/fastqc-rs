use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

use serde_json::{json, Value};

use crate::config::Config;
use crate::modules::QcModule;

/// Canonical module ids in the order mandated by the schema (prefixItems).
/// The module list produced by `modules::create_module_list` already follows
/// this order — we keep the names here only as a defensive assertion.
const CANONICAL_IDS: [&str; 12] = [
    "basic_statistics",
    "per_base_sequence_quality",
    "per_tile_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content",
    "kmer_content",
];

pub fn build_report(
    source_path: &Path,
    file_name: &str,
    detected_format: &str,
    modules: &mut [Box<dyn QcModule>],
    config: &Config,
) -> Value {
    assert_eq!(
        modules.len(),
        12,
        "json report expects 12 canonical modules"
    );

    let mut module_values: Vec<Value> = Vec::with_capacity(12);
    let mut pass = 0u64;
    let mut warn = 0u64;
    let mut fail = 0u64;
    let mut na = 0u64;

    let mut basic_statistics_data: Option<Value> = None;

    for (idx, module) in modules.iter_mut().enumerate() {
        let expected_id = CANONICAL_IDS[idx];
        let actual_id = module.module_id().to_string();
        let name = module.name().to_string();
        let description = module.description().to_string();
        assert_eq!(
            actual_id, expected_id,
            "module {} ordered unexpectedly (got {})",
            idx, actual_id
        );

        let in_report = !module.ignore_in_report();
        let (status, omitted_reason) = if !in_report {
            ("not_applicable", Some(Value::String("module ignored".into())))
        } else if module.raises_error() {
            ("fail", None)
        } else if module.raises_warning() {
            ("warn", None)
        } else {
            ("pass", None)
        };

        match status {
            "pass" => pass += 1,
            "warn" => warn += 1,
            "fail" => fail += 1,
            _ => na += 1,
        }

        let data = module.json_data(config);
        if idx == 0 {
            basic_statistics_data = Some(data.clone());
        }

        let mut obj = serde_json::Map::new();
        obj.insert("id".into(), Value::String(actual_id));
        obj.insert("name".into(), Value::String(name));
        obj.insert("description".into(), Value::String(description));
        obj.insert("status".into(), Value::String(status.to_string()));
        obj.insert("in_report".into(), Value::Bool(in_report));
        obj.insert(
            "omitted_reason".into(),
            omitted_reason.unwrap_or(Value::Null),
        );
        if let Some(t) = module.json_thresholds() {
            obj.insert("thresholds".into(), t);
        }
        obj.insert("data".into(), data);
        module_values.push(Value::Object(obj));
    }

    let overall_status = if fail > 0 {
        "fail"
    } else if warn > 0 {
        "warn"
    } else {
        "pass"
    };

    let bs = basic_statistics_data.unwrap_or_else(|| json!({}));
    let summary_basic = json!({
        "total_sequences": bs.get("total_sequences").cloned().unwrap_or(json!(0)),
        "filtered_sequences": bs.get("filtered_sequences").cloned().unwrap_or(json!(0)),
        "total_bases": bs.get("total_bases").cloned().unwrap_or(json!(0)),
        "min_sequence_length": bs
            .pointer("/sequence_length/min")
            .cloned()
            .unwrap_or(json!(0)),
        "max_sequence_length": bs
            .pointer("/sequence_length/max")
            .cloned()
            .unwrap_or(json!(0)),
        "gc_percent": bs.get("gc_percent").cloned().unwrap_or(json!(0)),
    });

    json!({
        "schema_version": "1.0.0",
        "generated_by": {
            "tool": "fastqc-rs",
            "version": crate::VERSION,
        },
        "generated_at": rfc3339_now(),
        "input": {
            "file_name": file_name,
            "source_path": source_path.to_string_lossy(),
            "detected_format": detected_format,
        },
        "summary": {
            "overall_status": overall_status,
            "module_counts": {
                "pass": pass,
                "warn": warn,
                "fail": fail,
                "not_applicable": na,
            },
            "basic_statistics": summary_basic,
        },
        "modules": module_values,
    })
}

fn rfc3339_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    format_rfc3339(secs)
}

fn format_rfc3339(unix_seconds: i64) -> String {
    let days = unix_seconds.div_euclid(86_400);
    let sod = unix_seconds.rem_euclid(86_400) as u32;
    let (year, month, day) = civil_from_days(days);
    let h = sod / 3600;
    let m = (sod % 3600) / 60;
    let s = sod % 60;
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, h, m, s
    )
}

/// Howard Hinnant's civil_from_days algorithm: converts days since Unix epoch
/// (1970-01-01) into a (year, month, day) triple in the proleptic Gregorian
/// calendar.
fn civil_from_days(z: i64) -> (i32, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let year = (y + if m <= 2 { 1 } else { 0 }) as i32;
    (year, m, d)
}
