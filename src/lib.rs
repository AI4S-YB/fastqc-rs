pub const VERSION: &str = "0.12.1";

pub mod config;
pub mod error;
pub mod sequence;
pub mod modules;
pub mod graphs;
pub mod stats;
pub mod report;
pub mod analysis;

/// Format f64 like Java's Double.toString(): always shows `.0` for integer values.
/// e.g., 9.0 -> "9.0", 3.14 -> "3.14", NaN -> "NaN"
pub fn format_double(v: f64) -> String {
    if v.is_nan() {
        return "NaN".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 {
            "Infinity".to_string()
        } else {
            "-Infinity".to_string()
        };
    }
    let s = format!("{}", v);
    if s.contains('.') || s.contains('E') || s.contains('e') {
        s
    } else {
        format!("{}.0", s)
    }
}
