pub const VERSION: &str = "0.14.0";

pub mod config;
pub mod error;
pub mod sequence;
pub mod modules;
pub mod graphs;
pub mod stats;
pub mod report;
pub mod analysis;
pub mod trim_galore;

/// Format f64 like Java's Double.toString(): always shows `.0` for integer values,
/// and uses scientific notation (uppercase E) for |v| < 1e-3 or |v| >= 1e7.
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
    let abs_v = v.abs();
    if abs_v != 0.0 && (abs_v < 1e-3 || abs_v >= 1e7) {
        return to_java_scientific(v);
    }
    let s = format!("{}", v);
    if s.contains('.') {
        s
    } else {
        format!("{}.0", s)
    }
}

fn to_java_scientific(v: f64) -> String {
    let s = format!("{}", v);
    let (sign, unsigned) = if let Some(stripped) = s.strip_prefix('-') {
        ("-", stripped)
    } else {
        ("", s.as_str())
    };

    if let Some(dot_pos) = unsigned.find('.') {
        let int_part = &unsigned[..dot_pos];
        let frac_part = &unsigned[dot_pos + 1..];

        if int_part == "0" {
            let leading_zeros = frac_part.len() - frac_part.trim_start_matches('0').len();
            let significant = frac_part.trim_start_matches('0');
            let exp = -(leading_zeros as i32 + 1);
            if significant.len() > 1 {
                format!("{}{}.{}E{}", sign, &significant[..1], &significant[1..], exp)
            } else {
                format!("{}{}.0E{}", sign, &significant[..1], exp)
            }
        } else {
            let exp = int_part.len() as i32 - 1;
            let all_digits = format!("{}{}", int_part, frac_part);
            let trimmed = all_digits.trim_end_matches('0');
            if trimmed.len() > 1 {
                format!("{}{}.{}E{}", sign, &trimmed[..1], &trimmed[1..], exp)
            } else {
                format!("{}{}.0E{}", sign, &trimmed[..1], exp)
            }
        }
    } else {
        let trimmed = unsigned.trim_end_matches('0');
        let exp = unsigned.len() as i32 - 1;
        if trimmed.len() > 1 {
            format!("{}{}.{}E{}", sign, &trimmed[..1], &trimmed[1..], exp)
        } else {
            format!("{}{}.0E{}", sign, &trimmed[..1], exp)
        }
    }
}
