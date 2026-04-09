use std::any::Any;

use crate::error::Result;
use crate::modules::QcModule;
use crate::report::ReportArchive;
use crate::sequence::phred::PhredEncoding;
use crate::sequence::Sequence;

pub struct BasicStats {
    name: Option<String>,
    actual_count: u64,
    filtered_count: u64,
    min_length: usize,
    max_length: usize,
    total_bases: u64,
    g_count: u64,
    c_count: u64,
    a_count: u64,
    t_count: u64,
    n_count: u64,
    lowest_char: u8,
    file_type: Option<String>,
}

impl BasicStats {
    pub fn new() -> Self {
        Self {
            name: None,
            actual_count: 0,
            filtered_count: 0,
            min_length: 0,
            max_length: 0,
            total_bases: 0,
            g_count: 0,
            c_count: 0,
            a_count: 0,
            t_count: 0,
            n_count: 0,
            lowest_char: 126,
            file_type: None,
        }
    }

    pub fn set_file_name(&mut self, name: &str) {
        let n = name.strip_prefix("stdin:").unwrap_or(name);
        self.name = Some(n.to_string());
    }

    pub fn format_length(original_length: u64) -> String {
        let (length, unit) = if original_length >= 1_000_000_000 {
            (original_length as f64 / 1_000_000_000.0, " Gbp")
        } else if original_length >= 1_000_000 {
            (original_length as f64 / 1_000_000.0, " Mbp")
        } else if original_length >= 1_000 {
            (original_length as f64 / 1_000.0, " kbp")
        } else {
            (original_length as f64, " bp")
        };

        // Match Java's formatting: show one decimal digit if non-zero, no trailing .0
        let raw = format!("{}", length);
        let chars: Vec<char> = raw.chars().collect();

        // Find the dot position
        let dot_pos = chars.iter().position(|&c| c == '.');
        let result = if let Some(dot) = dot_pos {
            if dot + 1 < chars.len() && chars[dot + 1] != '0' {
                // Keep one decimal place
                chars[..=dot + 1].iter().collect::<String>()
            } else {
                // No decimal
                chars[..dot].iter().collect::<String>()
            }
        } else {
            raw
        };

        format!("{}{}", result, unit)
    }
}

impl QcModule for BasicStats {
    fn name(&self) -> &str {
        "Basic Statistics"
    }

    fn description(&self) -> &str {
        "Calculates some basic statistics about the file"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        false
    }

    fn ignore_in_report(&self) -> bool {
        false
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if self.name.is_none() {
            self.set_file_name(&seq.file_name);
        }

        if seq.is_filtered {
            self.filtered_count += 1;
            return;
        }

        self.actual_count += 1;
        self.total_bases += seq.sequence.len() as u64;

        if self.file_type.is_none() {
            self.file_type = Some(if seq.colorspace.is_some() {
                "Colorspace converted to bases".to_string()
            } else {
                "Conventional base calls".to_string()
            });
        }

        let seq_len = seq.sequence.len();
        if self.actual_count == 1 {
            self.min_length = seq_len;
            self.max_length = seq_len;
        } else {
            if seq_len < self.min_length {
                self.min_length = seq_len;
            }
            if seq_len > self.max_length {
                self.max_length = seq_len;
            }
        }

        for c in seq.sequence.bytes() {
            match c {
                b'G' => self.g_count += 1,
                b'A' => self.a_count += 1,
                b'T' => self.t_count += 1,
                b'C' => self.c_count += 1,
                b'N' => self.n_count += 1,
                _ => {}
            }
        }

        for c in seq.quality.bytes() {
            if c < self.lowest_char {
                self.lowest_char = c;
            }
        }
    }

    fn reset(&mut self) {
        self.min_length = 0;
        self.max_length = 0;
        self.g_count = 0;
        self.c_count = 0;
        self.a_count = 0;
        self.t_count = 0;
        self.n_count = 0;
    }

    fn raises_error(&mut self) -> bool {
        false
    }

    fn raises_warning(&mut self) -> bool {
        false
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, other: Box<dyn Any + Send>) {
        if let Ok(other) = other.downcast::<Self>() {
            self.actual_count += other.actual_count;
            self.filtered_count += other.filtered_count;
            self.total_bases += other.total_bases;
            self.g_count += other.g_count;
            self.c_count += other.c_count;
            self.a_count += other.a_count;
            self.t_count += other.t_count;
            self.n_count += other.n_count;

            if other.min_length > 0 && (self.min_length == 0 || other.min_length < self.min_length)
            {
                self.min_length = other.min_length;
            }
            if other.max_length > self.max_length {
                self.max_length = other.max_length;
            }
            if other.lowest_char < self.lowest_char {
                self.lowest_char = other.lowest_char;
            }
            if self.name.is_none() {
                self.name = other.name;
            }
            if self.file_type.is_none() {
                self.file_type = other.file_type;
            }
        }
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        let encoding = PhredEncoding::from_lowest_char(self.lowest_char);
        let gc_percent = if self.a_count + self.t_count + self.g_count + self.c_count > 0 {
            (self.g_count + self.c_count) * 100
                / (self.a_count + self.t_count + self.g_count + self.c_count)
        } else {
            0
        };

        let seq_length = if self.min_length == self.max_length {
            format!("{}", self.min_length)
        } else {
            format!("{}-{}", self.min_length, self.max_length)
        };

        let name = self.name.clone().unwrap_or_else(|| "unknown".to_string());

        // Write table to data document
        let data = &mut report.data;
        data.push_str("#Measure\tValue\n");
        data.push_str(&format!("Filename\t{}\n", name));
        data.push_str(&format!(
            "File type\t{}\n",
            self.file_type
                .as_deref()
                .unwrap_or("Conventional base calls")
        ));
        data.push_str(&format!("Encoding\t{}\n", encoding));
        data.push_str(&format!("Total Sequences\t{}\n", self.actual_count));
        data.push_str(&format!(
            "Total Bases\t{}\n",
            Self::format_length(self.total_bases)
        ));
        data.push_str(&format!(
            "Sequences flagged as poor quality\t{}\n",
            self.filtered_count
        ));
        data.push_str(&format!("Sequence length\t{}\n", seq_length));
        data.push_str(&format!("%GC\t{}\n", gc_percent));

        // Write HTML table
        let html = &mut report.html_body;
        html.push_str("<table><thead><tr><th>Measure</th><th>Value</th></tr></thead><tbody>");
        for (measure, value) in &[
            ("Filename", name.clone()),
            (
                "File type",
                self.file_type
                    .clone()
                    .unwrap_or_else(|| "Conventional base calls".into()),
            ),
            ("Encoding", encoding.to_string()),
            ("Total Sequences", self.actual_count.to_string()),
            ("Total Bases", Self::format_length(self.total_bases)),
            (
                "Sequences flagged as poor quality",
                self.filtered_count.to_string(),
            ),
            ("Sequence length", seq_length.clone()),
            ("%GC", gc_percent.to_string()),
        ] {
            html.push_str(&format!("<tr><td>{}</td><td>{}</td></tr>", measure, value));
        }
        html.push_str("</tbody></table>");

        Ok(())
    }
}
