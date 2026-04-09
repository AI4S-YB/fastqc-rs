use std::any::Any;

use crate::error::Result;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;

pub struct SequenceLengthDistribution {
    length_counts: Vec<f64>,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl SequenceLengthDistribution {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            length_counts: Vec::new(),
            ignore: module_config.get_param("sequence_length", "ignore") > 0.0,
            warn_threshold: module_config.get_param("sequence_length", "warn"),
            error_threshold: module_config.get_param("sequence_length", "error"),
        }
    }

    fn min_length(&self) -> usize {
        for (i, &count) in self.length_counts.iter().enumerate() {
            if count > 0.0 {
                return i;
            }
        }
        0
    }

    fn max_length(&self) -> usize {
        for (i, &count) in self.length_counts.iter().enumerate().rev() {
            if count > 0.0 {
                return i;
            }
        }
        0
    }
}

impl QcModule for SequenceLengthDistribution {
    fn name(&self) -> &str {
        "Sequence Length Distribution"
    }

    fn description(&self) -> &str {
        "Shows the distribution of sequence lengths"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        let len = seq.sequence.len();
        if self.length_counts.len() <= len {
            self.length_counts.resize(len + 1, 0.0);
        }
        self.length_counts[len] += 1.0;
    }

    fn reset(&mut self) {
        self.length_counts.clear();
    }

    fn raises_error(&mut self) -> bool {
        // Error if any sequences with length 0
        if !self.length_counts.is_empty() && self.length_counts[0] > 0.0 {
            return true;
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        // Warn if variable length
        let min = self.min_length();
        let max = self.max_length();
        min != max
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, other: Box<dyn Any + Send>) {
        if let Ok(other) = other.downcast::<Self>() {
            let max_len = self.length_counts.len().max(other.length_counts.len());
            self.length_counts.resize(max_len, 0.0);

            for (i, &v) in other.length_counts.iter().enumerate() {
                self.length_counts[i] += v;
            }
        }
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        let min = self.min_length();
        let max = self.max_length();

        // Build display data - group lengths if range is too wide
        let range = max - min + 1;
        let interval = if range > 75 {
            let base_values = [2, 5, 10];
            let mut interval = 1;
            for multiplier in [1, 10, 100, 1000] {
                for &base in &base_values {
                    if range / (base * multiplier) <= 50 {
                        interval = base * multiplier;
                        break;
                    }
                }
                if range / interval <= 50 {
                    break;
                }
            }
            interval
        } else {
            1
        };

        let mut display_lengths: Vec<String> = Vec::new();
        let mut display_counts: Vec<f64> = Vec::new();

        let mut pos = min;
        while pos <= max {
            let end = (pos + interval - 1).min(max);
            let mut count = 0.0;
            for i in pos..=end {
                if i < self.length_counts.len() {
                    count += self.length_counts[i];
                }
            }
            if pos == end {
                display_lengths.push(format!("{}", pos));
            } else {
                display_lengths.push(format!("{}-{}", pos, end));
            }
            display_counts.push(count);
            pos = end + 1;
        }

        let max_count = display_counts.iter().cloned().fold(0.0f64, f64::max);

        let svg = crate::graphs::line_graph::render_line_graph(
            &[&display_counts],
            0.0,
            max_count.max(1.0),
            "Sequence Length (bp)",
            &["Sequence Length"],
            &display_lengths,
            "Distribution of sequence lengths over all sequences",
            800,
            600,
        );
        report.add_image("sequence_length_distribution.png", &svg);

        let data = &mut report.data;
        data.push_str("#Length\tCount\n");
        for i in 0..display_lengths.len() {
            data.push_str(&format!(
                "{}\t{}\n",
                display_lengths[i],
                crate::format_double(display_counts[i])
            ));
        }

        Ok(())
    }
}
