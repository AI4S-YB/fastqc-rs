use std::collections::HashMap;

use crate::error::Result;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::phred::PhredEncoding;
use crate::sequence::Sequence;

pub struct PerSequenceQualityScores {
    quality_distribution: HashMap<i32, f64>,
    lowest_char: u8,
    calculated: bool,
    most_frequent_score: i32,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl PerSequenceQualityScores {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            quality_distribution: HashMap::new(),
            lowest_char: 255,
            calculated: false,
            most_frequent_score: 0,
            ignore: module_config.get_param("quality_sequence", "ignore") > 0.0,
            warn_threshold: module_config.get_param("quality_sequence", "warn"),
            error_threshold: module_config.get_param("quality_sequence", "error"),
        }
    }

    fn calculate(&mut self) {
        if self.calculated {
            return;
        }
        let mut max_count = 0.0f64;
        self.most_frequent_score = 0;
        for (&score, &count) in &self.quality_distribution {
            if count > max_count {
                max_count = count;
                self.most_frequent_score = score;
            }
        }
        self.calculated = true;
    }
}

impl QcModule for PerSequenceQualityScores {
    fn name(&self) -> &str {
        "Per sequence quality scores"
    }

    fn description(&self) -> &str {
        "Shows the distribution of mean quality scores across all sequences"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore || self.quality_distribution.is_empty()
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.calculated = false;
        let qual = seq.quality.as_bytes();
        if qual.is_empty() {
            return;
        }

        for &c in qual {
            if c < self.lowest_char {
                self.lowest_char = c;
            }
        }

        let total: u64 = qual.iter().map(|&c| c as u64).sum();
        let avg = (total as f64 / qual.len() as f64).round() as i32;

        *self.quality_distribution.entry(avg).or_insert(0.0) += 1.0;
    }

    fn reset(&mut self) {
        self.quality_distribution.clear();
        self.lowest_char = 255;
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        self.calculate();
        let encoding = PhredEncoding::from_lowest_char(self.lowest_char);
        let adjusted = self.most_frequent_score - encoding.offset() as i32;
        // Java: error threshold is 27, fires when mode < 27
        adjusted < self.error_threshold as i32
    }

    fn raises_warning(&mut self) -> bool {
        self.calculate();
        let encoding = PhredEncoding::from_lowest_char(self.lowest_char);
        let adjusted = self.most_frequent_score - encoding.offset() as i32;
        // Java: warn threshold is 20, fires when mode < 20
        adjusted < self.warn_threshold as i32
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate();
        let encoding = PhredEncoding::from_lowest_char(self.lowest_char);
        let offset = encoding.offset() as i32;

        // Build sorted data
        let mut scores: Vec<i32> = self.quality_distribution.keys().copied().collect();
        scores.sort();

        let adjusted_scores: Vec<f64> = scores.iter().map(|&s| (s - offset) as f64).collect();
        let counts: Vec<f64> = scores.iter().map(|&s| self.quality_distribution[&s]).collect();
        let x_labels: Vec<String> = adjusted_scores.iter().map(|s| format!("{}", *s as i32)).collect();

        let max_count = counts.iter().cloned().fold(0.0f64, f64::max);

        let svg = crate::graphs::line_graph::render_line_graph(
            &[&counts],
            0.0,
            max_count,
            "Mean Sequence Quality (Phred Score)",
            &["Average Quality per read"],
            &x_labels,
            "Quality score distribution over all sequences",
            800,
            600,
        );
        report.add_image("per_sequence_quality.png", &svg);

        let data = &mut report.data;
        data.push_str("#Quality\tCount\n");
        for &score in &scores {
            let adjusted = score - offset;
            data.push_str(&format!("{}\t{}\n", adjusted, crate::format_double(self.quality_distribution[&score])));
        }

        Ok(())
    }
}
