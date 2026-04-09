use std::any::Any;

use crate::error::Result;
use crate::graphs::base_group;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;

pub struct NContent {
    n_counts: Vec<u64>,
    not_n_counts: Vec<u64>,
    calculated: bool,
    percentages: Vec<f64>,
    x_labels: Vec<String>,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl NContent {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            n_counts: Vec::new(),
            not_n_counts: Vec::new(),
            calculated: false,
            percentages: Vec::new(),
            x_labels: Vec::new(),
            ignore: module_config.get_param("n_content", "ignore") > 0.0,
            warn_threshold: module_config.get_param("n_content", "warn"),
            error_threshold: module_config.get_param("n_content", "error"),
        }
    }

    fn calculate(&mut self, config: &crate::config::Config) {
        if self.calculated {
            return;
        }

        let groups = base_group::make_base_groups(self.n_counts.len(), config);
        let mut percentages = Vec::with_capacity(groups.len());
        let mut x_labels = Vec::with_capacity(groups.len());

        for group in &groups {
            x_labels.push(group.to_string());
            let mut n_total = 0u64;
            let mut not_n_total = 0u64;

            for i in (group.lower_count() - 1)..group.upper_count().min(self.n_counts.len()) {
                n_total += self.n_counts[i];
                not_n_total += self.not_n_counts[i];
            }

            let total = n_total + not_n_total;
            if total > 0 {
                percentages.push(n_total as f64 * 100.0 / total as f64);
            } else {
                percentages.push(0.0);
            }
        }

        self.percentages = percentages;
        self.x_labels = x_labels;
        self.calculated = true;
    }

    fn check_threshold(&mut self, threshold: f64) -> bool {
        if !self.calculated {
            let config = crate::config::Config::default_config();
            self.calculate(&config);
        }
        self.percentages.iter().any(|&p| p > threshold)
    }
}

impl QcModule for NContent {
    fn name(&self) -> &str {
        "Per base N content"
    }

    fn description(&self) -> &str {
        "Shows the percentage of N bases at each position"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.calculated = false;
        let seq_bytes = seq.sequence.as_bytes();
        let len = seq_bytes.len();

        if self.n_counts.len() < len {
            self.n_counts.resize(len, 0);
            self.not_n_counts.resize(len, 0);
        }

        for (i, &c) in seq_bytes.iter().enumerate() {
            if c == b'N' {
                self.n_counts[i] += 1;
            } else {
                self.not_n_counts[i] += 1;
            }
        }
    }

    fn reset(&mut self) {
        self.n_counts.clear();
        self.not_n_counts.clear();
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        self.check_threshold(self.error_threshold)
    }

    fn raises_warning(&mut self) -> bool {
        self.check_threshold(self.warn_threshold)
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, other: Box<dyn Any + Send>) {
        if let Ok(other) = other.downcast::<Self>() {
            let max_len = self.n_counts.len().max(other.n_counts.len());
            self.n_counts.resize(max_len, 0);
            self.not_n_counts.resize(max_len, 0);

            for (i, &v) in other.n_counts.iter().enumerate() {
                self.n_counts[i] += v;
            }
            for (i, &v) in other.not_n_counts.iter().enumerate() {
                self.not_n_counts[i] += v;
            }

            self.calculated = false;
        }
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate(&report.config);

        let max = self.percentages.iter().cloned().fold(0.0f64, f64::max);
        let svg = crate::graphs::line_graph::render_line_graph(
            &[&self.percentages],
            0.0,
            max.max(1.0),
            "Position in read (bp)",
            &["%N"],
            &self.x_labels,
            "N content across all bases",
            800,
            600,
        );
        report.add_image("per_base_n_content.png", &svg);

        let data = &mut report.data;
        data.push_str("#Base\tN-Count\n");
        for i in 0..self.x_labels.len() {
            data.push_str(&format!(
                "{}\t{}\n",
                self.x_labels[i],
                crate::format_double(self.percentages[i])
            ));
        }

        Ok(())
    }
}
