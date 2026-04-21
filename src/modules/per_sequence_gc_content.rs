use std::any::Any;

use crate::error::Result;
use crate::modules::gc_model::GcModel;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;
use crate::stats::normal_distribution::NormalDistribution;

pub struct PerSequenceGCContent {
    gc_distribution: [f64; 101],
    theoretical_distribution: [f64; 101],
    deviation_percent: f64,
    calculated: bool,
    cached_models: Vec<Option<GcModel>>,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl PerSequenceGCContent {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            gc_distribution: [0.0; 101],
            theoretical_distribution: [0.0; 101],
            deviation_percent: 0.0,
            calculated: false,
            cached_models: vec![None; 200],
            ignore: module_config.get_param("gc_sequence", "ignore") > 0.0,
            warn_threshold: module_config.get_param("gc_sequence", "warn"),
            error_threshold: module_config.get_param("gc_sequence", "error"),
        }
    }

    fn truncate_sequence(seq: &str) -> &str {
        let len = seq.len();
        if len > 1000 {
            let truncated = (len / 1000) * 1000;
            &seq[..truncated]
        } else if len > 100 {
            let truncated = (len / 100) * 100;
            &seq[..truncated]
        } else {
            seq
        }
    }

    fn calculate(&mut self) {
        if self.calculated {
            return;
        }

        let mut max = 0.0f64;
        let mut total_count = 0.0f64;
        let mut first_mode = 0usize;
        let mut mode_count = 0.0f64;

        for i in 0..101 {
            total_count += self.gc_distribution[i];
            if self.gc_distribution[i] > mode_count {
                mode_count = self.gc_distribution[i];
                first_mode = i;
            }
            if self.gc_distribution[i] > max {
                max = self.gc_distribution[i];
            }
        }

        // Average mode over adjacent positions within 10% of modal count
        let mut mode = 0.0f64;
        let mut mode_duplicates = 0;
        let threshold =
            self.gc_distribution[first_mode] - (self.gc_distribution[first_mode] / 10.0);

        let mut fell_off_top = true;
        for i in first_mode..101 {
            if self.gc_distribution[i] > threshold {
                mode += i as f64;
                mode_duplicates += 1;
            } else {
                fell_off_top = false;
                break;
            }
        }

        let mut fell_off_bottom = true;
        if first_mode > 0 {
            for i in (0..first_mode).rev() {
                if self.gc_distribution[i] > threshold {
                    mode += i as f64;
                    mode_duplicates += 1;
                } else {
                    fell_off_bottom = false;
                    break;
                }
            }
        }

        if fell_off_bottom || fell_off_top {
            mode = first_mode as f64;
        } else {
            mode /= mode_duplicates as f64;
        }

        // Calculate standard deviation
        let mut stdev = 0.0f64;
        for i in 0..101 {
            stdev += (i as f64 - mode).powi(2) * self.gc_distribution[i];
        }
        if total_count > 1.0 {
            stdev /= total_count - 1.0;
        }
        stdev = stdev.sqrt();

        // Fit normal distribution
        let nd = NormalDistribution::new(mode, stdev);
        self.deviation_percent = 0.0;

        for i in 0..101 {
            let probability = nd.pdf(i as f64);
            self.theoretical_distribution[i] = probability * total_count;
            if self.theoretical_distribution[i] > max {
                max = self.theoretical_distribution[i];
            }
            self.deviation_percent +=
                (self.theoretical_distribution[i] - self.gc_distribution[i]).abs();
        }

        if total_count > 0.0 {
            self.deviation_percent /= total_count;
            self.deviation_percent *= 100.0;
        }

        self.calculated = true;
    }
}

impl QcModule for PerSequenceGCContent {
    fn name(&self) -> &str {
        "Per sequence GC content"
    }

    fn description(&self) -> &str {
        "Shows the distribution of GC contents for whole sequences"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.calculated = false;
        let truncated = Self::truncate_sequence(&seq.sequence);
        let len = truncated.len();
        if len == 0 {
            return;
        }

        let gc_count = truncated
            .bytes()
            .filter(|&c| c == b'G' || c == b'C')
            .count();

        // Ensure cached models capacity
        if len >= self.cached_models.len() {
            self.cached_models.resize_with(len + 1, || None);
        }

        if self.cached_models[len].is_none() {
            self.cached_models[len] = Some(GcModel::new(len));
        }

        let model = self.cached_models[len].as_ref().unwrap();
        let values = model.get_model_values(gc_count);

        for &(percentage, increment) in values {
            self.gc_distribution[percentage] += increment;
        }
    }

    fn reset(&mut self) {
        self.gc_distribution = [0.0; 101];
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        self.calculate();
        self.deviation_percent > self.error_threshold
    }

    fn raises_warning(&mut self) -> bool {
        self.calculate();
        self.deviation_percent > self.warn_threshold
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, other: Box<dyn Any + Send>) {
        if let Ok(other) = other.downcast::<Self>() {
            for i in 0..101 {
                self.gc_distribution[i] += other.gc_distribution[i];
            }
            self.calculated = false;
        }
    }

    fn module_id(&self) -> &str {
        "per_sequence_gc_content"
    }

    fn json_thresholds(&self) -> Option<serde_json::Value> {
        Some(serde_json::json!({
            "warn": self.warn_threshold,
            "error": self.error_threshold,
        }))
    }

    fn json_data(&mut self, _config: &crate::config::Config) -> serde_json::Value {
        self.calculate();
        let distribution: Vec<serde_json::Value> = (0..101)
            .map(|i| {
                serde_json::json!({
                    "gc_percent": i,
                    "observed_count": crate::modules::json_num(self.gc_distribution[i]),
                    "theoretical_count": crate::modules::json_num(self.theoretical_distribution[i]),
                })
            })
            .collect();
        serde_json::json!({
            "deviation_percent": crate::modules::json_num(self.deviation_percent),
            "distribution": distribution,
        })
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate();

        let max = self
            .gc_distribution
            .iter()
            .chain(self.theoretical_distribution.iter())
            .cloned()
            .fold(0.0f64, f64::max);

        let x_labels: Vec<String> = (0..101).map(|i| format!("{}", i)).collect();

        let svg = crate::graphs::line_graph::render_line_graph(
            &[&self.gc_distribution, &self.theoretical_distribution],
            0.0,
            max,
            "Mean GC content (%)",
            &["GC count per read", "Theoretical Distribution"],
            &x_labels,
            "GC distribution over all sequences",
            800,
            600,
        );
        report.add_image("per_sequence_gc_content.png", &svg);

        let data = &mut report.data;
        data.push_str("#GC Content\tCount\n");
        for i in 0..101 {
            data.push_str(&format!(
                "{}\t{}\n",
                i,
                crate::format_double(self.gc_distribution[i])
            ));
        }

        Ok(())
    }
}
