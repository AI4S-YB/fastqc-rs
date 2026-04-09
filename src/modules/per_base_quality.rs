use std::any::Any;

use crate::error::Result;
use crate::graphs::base_group;
use crate::graphs::quality_count::QualityCount;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::phred::PhredEncoding;
use crate::sequence::Sequence;

pub struct PerBaseQualityScores {
    quality_counts: Vec<QualityCount>,
    means: Option<Vec<f64>>,
    medians: Option<Vec<f64>>,
    lower_quartile: Option<Vec<f64>>,
    upper_quartile: Option<Vec<f64>>,
    lowest: Option<Vec<f64>>,
    highest: Option<Vec<f64>>,
    x_labels: Option<Vec<String>>,
    ignore: bool,
    warn_lower: f64,
    warn_median: f64,
    error_lower: f64,
    error_median: f64,
}

impl PerBaseQualityScores {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            quality_counts: Vec::new(),
            means: None,
            medians: None,
            lower_quartile: None,
            upper_quartile: None,
            lowest: None,
            highest: None,
            x_labels: None,
            ignore: module_config.get_param("quality_base", "ignore") > 0.0,
            warn_lower: module_config.get_param("quality_base_lower", "warn"),
            warn_median: module_config.get_param("quality_base_median", "warn"),
            error_lower: module_config.get_param("quality_base_lower", "error"),
            error_median: module_config.get_param("quality_base_median", "error"),
        }
    }

    fn calculate(&mut self, config: &crate::config::Config) {
        if self.means.is_some() {
            return;
        }

        let (min_char, _) = self.calculate_offsets();
        let encoding = PhredEncoding::from_lowest_char(min_char);
        let offset = encoding.offset();

        let groups = base_group::make_base_groups(self.quality_counts.len(), config);

        let mut means = Vec::with_capacity(groups.len());
        let mut medians = Vec::with_capacity(groups.len());
        let mut lower_quartile = Vec::with_capacity(groups.len());
        let mut upper_quartile = Vec::with_capacity(groups.len());
        let mut lowest = Vec::with_capacity(groups.len());
        let mut highest = Vec::with_capacity(groups.len());
        let mut x_labels = Vec::with_capacity(groups.len());

        for group in &groups {
            x_labels.push(group.to_string());
            let min_base = group.lower_count();
            let max_base = group.upper_count();
            lowest.push(self.get_percentile(min_base, max_base, offset, 10));
            highest.push(self.get_percentile(min_base, max_base, offset, 90));
            means.push(self.get_mean(min_base, max_base, offset));
            medians.push(self.get_percentile(min_base, max_base, offset, 50));
            lower_quartile.push(self.get_percentile(min_base, max_base, offset, 25));
            upper_quartile.push(self.get_percentile(min_base, max_base, offset, 75));
        }

        self.means = Some(means);
        self.medians = Some(medians);
        self.lower_quartile = Some(lower_quartile);
        self.upper_quartile = Some(upper_quartile);
        self.lowest = Some(lowest);
        self.highest = Some(highest);
        self.x_labels = Some(x_labels);
    }

    fn calculate_offsets(&self) -> (u8, u8) {
        let mut min_char: u8 = 255;
        let mut max_char: u8 = 0;
        for qc in &self.quality_counts {
            let mn = qc.min_char();
            let mx = qc.max_char();
            if mn < min_char {
                min_char = mn;
            }
            if mx > max_char {
                max_char = mx;
            }
        }
        (min_char, max_char)
    }

    fn get_percentile(&self, min_bp: usize, max_bp: usize, offset: usize, percentile: u32) -> f64 {
        let mut count = 0;
        let mut total = 0.0;
        for i in (min_bp - 1)..max_bp.min(self.quality_counts.len()) {
            if self.quality_counts[i].total_count() > 100 {
                count += 1;
                total += self.quality_counts[i].percentile(offset, percentile);
            }
        }
        if count > 0 {
            total / count as f64
        } else {
            f64::NAN
        }
    }

    fn get_mean(&self, min_bp: usize, max_bp: usize, offset: usize) -> f64 {
        let mut count = 0;
        let mut total = 0.0;
        for i in (min_bp - 1)..max_bp.min(self.quality_counts.len()) {
            if self.quality_counts[i].total_count() > 0 {
                count += 1;
                total += self.quality_counts[i].mean(offset);
            }
        }
        if count > 0 {
            total / count as f64
        } else {
            0.0
        }
    }

    fn ensure_calculated(&mut self) {
        if self.means.is_none() {
            let config = crate::config::Config::default_config();
            self.calculate(&config);
        }
    }
}

impl QcModule for PerBaseQualityScores {
    fn name(&self) -> &str {
        "Per base sequence quality"
    }

    fn description(&self) -> &str {
        "Shows the Quality scores of all bases at a given position in a sequencing run"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore || self.quality_counts.is_empty()
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.means = None; // Invalidate cache
        let qual = seq.quality.as_bytes();
        if self.quality_counts.len() < qual.len() {
            self.quality_counts
                .resize_with(qual.len(), QualityCount::new);
        }
        for (i, &c) in qual.iter().enumerate() {
            self.quality_counts[i].add_value(c);
        }
    }

    fn reset(&mut self) {
        self.quality_counts.clear();
        self.means = None;
    }

    fn raises_error(&mut self) -> bool {
        self.ensure_calculated();
        let lq = self.lower_quartile.as_ref().unwrap();
        let med = self.medians.as_ref().unwrap();
        for i in 0..lq.len() {
            if lq[i].is_nan() {
                continue;
            }
            if lq[i] < self.error_lower || med[i] < self.error_median {
                return true;
            }
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        self.ensure_calculated();
        let lq = self.lower_quartile.as_ref().unwrap();
        let med = self.medians.as_ref().unwrap();
        for i in 0..lq.len() {
            if lq[i].is_nan() {
                continue;
            }
            if lq[i] < self.warn_lower || med[i] < self.warn_median {
                return true;
            }
        }
        false
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, other: Box<dyn Any + Send>) {
        if let Ok(other) = other.downcast::<Self>() {
            let max_len = self.quality_counts.len().max(other.quality_counts.len());
            self.quality_counts.resize_with(max_len, QualityCount::new);

            for (i, other_qc) in other.quality_counts.iter().enumerate() {
                self.quality_counts[i].merge_from(other_qc);
            }

            // Invalidate all caches
            self.means = None;
            self.medians = None;
            self.lower_quartile = None;
            self.upper_quartile = None;
            self.lowest = None;
            self.highest = None;
            self.x_labels = None;
        }
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate(&report.config);

        let means = self.means.as_ref().unwrap();
        let medians = self.medians.as_ref().unwrap();
        let lq = self.lower_quartile.as_ref().unwrap();
        let uq = self.upper_quartile.as_ref().unwrap();
        let lowest = self.lowest.as_ref().unwrap();
        let highest = self.highest.as_ref().unwrap();
        let x_labels = self.x_labels.as_ref().unwrap();

        // Generate SVG
        let (min_char, max_char) = self.calculate_offsets();
        let encoding = PhredEncoding::from_lowest_char(min_char);
        let high = {
            let h = max_char as i32 - encoding.offset() as i32;
            if h < 35 {
                35.0
            } else {
                h as f64
            }
        };
        let width = 800.max(means.len() * 15);
        let svg = crate::graphs::box_plot::render_quality_box_plot(
            means,
            medians,
            lowest,
            highest,
            lq,
            uq,
            0.0,
            high,
            x_labels,
            &format!("Quality scores across all bases ({} encoding)", encoding),
            width,
            600,
        );
        report.add_image("per_base_quality.png", &svg);

        // Write data
        let data = &mut report.data;
        data.push_str("#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n");
        for i in 0..means.len() {
            data.push_str(&x_labels[i]);
            data.push('\t');
            data.push_str(&crate::format_double(means[i]));
            data.push('\t');
            data.push_str(&crate::format_double(medians[i]));
            data.push('\t');
            data.push_str(&crate::format_double(lq[i]));
            data.push('\t');
            data.push_str(&crate::format_double(uq[i]));
            data.push('\t');
            data.push_str(&crate::format_double(lowest[i]));
            data.push('\t');
            data.push_str(&crate::format_double(highest[i]));
            data.push('\n');
        }

        Ok(())
    }
}
