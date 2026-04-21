use std::any::Any;
use std::sync::{Arc, Mutex};

use ahash::AHashMap;

use crate::error::Result;
use crate::modules::overrepresented_seqs::SharedDuplicationData;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;

pub struct DuplicationLevel {
    shared: Arc<Mutex<SharedDuplicationData>>,
    total_percentages: Option<[f64; 16]>,
    percent_different_seqs: f64,
    labels: [&'static str; 16],
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl DuplicationLevel {
    pub fn new(shared: Arc<Mutex<SharedDuplicationData>>, module_config: &ModuleConfig) -> Self {
        Self {
            shared,
            total_percentages: None,
            percent_different_seqs: 0.0,
            labels: [
                "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k",
                ">5k", ">10k",
            ],
            ignore: module_config.get_param("duplication", "ignore") > 0.0,
            warn_threshold: module_config.get_param("duplication", "warn"),
            error_threshold: module_config.get_param("duplication", "error"),
        }
    }

    pub fn calculate_levels(&mut self) {
        if self.total_percentages.is_some() {
            return;
        }

        let data = self.shared.lock().unwrap();
        let mut total_percentages = [0.0f64; 16];

        // Collate counts: how many sequences have each duplication level
        let mut collated_counts: AHashMap<u64, u64> = AHashMap::new();
        for &count in data.sequences.values() {
            *collated_counts.entry(count).or_insert(0) += 1;
        }

        // Apply Bayesian correction
        let mut corrected_counts: AHashMap<u64, f64> = AHashMap::new();
        for (&dup_level, &num_obs) in &collated_counts {
            let corrected =
                get_corrected_count(data.count_at_unique_limit, data.count, dup_level, num_obs);
            corrected_counts.insert(dup_level, corrected);
        }

        // Calculate totals and bin into slots
        let mut dedup_total = 0.0f64;
        let mut raw_total = 0.0f64;

        for (&dup_level, &count) in &corrected_counts {
            dedup_total += count;
            raw_total += count * dup_level as f64;

            let temp_dup_slot = dup_level as i64 - 1;
            let dup_slot = if temp_dup_slot > 9999 || temp_dup_slot < 0 {
                15
            } else if temp_dup_slot > 4999 {
                14
            } else if temp_dup_slot > 999 {
                13
            } else if temp_dup_slot > 499 {
                12
            } else if temp_dup_slot > 99 {
                11
            } else if temp_dup_slot > 49 {
                10
            } else if temp_dup_slot > 9 {
                9
            } else {
                temp_dup_slot as usize
            };

            total_percentages[dup_slot] += count * dup_level as f64;
        }

        // Normalize to percentages
        if raw_total > 0.0 {
            for tp in total_percentages.iter_mut() {
                *tp /= raw_total;
                *tp *= 100.0;
            }
            self.percent_different_seqs = (dedup_total / raw_total) * 100.0;
        } else {
            self.percent_different_seqs = 100.0;
        }

        self.total_percentages = Some(total_percentages);
    }
}

fn get_corrected_count(
    count_at_limit: u64,
    total_count: u64,
    duplication_level: u64,
    number_of_observations: u64,
) -> f64 {
    // Early bail-out
    if count_at_limit == total_count {
        return number_of_observations as f64;
    }
    if total_count - number_of_observations < count_at_limit {
        return number_of_observations as f64;
    }

    let mut p_not_seeing_at_limit = 1.0f64;
    let limit_of_caring =
        1.0 - (number_of_observations as f64 / (number_of_observations as f64 + 0.01));

    for i in 0..count_at_limit {
        p_not_seeing_at_limit *=
            ((total_count - i) as f64 - duplication_level as f64) / (total_count - i) as f64;

        if p_not_seeing_at_limit < limit_of_caring {
            p_not_seeing_at_limit = 0.0;
            break;
        }
    }

    let p_seeing_at_limit = 1.0 - p_not_seeing_at_limit;
    number_of_observations as f64 / p_seeing_at_limit
}

impl QcModule for DuplicationLevel {
    fn name(&self) -> &str {
        "Sequence Duplication Levels"
    }

    fn description(&self) -> &str {
        "Plots the number of sequences which are duplicated to different levels"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        self.ignore
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, _seq: &Sequence) {
        // No-op: uses data from OverRepresentedSeqs
    }

    fn reset(&mut self) {
        self.total_percentages = None;
    }

    fn raises_error(&mut self) -> bool {
        self.calculate_levels();
        self.percent_different_seqs < self.error_threshold
    }

    fn raises_warning(&mut self) -> bool {
        self.calculate_levels();
        self.percent_different_seqs < self.warn_threshold
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, _other: Box<dyn Any + Send>) {
        // No-op: DuplicationLevel reads from shared OverRepresentedSeqs data at report time.
    }

    fn module_id(&self) -> &str {
        "sequence_duplication_levels"
    }

    fn json_thresholds(&self) -> Option<serde_json::Value> {
        Some(serde_json::json!({
            "warn": self.warn_threshold,
            "error": self.error_threshold,
        }))
    }

    fn json_data(&mut self, _config: &crate::config::Config) -> serde_json::Value {
        self.calculate_levels();
        let percentages = self.total_percentages.unwrap_or([0.0; 16]);
        let distribution: Vec<serde_json::Value> = (0..16)
            .map(|i| {
                let label = if i == 15 {
                    format!("{}+", self.labels[i])
                } else {
                    self.labels[i].to_string()
                };
                serde_json::json!({
                    "label": label,
                    "percentage_of_total": crate::modules::json_num(percentages[i]),
                })
            })
            .collect();
        serde_json::json!({
            "deduplicated_percentage": crate::modules::json_num(self.percent_different_seqs),
            "distribution": distribution,
        })
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate_levels();
        let percentages = self.total_percentages.as_ref().unwrap();

        let pcts_vec: Vec<f64> = percentages.to_vec();
        let x_labels: Vec<String> = self.labels.iter().map(|l| l.to_string()).collect();

        let svg = crate::graphs::line_graph::render_line_graph(
            &[&pcts_vec],
            0.0,
            100.0,
            "Sequence Duplication Level",
            &["% Total sequences"],
            &x_labels,
            &format!(
                "Percent of seqs remaining if deduplicated {:.2}%",
                self.percent_different_seqs
            ),
            800,
            600,
        );
        report.add_image("duplication_levels.png", &svg);

        let data = &mut report.data;
        data.push_str(&format!(
            "#Total Deduplicated Percentage\t{}\n",
            crate::format_double(self.percent_different_seqs)
        ));
        data.push_str("#Duplication Level\tPercentage of total\n");
        for i in 0..16 {
            data.push_str(self.labels[i]);
            if i == 15 {
                data.push('+');
            }
            data.push('\t');
            data.push_str(&crate::format_double(percentages[i]));
            data.push('\n');
        }

        Ok(())
    }
}
