use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use crate::config::Config;
use crate::error::Result;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::contaminant::ContaminantFinder;
use crate::sequence::Sequence;

/// Shared data between OverRepresentedSeqs and DuplicationLevel.
pub struct SharedDuplicationData {
    pub sequences: HashMap<String, u64>,
    pub count: u64,
    pub count_at_unique_limit: u64,
    pub frozen: bool,
    pub unique_sequence_count: usize,
}

impl SharedDuplicationData {
    pub fn new() -> Self {
        Self {
            sequences: HashMap::new(),
            count: 0,
            count_at_unique_limit: 0,
            frozen: false,
            unique_sequence_count: 0,
        }
    }
}

const OBSERVATION_CUTOFF: usize = 100000;

pub struct OverRepresentedSeqs {
    shared: Arc<Mutex<SharedDuplicationData>>,
    dup_length: usize,
    calculated: bool,
    overrepresented: Vec<OverrepresentedSeq>,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
    contaminant_finder: ContaminantFinder,
}

struct OverrepresentedSeq {
    seq: String,
    count: u64,
    percentage: f64,
    contaminant_hit: String,
}

impl OverRepresentedSeqs {
    pub fn new(config: &Config, module_config: &ModuleConfig) -> Self {
        Self {
            shared: Arc::new(Mutex::new(SharedDuplicationData::new())),
            dup_length: config.dup_length,
            calculated: false,
            overrepresented: Vec::new(),
            ignore: module_config.get_param("overrepresented", "ignore") > 0.0,
            warn_threshold: module_config.get_param("overrepresented", "warn"),
            error_threshold: module_config.get_param("overrepresented", "error"),
            contaminant_finder: ContaminantFinder::new(config),
        }
    }

    pub fn shared_data(&self) -> Arc<Mutex<SharedDuplicationData>> {
        self.shared.clone()
    }

    fn calculate_overrepresented(&mut self) {
        if self.calculated {
            return;
        }

        let data = self.shared.lock().unwrap();
        let count = data.count;
        let mut keepers: Vec<OverrepresentedSeq> = Vec::new();

        for (seq, &seq_count) in &data.sequences {
            let percentage = seq_count as f64 / count as f64 * 100.0;
            if percentage > self.warn_threshold {
                let hit = self.contaminant_finder.find_hit(seq);
                keepers.push(OverrepresentedSeq {
                    seq: seq.clone(),
                    count: seq_count,
                    percentage,
                    contaminant_hit: hit.unwrap_or_else(|| "No Hit".to_string()),
                });
            }
        }

        keepers.sort_by(|a, b| b.count.cmp(&a.count));
        drop(data);

        self.overrepresented = keepers;
        self.calculated = true;
    }
}

impl QcModule for OverRepresentedSeqs {
    fn name(&self) -> &str {
        "Overrepresented sequences"
    }

    fn description(&self) -> &str {
        "Identifies sequences which are overrepresented in the set"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.calculated = false;

        let mut data = self.shared.lock().unwrap();
        data.count += 1;

        // Truncate sequence for comparison
        let s = &seq.sequence;
        let truncated = if self.dup_length > 0 && s.len() > self.dup_length {
            &s[..self.dup_length]
        } else if s.len() > 50 {
            &s[..50]
        } else {
            s.as_str()
        };

        if let Some(count) = data.sequences.get_mut(truncated) {
            *count += 1;
            if !data.frozen {
                data.count_at_unique_limit = data.count;
            }
        } else if !data.frozen {
            data.sequences.insert(truncated.to_string(), 1);
            data.unique_sequence_count += 1;
            data.count_at_unique_limit = data.count;
            if data.unique_sequence_count == OBSERVATION_CUTOFF {
                data.frozen = true;
            }
        }
    }

    fn reset(&mut self) {
        let mut data = self.shared.lock().unwrap();
        data.count = 0;
        data.sequences.clear();
        data.frozen = false;
        data.unique_sequence_count = 0;
        data.count_at_unique_limit = 0;
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        self.calculate_overrepresented();
        if let Some(first) = self.overrepresented.first() {
            first.percentage > self.error_threshold
        } else {
            false
        }
    }

    fn raises_warning(&mut self) -> bool {
        self.calculate_overrepresented();
        !self.overrepresented.is_empty()
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate_overrepresented();

        if self.overrepresented.is_empty() {
            report.html_body.push_str("<p>No overrepresented sequences</p>");
        } else {
            report.html_body.push_str("<table><thead><tr><th>Sequence</th><th>Count</th><th>Percentage</th><th>Possible Source</th></tr></thead><tbody>");
            for seq in &self.overrepresented {
                let rounded_pct = (seq.percentage * 100.0).round() / 100.0;
                report.html_body.push_str(&format!(
                    "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                    seq.seq, seq.count, rounded_pct, seq.contaminant_hit
                ));
            }
            report.html_body.push_str("</tbody></table>");
        }

        let data = &mut report.data;
        if self.overrepresented.is_empty() {
            // No data section needed beyond the header
        } else {
            data.push_str("#Sequence\tCount\tPercentage\tPossible Source\n");
            for seq in &self.overrepresented {
                let rounded_pct = (seq.percentage * 100.0).round() / 100.0;
                data.push_str(&format!(
                    "{}\t{}\t{}\t{}\n",
                    seq.seq, seq.count, crate::format_double(rounded_pct), seq.contaminant_hit
                ));
            }
        }

        Ok(())
    }
}
