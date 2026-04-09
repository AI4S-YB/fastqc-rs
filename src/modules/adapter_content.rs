use crate::config::Config;
use crate::error::Result;
use crate::graphs::base_group;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;

struct Adapter {
    name: String,
    sequence: String,
    positions: Vec<u64>,
}

impl Adapter {
    fn new(name: String, sequence: String) -> Self {
        Self {
            name,
            sequence,
            positions: vec![0],
        }
    }

    fn expand_length_to(&mut self, new_length: usize) {
        if new_length <= self.positions.len() {
            return;
        }
        let last_val = *self.positions.last().unwrap_or(&0);
        self.positions.resize(new_length, last_val);
    }

    fn increment_count(&mut self, position: usize) {
        self.positions[position] += 1;
    }

    fn reset(&mut self) {
        self.positions.clear();
    }
}

pub struct AdapterContent {
    adapters: Vec<Adapter>,
    labels: Vec<String>,
    longest_sequence: usize,
    longest_adapter: usize,
    total_count: u64,
    calculated: bool,
    enrichments: Vec<Vec<f64>>,
    x_labels: Vec<String>,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl AdapterContent {
    pub fn new(config: &Config, module_config: &ModuleConfig) -> Self {
        let adapter_text = if let Some(ref path) = config.adapter_file {
            std::fs::read_to_string(path).unwrap_or_default()
        } else {
            include_str!("../data/adapter_list.txt").to_string()
        };

        let mut adapters = Vec::new();
        let mut labels = Vec::new();
        let mut longest_adapter = 0;

        for line in adapter_text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.splitn(2, '\t').collect();
            if parts.len() != 2 {
                // Try splitting on multiple tabs/spaces
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let name = parts[..parts.len() - 1].join(" ");
                    let seq = parts.last().unwrap().to_string();
                    if seq.len() > longest_adapter {
                        longest_adapter = seq.len();
                    }
                    labels.push(name.clone());
                    adapters.push(Adapter::new(name, seq));
                }
                continue;
            }
            let name = parts[0].trim().to_string();
            let seq = parts[1].trim().to_string();
            if seq.len() > longest_adapter {
                longest_adapter = seq.len();
            }
            labels.push(name.clone());
            adapters.push(Adapter::new(name, seq));
        }

        Self {
            adapters,
            labels,
            longest_sequence: 0,
            longest_adapter,
            total_count: 0,
            calculated: false,
            enrichments: Vec::new(),
            x_labels: Vec::new(),
            ignore: module_config.get_param("adapter", "ignore") > 0.0,
            warn_threshold: module_config.get_param("adapter", "warn"),
            error_threshold: module_config.get_param("adapter", "error"),
        }
    }

    fn calculate_enrichment(&mut self, config: &Config) {
        if self.calculated {
            return;
        }

        let max_length = self.adapters.iter().map(|a| a.positions.len()).max().unwrap_or(0);
        let groups = base_group::make_base_groups(max_length, config);

        self.x_labels = groups.iter().map(|g| g.to_string()).collect();
        self.enrichments = vec![vec![0.0; groups.len()]; self.adapters.len()];

        for (a, adapter) in self.adapters.iter().enumerate() {
            let positions = &adapter.positions;
            for (g, group) in groups.iter().enumerate() {
                let mut total = 0.0;
                for p in (group.lower_count() - 1)..group.upper_count().min(positions.len()) {
                    total += (positions[p] as f64 * 100.0) / self.total_count.max(1) as f64;
                }
                self.enrichments[a][g] =
                    total / (group.upper_count() - group.lower_count() + 1) as f64;
            }
        }

        self.calculated = true;
    }

    fn ensure_calculated(&mut self) {
        if !self.calculated {
            let config = crate::config::Config::default_config();
            self.calculate_enrichment(&config);
        }
    }
}

impl QcModule for AdapterContent {
    fn name(&self) -> &str {
        "Adapter Content"
    }

    fn description(&self) -> &str {
        "Searches for specific adapter sequences in a library"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.calculated = false;
        self.total_count += 1;

        let seq_len = seq.sequence.len();
        if seq_len > self.longest_sequence && seq_len > self.longest_adapter {
            self.longest_sequence = seq_len;
            let new_len = (seq_len - self.longest_adapter) + 1;
            for adapter in &mut self.adapters {
                adapter.expand_length_to(new_len);
            }
        }

        for adapter in &mut self.adapters {
            if let Some(index) = seq.sequence.find(&adapter.sequence) {
                let max_pos = if self.longest_sequence > self.longest_adapter {
                    self.longest_sequence - self.longest_adapter
                } else {
                    0
                };
                for i in index..=max_pos.min(adapter.positions.len() - 1) {
                    adapter.increment_count(i);
                }
            }
        }
    }

    fn reset(&mut self) {
        self.calculated = false;
        self.total_count = 0;
        self.longest_sequence = 0;
        for adapter in &mut self.adapters {
            adapter.reset();
        }
    }

    fn raises_error(&mut self) -> bool {
        self.ensure_calculated();
        for row in &self.enrichments {
            for &val in row {
                if val > self.error_threshold {
                    return true;
                }
            }
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        if self.longest_adapter > self.longest_sequence {
            return true;
        }
        self.ensure_calculated();
        for row in &self.enrichments {
            for &val in row {
                if val > self.warn_threshold {
                    return true;
                }
            }
        }
        false
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        if self.longest_adapter > self.longest_sequence {
            report.html_body.push_str(&format!(
                "<p>Can't analyse adapters as read length is too short ({} vs {})</p>",
                self.longest_adapter, self.longest_sequence
            ));
            return Ok(());
        }

        self.calculate_enrichment(&report.config);

        let refs: Vec<&[f64]> = self.enrichments.iter().map(|v| v.as_slice()).collect();
        let label_refs: Vec<&str> = self.labels.iter().map(|s| s.as_str()).collect();
        let width = 800.max(self.x_labels.len() * 15);

        let svg = crate::graphs::line_graph::render_line_graph(
            &refs,
            0.0,
            100.0,
            "Position in read (bp)",
            &label_refs,
            &self.x_labels,
            "% Adapter",
            width,
            600,
        );
        report.add_image("adapter_content.png", &svg);

        let data = &mut report.data;
        // Header
        data.push('#');
        data.push_str("Position");
        for label in &self.labels {
            data.push('\t');
            data.push_str(label);
        }
        data.push('\n');

        // Data rows
        for g in 0..self.x_labels.len() {
            data.push_str(&self.x_labels[g]);
            for a in 0..self.adapters.len() {
                data.push('\t');
                data.push_str(&crate::format_double(self.enrichments[a][g]));
            }
            data.push('\n');
        }

        Ok(())
    }
}
