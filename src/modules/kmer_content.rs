use std::any::Any;

use ahash::AHashMap;

use crate::config::Config;
use crate::error::Result;
use crate::graphs::base_group;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;

struct Kmer {
    sequence: String,
    count: u64,
    positions: Vec<u64>,
    obs_exp_positions: Vec<f32>,
    lowest_p_value: f32,
}

impl Kmer {
    fn new(sequence: String, position: usize, seq_length: usize) -> Self {
        let mut positions = vec![0u64; seq_length];
        positions[position] = 1;
        Self {
            sequence,
            count: 1,
            positions,
            obs_exp_positions: Vec::new(),
            lowest_p_value: 0.0,
        }
    }

    fn increment_count(&mut self, position: usize) {
        self.count += 1;
        if position >= self.positions.len() {
            self.positions.resize(position + 1, 0);
        }
        self.positions[position] += 1;
    }

    fn max_obs_exp(&self) -> f32 {
        self.obs_exp_positions
            .iter()
            .cloned()
            .fold(0.0f32, f32::max)
    }

    fn max_position(&self) -> usize {
        let mut max = 0.0f32;
        let mut pos = 0;
        for (i, &v) in self.obs_exp_positions.iter().enumerate() {
            if v > max {
                max = v;
                pos = i + 1;
            }
        }
        if pos == 0 {
            1
        } else {
            pos
        }
    }
}

pub struct KmerContent {
    kmers: AHashMap<String, Kmer>,
    total_kmer_counts: Vec<Vec<u64>>,
    longest_sequence: usize,
    skip_count: u64,
    kmer_size: usize,
    calculated: bool,
    enriched_kmers: Vec<usize>, // indices into sorted_kmers
    sorted_kmer_data: Vec<KmerResult>,
    enrichments: Vec<Vec<f64>>,
    x_categories: Vec<String>,
    x_labels: Vec<String>,
    max_graph_value: f64,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

struct KmerResult {
    sequence: String,
    count: u64,
    p_value: f32,
    max_obs_exp: f32,
    max_position: String,
}

impl KmerContent {
    pub fn new(config: &Config, module_config: &ModuleConfig) -> Self {
        let kmer_size = config.kmer_size.unwrap_or(7);
        Self {
            kmers: AHashMap::new(),
            total_kmer_counts: Vec::new(),
            longest_sequence: 0,
            skip_count: 0,
            kmer_size,
            calculated: false,
            enriched_kmers: Vec::new(),
            sorted_kmer_data: Vec::new(),
            enrichments: Vec::new(),
            x_categories: Vec::new(),
            x_labels: Vec::new(),
            max_graph_value: 0.0,
            ignore: module_config.get_param("kmer", "ignore") > 0.0,
            warn_threshold: module_config.get_param("kmer", "warn"),
            error_threshold: module_config.get_param("kmer", "error"),
        }
    }

    fn calculate_enrichment(&mut self, config: &Config) {
        if self.calculated {
            return;
        }

        let groups = base_group::make_base_groups(
            if self.longest_sequence >= self.kmer_size {
                (self.longest_sequence - self.kmer_size) + 1
            } else {
                0
            },
            config,
        );

        self.x_categories = groups.iter().map(|g| g.to_string()).collect();

        let mut enriched: Vec<(f32, f32, String, u64, Vec<f32>)> = Vec::new();

        let kmer_len_idx = self.kmer_size - 1;

        for kmer in self.kmers.values() {
            // Total kmer count of this length across all positions
            let total_kmer_count: u64 = self
                .total_kmer_counts
                .iter()
                .filter(|row| row.len() > kmer_len_idx)
                .map(|row| row[kmer_len_idx])
                .sum();

            if total_kmer_count == 0 {
                continue;
            }

            let expected_proportion = kmer.count as f32 / total_kmer_count as f32;

            let mut obs_exp_positions = vec![0.0f32; groups.len()];
            let mut has_significant = false;
            let mut lowest_p = 1.0f32;

            for (g, group) in groups.iter().enumerate() {
                let mut total_group_count = 0u64;
                let mut total_group_hits = 0u64;

                for p in (group.lower_count() - 1)..group.upper_count().min(kmer.positions.len()) {
                    if p < self.total_kmer_counts.len()
                        && self.total_kmer_counts[p].len() > kmer_len_idx
                    {
                        total_group_count += self.total_kmer_counts[p][kmer_len_idx];
                    }
                    total_group_hits += kmer.positions[p];
                }

                let predicted = expected_proportion * total_group_count as f32;
                obs_exp_positions[g] = if predicted > 0.0 {
                    total_group_hits as f32 / predicted
                } else {
                    0.0
                };

                // Simplified p-value estimation (skip full binomial for now)
                if total_group_hits as f32 > predicted && obs_exp_positions[g] > 5.0 {
                    let p_value = approximate_binomial_p(
                        total_group_count as usize,
                        expected_proportion as f64,
                        total_group_hits as usize,
                    ) * 4.0f64.powi(kmer.sequence.len() as i32);

                    if (p_value as f32) < 0.01 {
                        has_significant = true;
                        if (p_value as f32) < lowest_p {
                            lowest_p = p_value as f32;
                        }
                    }
                }
            }

            if has_significant {
                enriched.push((
                    lowest_p,
                    obs_exp_positions.iter().cloned().fold(0.0f32, f32::max),
                    kmer.sequence.clone(),
                    kmer.count,
                    obs_exp_positions,
                ));
            }
        }

        // Sort by max obs/exp descending
        enriched.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        // Take top 20
        enriched.truncate(20);

        // Build result data
        self.sorted_kmer_data.clear();
        self.enrichments.clear();
        self.x_labels.clear();
        self.max_graph_value = 0.0;

        for (i, (p_val, max_oe, seq, count, obs_exp)) in enriched.iter().enumerate() {
            let max_pos_idx = obs_exp
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, _)| i)
                .unwrap_or(0);

            self.sorted_kmer_data.push(KmerResult {
                sequence: seq.clone(),
                count: count * 50, // Correct for 2% sampling
                p_value: *p_val,
                max_obs_exp: *max_oe,
                max_position: if max_pos_idx < self.x_categories.len() {
                    self.x_categories[max_pos_idx].clone()
                } else {
                    "1".to_string()
                },
            });

            // Top 6 for graph
            if i < 6 {
                let enrichment: Vec<f64> = obs_exp.iter().map(|&v| v as f64).collect();
                for &v in &enrichment {
                    if v > self.max_graph_value {
                        self.max_graph_value = v;
                    }
                }
                self.x_labels.push(seq.clone());
                self.enrichments.push(enrichment);
            }
        }

        self.kmers.clear();
        self.calculated = true;
    }

    fn ensure_calculated(&mut self) {
        if !self.calculated {
            let config = crate::config::Config::default_config();
            self.calculate_enrichment(&config);
        }
    }
}

/// Simplified binomial p-value approximation using normal approximation.
fn approximate_binomial_p(n: usize, p: f64, k: usize) -> f64 {
    if n == 0 || p <= 0.0 || p >= 1.0 {
        return 1.0;
    }
    let mean = n as f64 * p;
    let std = (n as f64 * p * (1.0 - p)).sqrt();
    if std < f64::EPSILON {
        return 1.0;
    }
    let z = (k as f64 - mean) / std;
    // Approximate upper tail probability using error function approximation
    0.5 * erfc(z / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    // Abramowitz and Stegun approximation
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = poly * (-x * x).exp();
    if x >= 0.0 {
        result
    } else {
        2.0 - result
    }
}

impl QcModule for KmerContent {
    fn name(&self) -> &str {
        "Kmer Content"
    }

    fn description(&self) -> &str {
        "Identifies short sequences which have uneven representation"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.calculated = false;
        self.skip_count += 1;
        if self.skip_count % 50 != 0 {
            return;
        }

        let s = if seq.sequence.len() > 500 {
            &seq.sequence[..500]
        } else {
            &seq.sequence
        };

        if s.len() > self.longest_sequence {
            self.longest_sequence = s.len();
        }

        let kmer_size = self.kmer_size;
        if s.len() < kmer_size {
            return;
        }

        for i in 0..=(s.len() - kmer_size) {
            let kmer = &s[i..i + kmer_size];

            // Expand total kmer counts
            if i >= self.total_kmer_counts.len() {
                self.total_kmer_counts
                    .resize_with(i + 1, || vec![0u64; kmer_size]);
            }

            // Count total regardless of N content
            if !kmer.contains('N') {
                self.total_kmer_counts[i][kmer_size - 1] += 1;
            }

            // Skip kmers with N
            if kmer.contains('N') {
                continue;
            }

            if let Some(k) = self.kmers.get_mut(kmer) {
                k.increment_count(i);
            } else {
                self.kmers.insert(
                    kmer.to_string(),
                    Kmer::new(kmer.to_string(), i, (s.len() - kmer_size) + 1),
                );
            }
        }
    }

    fn reset(&mut self) {
        self.calculated = false;
        self.total_kmer_counts.clear();
        self.longest_sequence = 0;
        self.skip_count = 0;
        self.kmers.clear();
    }

    fn raises_error(&mut self) -> bool {
        self.ensure_calculated();
        if let Some(first) = self.sorted_kmer_data.first() {
            if first.p_value > 0.0 {
                return (0.0 - (first.p_value as f64).log10()) > self.error_threshold;
            }
        }
        false
    }

    fn raises_warning(&mut self) -> bool {
        self.ensure_calculated();
        if let Some(first) = self.sorted_kmer_data.first() {
            if first.p_value > 0.0 {
                return (0.0 - (first.p_value as f64).log10()) > self.warn_threshold;
            }
        }
        false
    }

    fn into_any(self: Box<Self>) -> Box<dyn Any + Send> {
        self
    }

    fn merge(&mut self, other: Box<dyn Any + Send>) {
        if let Ok(other) = other.downcast::<Self>() {
            // Merge kmers
            for (seq, other_kmer) in other.kmers {
                if let Some(self_kmer) = self.kmers.get_mut(&seq) {
                    self_kmer.count += other_kmer.count;
                    let max_len = self_kmer.positions.len().max(other_kmer.positions.len());
                    self_kmer.positions.resize(max_len, 0);
                    for (i, &v) in other_kmer.positions.iter().enumerate() {
                        self_kmer.positions[i] += v;
                    }
                } else {
                    self.kmers.insert(seq, other_kmer);
                }
            }

            // Merge total_kmer_counts
            let max_len = self
                .total_kmer_counts
                .len()
                .max(other.total_kmer_counts.len());
            self.total_kmer_counts
                .resize_with(max_len, || vec![0u64; self.kmer_size]);
            for (i, other_row) in other.total_kmer_counts.iter().enumerate() {
                let self_row = &mut self.total_kmer_counts[i];
                if self_row.len() < other_row.len() {
                    self_row.resize(other_row.len(), 0);
                }
                for (j, &v) in other_row.iter().enumerate() {
                    self_row[j] += v;
                }
            }

            if other.longest_sequence > self.longest_sequence {
                self.longest_sequence = other.longest_sequence;
            }
            self.skip_count += other.skip_count;
            self.calculated = false;
        }
    }

    fn module_id(&self) -> &str {
        "kmer_content"
    }

    fn json_thresholds(&self) -> Option<serde_json::Value> {
        Some(serde_json::json!({
            "warn": self.warn_threshold,
            "error": self.error_threshold,
        }))
    }

    fn json_data(&mut self, config: &crate::config::Config) -> serde_json::Value {
        self.calculate_enrichment(config);
        let plot_series: Vec<serde_json::Value> = self
            .enrichments
            .iter()
            .enumerate()
            .map(|(i, row)| {
                let obs_exp: Vec<serde_json::Value> = row
                    .iter()
                    .map(|&v| crate::modules::json_num(v))
                    .collect();
                serde_json::json!({
                    "sequence": self.x_labels.get(i).cloned().unwrap_or_default(),
                    "obs_exp": obs_exp,
                })
            })
            .collect();
        let records: Vec<serde_json::Value> = self
            .sorted_kmer_data
            .iter()
            .map(|k| {
                serde_json::json!({
                    "sequence": k.sequence,
                    "estimated_count": k.count,
                    "p_value": crate::modules::json_num(k.p_value as f64),
                    "obs_exp_max": crate::modules::json_num(k.max_obs_exp as f64),
                    "max_obs_exp_position": k.max_position,
                })
            })
            .collect();
        serde_json::json!({
            "kmer_size": self.kmer_size,
            "sampled_every_n_reads": 50,
            "estimated_count_multiplier": 50,
            "position_groups": self.x_categories,
            "plot_series": plot_series,
            "records": records,
        })
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate_enrichment(&report.config);

        if !self.enrichments.is_empty() {
            let refs: Vec<&[f64]> = self.enrichments.iter().map(|v| v.as_slice()).collect();
            let label_refs: Vec<&str> = self.x_labels.iter().map(|s| s.as_str()).collect();
            let width = 800.max(self.x_categories.len() * 15);

            let svg = crate::graphs::line_graph::render_line_graph(
                &refs,
                0.0,
                self.max_graph_value,
                "Position in read (bp)",
                &label_refs,
                &self.x_categories,
                "Log2 Obs/Exp",
                width,
                600,
            );
            report.add_image("kmer_profiles.png", &svg);
        }

        if self.sorted_kmer_data.is_empty() {
            report.html_body.push_str("<p>No overrepresented Kmers</p>");
        } else {
            report.html_body.push_str("<table><thead><tr><th>Sequence</th><th>Count</th><th>PValue</th><th>Obs/Exp Max</th><th>Max Obs/Exp Position</th></tr></thead><tbody>");
            for kmer in &self.sorted_kmer_data {
                report.html_body.push_str(&format!(
                    "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                    kmer.sequence, kmer.count, kmer.p_value, kmer.max_obs_exp, kmer.max_position
                ));
            }
            report.html_body.push_str("</tbody></table>");
        }

        let data = &mut report.data;
        if !self.sorted_kmer_data.is_empty() {
            data.push_str("#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position\n");
            for kmer in &self.sorted_kmer_data {
                data.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    kmer.sequence, kmer.count, kmer.p_value, kmer.max_obs_exp, kmer.max_position
                ));
            }
        }

        Ok(())
    }
}
