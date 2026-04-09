use crate::error::Result;
use crate::graphs::base_group;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::Sequence;

pub struct PerBaseSequenceContent {
    g_counts: Vec<u64>,
    a_counts: Vec<u64>,
    t_counts: Vec<u64>,
    c_counts: Vec<u64>,
    calculated: bool,
    g_pcts: Vec<f64>,
    a_pcts: Vec<f64>,
    t_pcts: Vec<f64>,
    c_pcts: Vec<f64>,
    x_labels: Vec<String>,
    ignore: bool,
    warn_threshold: f64,
    error_threshold: f64,
}

impl PerBaseSequenceContent {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            g_counts: Vec::new(),
            a_counts: Vec::new(),
            t_counts: Vec::new(),
            c_counts: Vec::new(),
            calculated: false,
            g_pcts: Vec::new(),
            a_pcts: Vec::new(),
            t_pcts: Vec::new(),
            c_pcts: Vec::new(),
            x_labels: Vec::new(),
            ignore: module_config.get_param("sequence", "ignore") > 0.0,
            warn_threshold: module_config.get_param("sequence", "warn"),
            error_threshold: module_config.get_param("sequence", "error"),
        }
    }

    fn calculate(&mut self, config: &crate::config::Config) {
        if self.calculated {
            return;
        }

        let groups = base_group::make_base_groups(self.g_counts.len(), config);
        let mut g_pcts = Vec::with_capacity(groups.len());
        let mut a_pcts = Vec::with_capacity(groups.len());
        let mut t_pcts = Vec::with_capacity(groups.len());
        let mut c_pcts = Vec::with_capacity(groups.len());
        let mut x_labels = Vec::with_capacity(groups.len());

        for group in &groups {
            x_labels.push(group.to_string());
            let mut g_total = 0u64;
            let mut a_total = 0u64;
            let mut t_total = 0u64;
            let mut c_total = 0u64;

            for i in (group.lower_count() - 1)..group.upper_count().min(self.g_counts.len()) {
                g_total += self.g_counts[i];
                a_total += self.a_counts[i];
                t_total += self.t_counts[i];
                c_total += self.c_counts[i];
            }

            let total = g_total + a_total + t_total + c_total;
            if total > 0 {
                g_pcts.push(g_total as f64 * 100.0 / total as f64);
                a_pcts.push(a_total as f64 * 100.0 / total as f64);
                t_pcts.push(t_total as f64 * 100.0 / total as f64);
                c_pcts.push(c_total as f64 * 100.0 / total as f64);
            } else {
                g_pcts.push(0.0);
                a_pcts.push(0.0);
                t_pcts.push(0.0);
                c_pcts.push(0.0);
            }
        }

        self.g_pcts = g_pcts;
        self.a_pcts = a_pcts;
        self.t_pcts = t_pcts;
        self.c_pcts = c_pcts;
        self.x_labels = x_labels;
        self.calculated = true;
    }

    fn check_threshold(&mut self, threshold: f64) -> bool {
        if !self.calculated {
            let config = crate::config::Config::default_config();
            self.calculate(&config);
        }
        for i in 0..self.g_pcts.len() {
            let gc_diff = (self.g_pcts[i] - self.c_pcts[i]).abs();
            let at_diff = (self.a_pcts[i] - self.t_pcts[i]).abs();
            if gc_diff > threshold || at_diff > threshold {
                return true;
            }
        }
        false
    }
}

impl QcModule for PerBaseSequenceContent {
    fn name(&self) -> &str {
        "Per base sequence content"
    }

    fn description(&self) -> &str {
        "Shows the proportion of each base at each position in a sequence"
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

        if self.g_counts.len() < len {
            self.g_counts.resize(len, 0);
            self.a_counts.resize(len, 0);
            self.t_counts.resize(len, 0);
            self.c_counts.resize(len, 0);
        }

        for (i, &c) in seq_bytes.iter().enumerate() {
            match c {
                b'G' => self.g_counts[i] += 1,
                b'A' => self.a_counts[i] += 1,
                b'T' => self.t_counts[i] += 1,
                b'C' => self.c_counts[i] += 1,
                _ => {}
            }
        }
    }

    fn reset(&mut self) {
        self.g_counts.clear();
        self.a_counts.clear();
        self.t_counts.clear();
        self.c_counts.clear();
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        self.check_threshold(self.error_threshold)
    }

    fn raises_warning(&mut self) -> bool {
        self.check_threshold(self.warn_threshold)
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate(&report.config);

        let svg = crate::graphs::line_graph::render_line_graph(
            &[&self.g_pcts, &self.a_pcts, &self.t_pcts, &self.c_pcts],
            0.0,
            100.0,
            "Position in read (bp)",
            &["%G", "%A", "%T", "%C"],
            &self.x_labels,
            "Sequence content across all bases",
            800,
            600,
        );
        report.add_image("per_base_sequence_content.png", &svg);

        let data = &mut report.data;
        data.push_str("#Base\tG\tA\tT\tC\n");
        for i in 0..self.x_labels.len() {
            data.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                self.x_labels[i],
                crate::format_double(self.g_pcts[i]),
                crate::format_double(self.a_pcts[i]),
                crate::format_double(self.t_pcts[i]),
                crate::format_double(self.c_pcts[i])
            ));
        }

        Ok(())
    }
}
