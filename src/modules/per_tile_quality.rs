use std::collections::HashMap;

use crate::error::Result;
use crate::graphs::base_group;
use crate::graphs::quality_count::QualityCount;
use crate::modules::{ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::sequence::phred::PhredEncoding;
use crate::sequence::Sequence;

pub struct PerTileQualityScores {
    per_tile_quality_counts: HashMap<i32, Vec<QualityCount>>,
    current_length: usize,
    total_count: u64,
    split_position: Option<usize>,
    ignore_in_report: bool,
    calculated: bool,
    means: Vec<Vec<f64>>,
    x_labels: Vec<String>,
    tiles: Vec<i32>,
    max_deviation: f64,
    warn_threshold: f64,
    error_threshold: f64,
    ignore: bool,
}

impl PerTileQualityScores {
    pub fn new(module_config: &ModuleConfig) -> Self {
        Self {
            per_tile_quality_counts: HashMap::new(),
            current_length: 0,
            total_count: 0,
            split_position: None,
            ignore_in_report: false,
            calculated: false,
            means: Vec::new(),
            x_labels: Vec::new(),
            tiles: Vec::new(),
            max_deviation: 0.0,
            warn_threshold: module_config.get_param("tile", "warn"),
            error_threshold: module_config.get_param("tile", "error"),
            ignore: module_config.get_param("tile", "ignore") > 0.0,
        }
    }

    fn calculate(&mut self, config: &crate::config::Config) {
        if self.calculated {
            return;
        }

        let (min_char, _max_char) = self.calculate_offsets();
        let encoding = PhredEncoding::from_lowest_char(min_char);
        let offset = encoding.offset();

        let groups = base_group::make_base_groups(self.current_length, config);

        let mut tile_numbers: Vec<i32> = self.per_tile_quality_counts.keys().copied().collect();
        tile_numbers.sort();

        self.tiles = tile_numbers.clone();
        self.means = vec![vec![0.0; groups.len()]; tile_numbers.len()];
        self.x_labels = groups.iter().map(|g| g.to_string()).collect();

        for (t, &tile) in tile_numbers.iter().enumerate() {
            for (i, group) in groups.iter().enumerate() {
                let min_base = group.lower_count();
                let max_base = group.upper_count();
                self.means[t][i] = self.get_mean(tile, min_base, max_base, offset);
            }
        }

        // Normalize: subtract per-group average
        let mut avg_per_group = vec![0.0f64; groups.len()];
        for t in 0..tile_numbers.len() {
            for i in 0..groups.len() {
                avg_per_group[i] += self.means[t][i];
            }
        }
        for avg in avg_per_group.iter_mut() {
            *avg /= tile_numbers.len() as f64;
        }

        let mut max_dev = 0.0f64;
        for i in 0..groups.len() {
            for t in 0..tile_numbers.len() {
                self.means[t][i] -= avg_per_group[i];
                if self.means[t][i].abs() > max_dev {
                    max_dev = self.means[t][i].abs();
                }
            }
        }
        self.max_deviation = max_dev;
        self.calculated = true;
    }

    fn calculate_offsets(&self) -> (u8, u8) {
        let mut min_char: u8 = 255;
        let mut max_char: u8 = 0;
        for counts in self.per_tile_quality_counts.values() {
            for qc in counts {
                let mn = qc.min_char();
                let mx = qc.max_char();
                if min_char == 255 || mn < min_char {
                    min_char = mn;
                }
                if mx > max_char {
                    max_char = mx;
                }
            }
        }
        (min_char, max_char)
    }

    fn get_mean(&self, tile: i32, min_bp: usize, max_bp: usize, offset: usize) -> f64 {
        let quality_counts = match self.per_tile_quality_counts.get(&tile) {
            Some(qc) => qc,
            None => return 0.0,
        };
        let mut count = 0;
        let mut total = 0.0;
        for i in (min_bp - 1)..max_bp.min(quality_counts.len()) {
            if quality_counts[i].total_count() > 0 {
                count += 1;
                total += quality_counts[i].mean(offset);
            }
        }
        if count > 0 {
            total / count as f64
        } else {
            0.0
        }
    }

    fn ensure_calculated(&mut self) {
        if !self.calculated {
            let config = crate::config::Config::default_config();
            self.calculate(&config);
        }
    }
}

impl QcModule for PerTileQualityScores {
    fn name(&self) -> &str {
        "Per tile sequence quality"
    }

    fn description(&self) -> &str {
        "Shows the per tile Quality scores of all bases at a given position in a sequencing run"
    }

    fn ignore_filtered_sequences(&self) -> bool {
        true
    }

    fn ignore_in_report(&self) -> bool {
        self.ignore_in_report || self.ignore || self.current_length == 0
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if self.total_count == 0 && self.ignore {
            self.ignore_in_report = true;
        }
        if self.ignore_in_report {
            return;
        }
        if seq.quality.is_empty() {
            return;
        }

        self.calculated = false;
        self.total_count += 1;

        // Sampling: first 10K + 10% of rest
        if self.total_count > 10000 && self.total_count % 10 != 0 {
            return;
        }

        // Extract tile ID
        let split_id: Vec<&str> = seq.id.split(':').collect();
        let tile: i32;

        if let Some(pos) = self.split_position {
            if split_id.len() <= pos {
                self.ignore_in_report = true;
                return;
            }
            match split_id[pos].parse() {
                Ok(t) => tile = t,
                Err(_) => {
                    self.ignore_in_report = true;
                    return;
                }
            }
        } else if split_id.len() >= 7 {
            match split_id[4].parse() {
                Ok(t) => {
                    self.split_position = Some(4);
                    tile = t;
                }
                Err(_) => {
                    self.ignore_in_report = true;
                    return;
                }
            }
        } else if split_id.len() >= 5 {
            match split_id[2].parse() {
                Ok(t) => {
                    self.split_position = Some(2);
                    tile = t;
                }
                Err(_) => {
                    self.ignore_in_report = true;
                    return;
                }
            }
        } else {
            self.ignore_in_report = true;
            return;
        }

        let qual = seq.quality.as_bytes();
        if self.current_length < qual.len() {
            // Expand all existing tile arrays
            for counts in self.per_tile_quality_counts.values_mut() {
                counts.resize_with(qual.len(), QualityCount::new);
            }
            self.current_length = qual.len();
        }

        if !self.per_tile_quality_counts.contains_key(&tile) {
            if self.per_tile_quality_counts.len() > 2500 {
                eprintln!("Too many tiles (>2500) so giving up trying to do per-tile qualities");
                self.ignore_in_report = true;
                self.per_tile_quality_counts.clear();
                return;
            }
            let mut counts = Vec::with_capacity(self.current_length);
            for _ in 0..self.current_length {
                counts.push(QualityCount::new());
            }
            self.per_tile_quality_counts.insert(tile, counts);
        }

        let counts = self.per_tile_quality_counts.get_mut(&tile).unwrap();
        for (i, &c) in qual.iter().enumerate() {
            counts[i].add_value(c);
        }
    }

    fn reset(&mut self) {
        self.total_count = 0;
        self.per_tile_quality_counts.clear();
        self.calculated = false;
    }

    fn raises_error(&mut self) -> bool {
        self.ensure_calculated();
        self.max_deviation > self.error_threshold
    }

    fn raises_warning(&mut self) -> bool {
        self.ensure_calculated();
        self.max_deviation > self.warn_threshold
    }

    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()> {
        self.calculate(&report.config);

        let width = 800.max(self.x_labels.len() * 15);
        let svg = crate::graphs::tile_graph::render_tile_graph(
            &self.x_labels,
            &self.tiles,
            &self.means,
            width,
            600,
        );
        report.add_image("per_tile_quality.png", &svg);

        let data = &mut report.data;
        data.push_str("#Tile\tBase\tMean\n");
        for (t, &tile) in self.tiles.iter().enumerate() {
            for i in 0..self.means[t].len() {
                data.push_str(&format!("{}\t{}\t{}\n", tile, self.x_labels[i], crate::format_double(self.means[t][i])));
            }
        }

        Ok(())
    }
}
