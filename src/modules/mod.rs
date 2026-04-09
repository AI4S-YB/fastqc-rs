pub mod adapter_content;
pub mod basic_stats;
pub mod duplication_level;
pub mod gc_model;
pub mod kmer_content;
pub mod n_content;
pub mod overrepresented_seqs;
pub mod per_base_quality;
pub mod per_base_sequence_content;
pub mod per_sequence_gc_content;
pub mod per_sequence_quality;
pub mod per_tile_quality;
pub mod sequence_length_distribution;

use crate::config::Config;
use crate::error::Result;
use crate::report::ReportArchive;
use crate::sequence::Sequence;

/// Quality control module trait, equivalent to Java's QCModule interface.
/// Requires Send for multi-threaded data-parallel processing.
pub trait QcModule: Send {
    fn name(&self) -> &str;
    fn description(&self) -> &str;
    fn process_sequence(&mut self, seq: &Sequence);
    fn reset(&mut self);
    fn raises_error(&mut self) -> bool;
    fn raises_warning(&mut self) -> bool;
    fn ignore_filtered_sequences(&self) -> bool;
    fn ignore_in_report(&self) -> bool;
    fn make_report(&mut self, report: &mut ReportArchive) -> Result<()>;

    /// Convert self into Any for downcasting during merge.
    fn into_any(self: Box<Self>) -> Box<dyn std::any::Any + Send>;

    /// Merge another module's accumulated state into this one.
    /// The `other` must be the same concrete type (obtained via `into_any`).
    fn merge(&mut self, other: Box<dyn std::any::Any + Send>);
}

/// Module threshold configuration, loaded from limits.txt.
pub struct ModuleConfig {
    params: std::collections::HashMap<String, f64>,
}

impl ModuleConfig {
    pub fn load(config: &Config) -> Self {
        let mut params = std::collections::HashMap::new();

        // Set defaults (matching Java's ModuleConfig)
        params.insert("duplication:warn".into(), 70.0);
        params.insert("duplication:error".into(), 50.0);
        params.insert("kmer:warn".into(), 2.0);
        params.insert("kmer:error".into(), 5.0);
        params.insert("n_content:warn".into(), 5.0);
        params.insert("n_content:error".into(), 20.0);
        params.insert("overrepresented:warn".into(), 0.1);
        params.insert("overrepresented:error".into(), 1.0);
        params.insert("quality_base_lower:warn".into(), 10.0);
        params.insert("quality_base_lower:error".into(), 5.0);
        params.insert("quality_base_median:warn".into(), 25.0);
        params.insert("quality_base_median:error".into(), 20.0);
        params.insert("sequence:warn".into(), 10.0);
        params.insert("sequence:error".into(), 20.0);
        params.insert("gc_sequence:warn".into(), 15.0);
        params.insert("gc_sequence:error".into(), 30.0);
        params.insert("quality_sequence:warn".into(), 20.0);
        params.insert("quality_sequence:error".into(), 27.0);
        params.insert("tile:warn".into(), 5.0);
        params.insert("tile:error".into(), 10.0);
        params.insert("sequence_length:warn".into(), 1.0);
        params.insert("sequence_length:error".into(), 1.0);
        params.insert("adapter:warn".into(), 5.0);
        params.insert("adapter:error".into(), 10.0);

        // Ignore flags default to 0
        for module in &[
            "duplication",
            "kmer",
            "n_content",
            "overrepresented",
            "quality_base",
            "sequence",
            "gc_sequence",
            "quality_sequence",
            "tile",
            "sequence_length",
            "adapter",
        ] {
            params.insert(format!("{}:ignore", module), 0.0);
        }

        // Load from limits file
        let limits_text = if let Some(ref path) = config.limits_file {
            std::fs::read_to_string(path).unwrap_or_default()
        } else {
            include_str!("../data/limits.txt").to_string()
        };

        for line in limits_text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                eprintln!(
                    "Config line '{}' didn't contain the 3 required sections",
                    line
                );
                continue;
            }
            if !["warn", "error", "ignore"].contains(&parts[1]) {
                eprintln!(
                    "Second config field must be error, warn or ignore, not '{}'",
                    parts[1]
                );
                continue;
            }
            let value: f64 = match parts[2].parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("Value {} didn't look like a number", parts[2]);
                    continue;
                }
            };
            params.insert(format!("{}:{}", parts[0], parts[1]), value);
        }

        Self { params }
    }

    pub fn get_param(&self, module: &str, level: &str) -> f64 {
        let key = format!("{}:{}", module, level);
        *self.params.get(&key).unwrap_or(&0.0)
    }
}

/// Create the standard module list in the same order as Java FastQC.
pub fn create_module_list(config: &Config, module_config: &ModuleConfig) -> Vec<Box<dyn QcModule>> {
    let overrep = overrepresented_seqs::OverRepresentedSeqs::new(config, module_config);
    let dup_data = overrep.shared_data();
    let dup = duplication_level::DuplicationLevel::new(dup_data.clone(), module_config);

    let mut modules: Vec<Box<dyn QcModule>> = Vec::new();
    modules.push(Box::new(basic_stats::BasicStats::new()));
    modules.push(Box::new(per_base_quality::PerBaseQualityScores::new(
        module_config,
    )));
    modules.push(Box::new(per_tile_quality::PerTileQualityScores::new(
        module_config,
    )));
    modules.push(Box::new(
        per_sequence_quality::PerSequenceQualityScores::new(module_config),
    ));
    modules.push(Box::new(
        per_base_sequence_content::PerBaseSequenceContent::new(module_config),
    ));
    modules.push(Box::new(
        per_sequence_gc_content::PerSequenceGCContent::new(module_config),
    ));
    modules.push(Box::new(n_content::NContent::new(module_config)));
    modules.push(Box::new(
        sequence_length_distribution::SequenceLengthDistribution::new(module_config),
    ));
    modules.push(Box::new(dup));
    modules.push(Box::new(overrep));
    modules.push(Box::new(adapter_content::AdapterContent::new(
        config,
        module_config,
    )));
    modules.push(Box::new(kmer_content::KmerContent::new(
        config,
        module_config,
    )));
    modules
}
