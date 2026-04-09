use std::path::Path;

use crate::config::Config;
use crate::error::Result;
use crate::sequence::{Sequence, SequenceFile};

/// BAM/SAM file reader (stub - requires noodles dependency for full implementation).
pub struct BamReader {
    name: String,
    finished: bool,
}

impl BamReader {
    pub fn new(path: &Path, _config: &Config, _only_mapped: bool) -> Result<Self> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        // TODO: Implement BAM/SAM parsing using noodles
        eprintln!("BAM/SAM support is not yet fully implemented. File: {}", name);

        Ok(Self {
            name,
            finished: true,
        })
    }
}

impl Iterator for BamReader {
    type Item = Result<Sequence>;

    fn next(&mut self) -> Option<Self::Item> {
        None // Stub
    }
}

impl SequenceFile for BamReader {
    fn name(&self) -> &str {
        &self.name
    }

    fn percent_complete(&self) -> u8 {
        100
    }

    fn is_colorspace(&self) -> bool {
        false
    }
}
