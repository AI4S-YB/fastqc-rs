use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use bzip2::read::BzDecoder;

use crate::config::Config;
use crate::error::{FastqcError, Result};
use crate::sequence::{Sequence, SequenceFile};
use crate::sequence::colorspace;

/// FASTQ file reader supporting plain text, gzip, and bzip2 compression.
pub struct FastqReader {
    reader: Box<dyn BufRead>,
    name: String,
    file_size: u64,
    bytes_read: u64,
    line_number: usize,
    is_colorspace: bool,
    colorspace_checked: bool,
    casava: bool,
    nofilter: bool,
    finished: bool,
}

impl FastqReader {
    pub fn new(path: &Path, config: &Config) -> Result<Self> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        if name == "stdin" || path.to_string_lossy() == "stdin" {
            return Ok(Self {
                reader: Box::new(BufReader::new(io::stdin())),
                name: "stdin".to_string(),
                file_size: u64::MAX,
                bytes_read: 0,
                line_number: 0,
                is_colorspace: false,
                colorspace_checked: false,
                casava: config.casava,
                nofilter: config.nofilter,
                finished: false,
            });
        }

        let file_size = std::fs::metadata(path)?.len();
        let file = File::open(path)?;

        let path_str = path.to_string_lossy().to_lowercase();
        let reader: Box<dyn BufRead> = if path_str.ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else if path_str.ends_with(".bz2") {
            Box::new(BufReader::new(BzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(Self {
            reader,
            name,
            file_size,
            bytes_read: 0,
            line_number: 0,
            is_colorspace: false,
            colorspace_checked: false,
            casava: config.casava,
            nofilter: config.nofilter,
            finished: false,
        })
    }

    fn read_line(&mut self) -> Result<Option<String>> {
        let mut line = String::new();
        let n = self.reader.read_line(&mut line)?;
        if n == 0 {
            return Ok(None);
        }
        self.bytes_read += n as u64;
        self.line_number += 1;

        // Trim trailing newline/carriage return
        while line.ends_with('\n') || line.ends_with('\r') {
            line.pop();
        }
        Ok(Some(line))
    }

    fn read_sequence(&mut self) -> Result<Option<Sequence>> {
        // Skip blank lines to find ID line
        let id_line = loop {
            match self.read_line()? {
                None => return Ok(None),
                Some(line) if line.is_empty() => continue,
                Some(line) => break line,
            }
        };

        if !id_line.starts_with('@') {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {} expected to start with @, got: {}",
                self.line_number, id_line
            )));
        }

        let id = id_line[1..].to_string();

        // Read sequence line
        let seq_line = self
            .read_line()?
            .ok_or_else(|| FastqcError::SequenceFormat("Unexpected end of file reading sequence".into()))?;

        // Read separator line (+)
        let sep_line = self
            .read_line()?
            .ok_or_else(|| FastqcError::SequenceFormat("Unexpected end of file reading separator".into()))?;

        if !sep_line.starts_with('+') {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {} expected to start with +, got: {}",
                self.line_number, sep_line
            )));
        }

        // Read quality line
        let qual_line = self
            .read_line()?
            .ok_or_else(|| FastqcError::SequenceFormat("Unexpected end of file reading quality".into()))?;

        if qual_line.len() != seq_line.len() {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {}: quality length ({}) != sequence length ({})",
                self.line_number,
                qual_line.len(),
                seq_line.len()
            )));
        }

        // Handle colorspace
        let mut seq = if !self.colorspace_checked {
            self.colorspace_checked = true;
            if colorspace::is_colorspace(&seq_line) {
                self.is_colorspace = true;
                let converted = colorspace::colorspace_to_bases(&seq_line);
                let mut s = Sequence::with_colorspace(
                    id,
                    converted,
                    seq_line,
                    qual_line[1..].to_string(), // Colorspace quality is 1 shorter
                    self.name.clone(),
                );
                // CASAVA filtering
                if self.casava && !self.nofilter && s.id.contains(":Y:") {
                    s.is_filtered = true;
                }
                return Ok(Some(s));
            }
            let mut s = Sequence::new(id, seq_line, qual_line, self.name.clone());
            if self.casava && !self.nofilter && s.id.contains(":Y:") {
                s.is_filtered = true;
            }
            s
        } else if self.is_colorspace {
            let converted = colorspace::colorspace_to_bases(&seq_line);
            let mut s = Sequence::with_colorspace(
                id,
                converted,
                seq_line,
                qual_line[1..].to_string(),
                self.name.clone(),
            );
            if self.casava && !self.nofilter && s.id.contains(":Y:") {
                s.is_filtered = true;
            }
            s
        } else {
            let mut s = Sequence::new(id, seq_line, qual_line, self.name.clone());
            if self.casava && !self.nofilter && s.id.contains(":Y:") {
                s.is_filtered = true;
            }
            s
        };

        Ok(Some(seq))
    }
}

impl Iterator for FastqReader {
    type Item = Result<Sequence>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.read_sequence() {
            Ok(Some(seq)) => Some(Ok(seq)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(e) => {
                self.finished = true;
                Some(Err(e))
            }
        }
    }
}

impl SequenceFile for FastqReader {
    fn name(&self) -> &str {
        &self.name
    }

    fn percent_complete(&self) -> u8 {
        if self.file_size == 0 || self.file_size == u64::MAX {
            return 0;
        }
        ((self.bytes_read * 100) / self.file_size).min(100) as u8
    }

    fn is_colorspace(&self) -> bool {
        self.is_colorspace
    }
}
