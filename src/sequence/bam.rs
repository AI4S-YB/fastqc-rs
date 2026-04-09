use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind as CigarOp;

use crate::config::Config;
use crate::error::{FastqcError, Result};
use crate::sequence::{Sequence, SequenceFile};

enum BamReaderInner {
    Bam {
        reader: bam::io::Reader<noodles::bgzf::io::Reader<File>>,
        record: bam::Record,
    },
    Sam {
        lines: Vec<sam::Record>,
        index: usize,
    },
}

pub struct BamReader {
    inner: BamReaderInner,
    name: String,
    file_size: u64,
    bytes_read_approx: u64,
    only_mapped: bool,
    finished: bool,
}

impl BamReader {
    pub fn new(path: &Path, _config: &Config, only_mapped: bool) -> Result<Self> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let file_size = std::fs::metadata(path)?.len();
        let path_lower = path.to_string_lossy().to_lowercase();
        let is_sam = path_lower.ends_with(".sam");

        let inner = if is_sam {
            let file = File::open(path)?;
            let buf = BufReader::new(file);
            let mut reader = sam::io::Reader::new(buf);
            let _header = reader.read_header().map_err(|e| {
                FastqcError::SequenceFormat(format!("Failed to read SAM header: {}", e))
            })?;
            let records: Vec<sam::Record> = reader.records().filter_map(|r| r.ok()).collect();
            BamReaderInner::Sam {
                lines: records,
                index: 0,
            }
        } else {
            let file = File::open(path)?;
            let mut reader = bam::io::Reader::new(file);
            let _header = reader.read_header().map_err(|e| {
                FastqcError::SequenceFormat(format!("Failed to read BAM header: {}", e))
            })?;
            BamReaderInner::Bam {
                reader,
                record: bam::Record::default(),
            }
        };

        Ok(Self {
            inner,
            name,
            file_size,
            bytes_read_approx: 0,
            only_mapped,
            finished: false,
        })
    }

    fn process_bam_record(
        record: &bam::Record,
        only_mapped: bool,
        file_name: &str,
    ) -> Result<Option<Sequence>> {
        let flags = record.flags();

        // Skip unmapped if only_mapped
        if only_mapped && flags.is_unmapped() {
            return Ok(None);
        }

        // Read name
        let read_name = record
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "unknown".to_string());

        // Sequence: bam sequence iter yields u8 base codes (A/C/G/T/N as ASCII)
        let bam_seq = record.sequence();
        let mut sequence: String = bam_seq.iter().map(|b| b as char).collect();

        // Quality scores: bam quality iter yields raw phred values (0-based)
        let bam_qual = record.quality_scores();
        let mut qualities: String = bam_qual.iter().map(|q| (q + 33) as char).collect();

        // If no quality data, fill with '!'
        if qualities.is_empty() && !sequence.is_empty() {
            qualities = "!".repeat(sequence.len());
        }

        // CIGAR soft clipping for mapped-only mode
        if only_mapped && !flags.is_unmapped() {
            let cigar = record.cigar();
            let ops: Vec<(CigarOp, usize)> = cigar
                .iter()
                .filter_map(|r| r.ok())
                .map(|op| (op.kind(), op.len()))
                .collect();

            if !ops.is_empty() {
                // Clip 3' end first
                if let Some(&(CigarOp::SoftClip, clip_len)) = ops.last() {
                    let new_len = sequence.len().saturating_sub(clip_len);
                    sequence.truncate(new_len);
                    qualities.truncate(new_len);
                }
                // Then clip 5' end
                if let Some(&(CigarOp::SoftClip, clip_len)) = ops.first() {
                    if clip_len <= sequence.len() {
                        sequence = sequence[clip_len..].to_string();
                        qualities = qualities[clip_len..].to_string();
                    }
                }
            }
        }

        // Reverse complement if on negative strand
        if flags.is_reverse_complemented() {
            sequence = reverse_complement(&sequence);
            qualities = reverse_string(&qualities);
        }

        Ok(Some(Sequence::new(
            read_name,
            sequence,
            qualities,
            file_name.to_string(),
        )))
    }

    fn process_sam_record(
        record: &sam::Record,
        only_mapped: bool,
        file_name: &str,
    ) -> Result<Option<Sequence>> {
        let flags = record
            .flags()
            .map_err(|e| FastqcError::SequenceFormat(format!("SAM flags error: {}", e)))?;

        if only_mapped && flags.is_unmapped() {
            return Ok(None);
        }

        let read_name = record
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "unknown".to_string());

        // SAM sequence - need trait import for .len() and .iter()
        let sam_seq = record.sequence();
        let seq_len = noodles::sam::alignment::record::Sequence::len(&sam_seq);
        let mut sequence: String = noodles::sam::alignment::record::Sequence::iter(&sam_seq)
            .map(|b| b as char)
            .collect();

        // SAM quality scores
        let sam_qual = record.quality_scores();
        let mut qualities: String =
            noodles::sam::alignment::record::QualityScores::iter(&sam_qual)
                .filter_map(|r| r.ok())
                .map(|q| (q + 33) as char)
                .collect();

        if qualities.is_empty() && seq_len > 0 {
            qualities = "!".repeat(seq_len);
        }

        if only_mapped && !flags.is_unmapped() {
            let cigar = record.cigar();
            let ops: Vec<(CigarOp, usize)> =
                noodles::sam::alignment::record::Cigar::iter(&cigar)
                    .filter_map(|r| r.ok())
                    .map(|op| (op.kind(), op.len()))
                    .collect();

            if !ops.is_empty() {
                if let Some(&(CigarOp::SoftClip, clip_len)) = ops.last() {
                    let new_len = sequence.len().saturating_sub(clip_len);
                    sequence.truncate(new_len);
                    qualities.truncate(new_len);
                }
                if let Some(&(CigarOp::SoftClip, clip_len)) = ops.first() {
                    if clip_len <= sequence.len() {
                        sequence = sequence[clip_len..].to_string();
                        qualities = qualities[clip_len..].to_string();
                    }
                }
            }
        }

        if flags.is_reverse_complemented() {
            sequence = reverse_complement(&sequence);
            qualities = reverse_string(&qualities);
        }

        Ok(Some(Sequence::new(
            read_name,
            sequence,
            qualities,
            file_name.to_string(),
        )))
    }
}

impl Iterator for BamReader {
    type Item = Result<Sequence>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        loop {
            match &mut self.inner {
                BamReaderInner::Bam { reader, record } => match reader.read_record(record) {
                    Ok(0) => {
                        self.finished = true;
                        return None;
                    }
                    Ok(n) => {
                        self.bytes_read_approx += n as u64;
                        match Self::process_bam_record(record, self.only_mapped, &self.name) {
                            Ok(Some(seq)) => return Some(Ok(seq)),
                            Ok(None) => continue,
                            Err(e) => {
                                self.finished = true;
                                return Some(Err(e));
                            }
                        }
                    }
                    Err(e) => {
                        self.finished = true;
                        return Some(Err(FastqcError::SequenceFormat(format!(
                            "BAM read error: {}",
                            e
                        ))));
                    }
                },
                BamReaderInner::Sam { lines, index } => {
                    if *index >= lines.len() {
                        self.finished = true;
                        return None;
                    }
                    let record = &lines[*index];
                    *index += 1;
                    match Self::process_sam_record(record, self.only_mapped, &self.name) {
                        Ok(Some(seq)) => return Some(Ok(seq)),
                        Ok(None) => continue,
                        Err(e) => {
                            self.finished = true;
                            return Some(Err(e));
                        }
                    }
                }
            }
        }
    }
}

impl SequenceFile for BamReader {
    fn name(&self) -> &str {
        &self.name
    }

    fn percent_complete(&self) -> u8 {
        if self.file_size == 0 {
            return 0;
        }
        match &self.inner {
            BamReaderInner::Sam { lines, index } => {
                if lines.is_empty() {
                    100
                } else {
                    ((*index as u64 * 100) / lines.len() as u64).min(100) as u8
                }
            }
            BamReaderInner::Bam { .. } => {
                ((self.bytes_read_approx * 100) / self.file_size).min(100) as u8
            }
        }
    }

    fn is_colorspace(&self) -> bool {
        false
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            other => other,
        })
        .collect()
}

fn reverse_string(s: &str) -> String {
    s.chars().rev().collect()
}
