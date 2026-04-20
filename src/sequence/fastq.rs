use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

use bzip2::read::BzDecoder;
use flate2::read::MultiGzDecoder;

use crate::config::Config;
use crate::error::{FastqcError, Result};
use crate::sequence::colorspace;
use crate::sequence::{Sequence, SequenceFile};

/// Buffer size for underlying I/O (256 KB).
const IO_BUF_SIZE: usize = 256 * 1024;

/// Block size for the record buffer (1 MB).
const BLOCK_SIZE: usize = 1024 * 1024;

/// Zero-copy FASTQ reader.
///
/// Instead of 4 `String` allocations per record (id, seq, sep, qual),
/// this reader fills a large byte buffer and locates newlines with `memchr`.
/// Records are constructed from `Vec<u8>` → `String::from_utf8_unchecked`
/// instead of per-char `read_line` + `String::new` + trim.
pub struct FastqReader {
    reader: Box<dyn Read + Send>,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
    eof: bool,

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
                reader: Box::new(BufReader::with_capacity(IO_BUF_SIZE, io::stdin())),
                buf: Vec::with_capacity(BLOCK_SIZE * 2),
                pos: 0,
                len: 0,
                eof: false,
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
        let reader: Box<dyn Read + Send> = if path_str.ends_with(".gz") {
            Box::new(BufReader::with_capacity(
                IO_BUF_SIZE,
                MultiGzDecoder::new(file),
            ))
        } else if path_str.ends_with(".bz2") {
            Box::new(BufReader::with_capacity(IO_BUF_SIZE, BzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(IO_BUF_SIZE, file))
        };

        Ok(Self {
            reader,
            buf: Vec::with_capacity(BLOCK_SIZE * 2),
            pos: 0,
            len: 0,
            eof: false,
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

    /// Ensure at least `needed` bytes are available from `pos`.
    fn fill_buf(&mut self, needed: usize) -> io::Result<()> {
        if self.len - self.pos >= needed {
            return Ok(());
        }
        if self.eof {
            return Ok(());
        }
        if self.pos > 0 {
            self.buf.copy_within(self.pos..self.len, 0);
            self.len -= self.pos;
            self.pos = 0;
        }
        let target = (self.len + BLOCK_SIZE).max(self.len + needed);
        if self.buf.len() < target {
            self.buf.resize(target, 0);
        }
        loop {
            let n = self.reader.read(&mut self.buf[self.len..target])?;
            if n == 0 {
                self.eof = true;
                break;
            }
            self.len += n;
            self.bytes_read += n as u64;
            if self.len - self.pos >= needed {
                break;
            }
        }
        Ok(())
    }

    /// Read one line, returning an *owned* `Vec<u8>` (no trailing newline/CR).
    /// This avoids the borrow-checker issues of returning `&[u8]`.
    /// The copy is from a contiguous buffer slice → memcpy, which is very fast.
    fn next_line_owned(&mut self) -> Result<Option<Vec<u8>>> {
        self.fill_buf(1).map_err(FastqcError::Io)?;
        if self.pos == self.len {
            return Ok(None);
        }

        loop {
            if let Some(offset) = memchr::memchr(b'\n', &self.buf[self.pos..self.len]) {
                let start = self.pos;
                let mut end = self.pos + offset;
                if end > start && self.buf[end - 1] == b'\r' {
                    end -= 1;
                }
                let line = self.buf[start..end].to_vec();
                self.pos = start + offset + 1;
                self.line_number += 1;
                return Ok(Some(line));
            }
            if self.eof {
                if self.pos < self.len {
                    let start = self.pos;
                    let mut end = self.len;
                    if end > start && self.buf[end - 1] == b'\r' {
                        end -= 1;
                    }
                    let line = self.buf[start..end].to_vec();
                    self.pos = self.len;
                    self.line_number += 1;
                    return Ok(Some(line));
                }
                return Ok(None);
            }
            self.fill_buf(self.len - self.pos + BLOCK_SIZE)
                .map_err(FastqcError::Io)?;
        }
    }

    /// Read one line but only check and discard it (for the '+' separator).
    /// Returns true if the line starts with '+', false otherwise.
    fn skip_separator_line(&mut self) -> Result<bool> {
        self.fill_buf(1).map_err(FastqcError::Io)?;
        if self.pos == self.len {
            return Err(FastqcError::SequenceFormat(
                "Unexpected end of file reading separator".into(),
            ));
        }

        loop {
            if let Some(offset) = memchr::memchr(b'\n', &self.buf[self.pos..self.len]) {
                let starts_with_plus = self.buf[self.pos] == b'+';
                if !starts_with_plus {
                    let mut end = self.pos + offset;
                    if end > self.pos && self.buf[end - 1] == b'\r' {
                        end -= 1;
                    }
                    let display = String::from_utf8_lossy(&self.buf[self.pos..end]).into_owned();
                    self.pos = self.pos + offset + 1;
                    self.line_number += 1;
                    return Err(FastqcError::SequenceFormat(format!(
                        "Line {} expected to start with +, got: {}",
                        self.line_number, display
                    )));
                }
                self.pos = self.pos + offset + 1;
                self.line_number += 1;
                return Ok(true);
            }
            if self.eof {
                if self.pos < self.len {
                    let starts_with_plus = self.buf[self.pos] == b'+';
                    self.pos = self.len;
                    self.line_number += 1;
                    return Ok(starts_with_plus);
                }
                return Err(FastqcError::SequenceFormat(
                    "Unexpected end of file reading separator".into(),
                ));
            }
            self.fill_buf(self.len - self.pos + BLOCK_SIZE)
                .map_err(FastqcError::Io)?;
        }
    }

    pub fn read_sequence(&mut self) -> Result<Option<Sequence>> {
        // --- Line 1: ID ---
        let id_bytes = loop {
            match self.next_line_owned()? {
                None => return Ok(None),
                Some(line) if line.is_empty() => continue,
                Some(line) => break line,
            }
        };

        if id_bytes.first() != Some(&b'@') {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {} expected to start with @, got: {}",
                self.line_number,
                String::from_utf8_lossy(&id_bytes)
            )));
        }

        let id = String::from_utf8(id_bytes[1..].to_vec()).map_err(|e| {
            FastqcError::SequenceFormat(format!(
                "Line {}: sequence ID contains invalid UTF-8: {}",
                self.line_number, e
            ))
        })?;

        // --- Line 2: Sequence ---
        let seq_bytes = self.next_line_owned()?.ok_or_else(|| {
            FastqcError::SequenceFormat("Unexpected end of file reading sequence".into())
        })?;

        // --- Line 3: Separator (+) — check and discard ---
        self.skip_separator_line()?;

        // --- Line 4: Quality ---
        let qual_bytes = self.next_line_owned()?.ok_or_else(|| {
            FastqcError::SequenceFormat("Unexpected end of file reading quality".into())
        })?;

        if qual_bytes.len() != seq_bytes.len() {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {}: quality length ({}) != sequence length ({})",
                self.line_number,
                qual_bytes.len(),
                seq_bytes.len()
            )));
        }

        // FASTQ spec requires ASCII for sequence and quality lines, but we must
        // not assume that on untrusted input: corrupted files or files written
        // with legacy non-UTF-8 encodings would trigger undefined behaviour if
        // wrapped with `from_utf8_unchecked`. Validate ASCII (cheap, tighter
        // than UTF-8) and only then perform the unchecked conversion.
        if !seq_bytes.is_ascii() {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {}: sequence line contains non-ASCII bytes",
                self.line_number.saturating_sub(2)
            )));
        }
        if !qual_bytes.is_ascii() {
            return Err(FastqcError::SequenceFormat(format!(
                "Line {}: quality line contains non-ASCII bytes",
                self.line_number
            )));
        }

        // SAFETY: both byte slices are pure ASCII (checked above), therefore
        // they are also valid UTF-8. Avoiding `String::from_utf8` here saves
        // the redundant re-validation on the hot path.
        let seq_line = unsafe { String::from_utf8_unchecked(seq_bytes) };
        let qual_line = unsafe { String::from_utf8_unchecked(qual_bytes) };

        // Handle colorspace
        let seq = if !self.colorspace_checked {
            self.colorspace_checked = true;
            if colorspace::is_colorspace(&seq_line) {
                self.is_colorspace = true;
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
