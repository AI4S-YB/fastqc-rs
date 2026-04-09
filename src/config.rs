use std::path::PathBuf;

use clap::Parser;

/// A quality control tool for high throughput sequence data
#[derive(Parser, Debug, Clone)]
#[command(name = "fastqc-rs", version, about)]
pub struct Config {
    /// Input files to process
    #[arg(required = true)]
    pub files: Vec<PathBuf>,

    /// Create an output directory for the results
    #[arg(short = 'o', long = "outdir")]
    pub output_dir: Option<PathBuf>,

    /// Unzip the output file after creating it
    #[arg(long)]
    pub extract: bool,

    /// Do not uncompress the output file after creating it
    #[arg(long)]
    pub noextract: bool,

    /// Delete the zip file after extraction
    #[arg(long)]
    pub delete: bool,

    /// Bypass the normal sequence file format detection and force the
    /// program to use the specified format
    #[arg(short = 'f', long = "format")]
    pub sequence_format: Option<String>,

    /// Specifies a non-default file which contains the list of contaminants
    #[arg(short = 'c', long = "contaminants")]
    pub contaminant_file: Option<PathBuf>,

    /// Specifies a non-default file which contains the list of adapter sequences
    #[arg(short = 'a', long = "adapters")]
    pub adapter_file: Option<PathBuf>,

    /// Specifies a non-default file which contains a set of criteria for
    /// pass/warn/fail for each module
    #[arg(short = 'l', long = "limits")]
    pub limits_file: Option<PathBuf>,

    /// Specifies the number of files which can be processed simultaneously
    #[arg(short = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Sets the base amount of memory, in Megabytes, used to process each file
    #[arg(long = "memory", default_value = "512")]
    pub memory: usize,

    /// Specifies the length of Kmer to look for
    #[arg(short = 'k', long = "kmers")]
    pub kmer_size: Option<usize>,

    /// Suppress all progress messages on stdout
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Files come from raw casava output
    #[arg(long)]
    pub casava: bool,

    /// Files come from nanopore sequences and are in fast5 format
    #[arg(long)]
    pub nano: bool,

    /// If running with --casava then don't remove read flagged by casava as poor quality
    #[arg(long)]
    pub nofilter: bool,

    /// Disable grouping of bases for reads >50bp
    #[arg(long)]
    pub nogroup: bool,

    /// Use exponential base grouping
    #[arg(long)]
    pub expgroup: bool,

    /// Sets an artificial lower limit on the length of the sequence to be shown
    #[arg(long, default_value = "0")]
    pub min_length: usize,

    /// Sets a length to which sequences are truncated for duplication analysis
    #[arg(long, default_value = "0")]
    pub dup_length: usize,

    /// Output SVG graphs instead of PNG
    #[arg(long)]
    pub svg: bool,

    /// Selects a directory to be used for temporary files
    #[arg(short = 'd', long = "dir")]
    pub temp_dir: Option<PathBuf>,
}

impl Config {
    /// Create a default config (for internal use when no CLI args available).
    pub fn default_config() -> Self {
        Self {
            files: vec![],
            output_dir: None,
            extract: false,
            noextract: false,
            delete: false,
            sequence_format: None,
            contaminant_file: None,
            adapter_file: None,
            limits_file: None,
            threads: 1,
            memory: 512,
            kmer_size: None,
            quiet: false,
            casava: false,
            nano: false,
            nofilter: false,
            nogroup: false,
            expgroup: false,
            min_length: 0,
            dup_length: 0,
            svg: false,
            temp_dir: None,
        }
    }

    pub fn validate(&self) -> crate::error::Result<()> {
        if let Some(ref dir) = self.output_dir {
            if !dir.exists() || !dir.is_dir() {
                return Err(crate::error::FastqcError::Config(format!(
                    "Output directory '{}' does not exist or is not a directory",
                    dir.display()
                )));
            }
        }

        if self.threads < 1 {
            return Err(crate::error::FastqcError::Config(
                "Thread count must be >= 1".into(),
            ));
        }

        if let Some(ref fmt) = self.sequence_format {
            match fmt.as_str() {
                "fastq" | "bam" | "sam" | "bam_mapped" | "sam_mapped" => {}
                _ => {
                    return Err(crate::error::FastqcError::Config(format!(
                        "Unknown sequence format '{}'. Valid options: fastq, bam, sam, bam_mapped, sam_mapped",
                        fmt
                    )));
                }
            }
        }

        if let Some(ref f) = self.contaminant_file {
            if !f.exists() {
                return Err(crate::error::FastqcError::FileNotFound(f.clone()));
            }
        }
        if let Some(ref f) = self.adapter_file {
            if !f.exists() {
                return Err(crate::error::FastqcError::FileNotFound(f.clone()));
            }
        }
        if let Some(ref f) = self.limits_file {
            if !f.exists() {
                return Err(crate::error::FastqcError::FileNotFound(f.clone()));
            }
        }

        Ok(())
    }

    /// Should we unzip the output?
    pub fn do_unzip(&self) -> bool {
        if self.extract {
            return true;
        }
        if self.noextract {
            return false;
        }
        false
    }
}
