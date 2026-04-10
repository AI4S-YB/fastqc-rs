use std::path::PathBuf;

use clap::Parser;

pub const TRIM_GALORE_VERSION: &str = "0.1.0";

/// Quality/adapter/RRBS trimming powered by Cutadapt
#[derive(Parser, Debug, Clone)]
#[command(name = "trim-galore", version = TRIM_GALORE_VERSION)]
pub struct TrimGaloreConfig {
    /// Input FASTQ files
    #[arg(required = true)]
    pub files: Vec<PathBuf>,

    /// Quality Phred score cutoff
    #[arg(short = 'q', long = "quality", default_value = "20")]
    pub quality: u32,

    /// Adapter sequence to trim
    #[arg(short = 'a', long = "adapter")]
    pub adapter: Option<String>,

    /// Optional adapter for read 2 (paired-end only)
    #[arg(long = "adapter2", alias = "a2")]
    pub adapter2: Option<String>,

    /// Illumina universal adapter (AGATCGGAAGAGC)
    #[arg(long)]
    pub illumina: bool,

    /// Nextera transposase adapter (CTGTCTCTTATA)
    #[arg(long)]
    pub nextera: bool,

    /// Illumina small RNA adapter (TGGAATTCTCGG)
    #[arg(long)]
    pub small_rna: bool,

    /// Illumina stranded mRNA adapter (ACTGTCTCTTATA)
    #[arg(long)]
    pub stranded_illumina: bool,

    /// BGISEQ/DNBSEQ/MGISEQ adapter sequences
    #[arg(long)]
    pub bgiseq: bool,

    /// Minimum overlap with adapter required to trim (default: 1)
    #[arg(long, default_value = "1")]
    pub stringency: u32,

    /// Maximum allowed error rate (default: 0.1)
    #[arg(short = 'e', long = "error_rate", default_value = "0.1")]
    pub error_rate: f64,

    /// Phred+33 quality encoding (default)
    #[arg(long)]
    pub phred33: bool,

    /// Phred+64 quality encoding
    #[arg(long)]
    pub phred64: bool,

    /// Run FastQC after trimming
    #[arg(long)]
    pub fastqc: bool,

    /// Extra FastQC arguments (e.g. "--nogroup --outdir /home/")
    #[arg(long)]
    pub fastqc_args: Option<String>,

    /// Paired-end trimming mode
    #[arg(long, alias = "paired_end")]
    pub paired: bool,

    /// Minimum read length after trimming (default: 20)
    #[arg(long, default_value = "20")]
    pub length: usize,

    /// Maximum read length (for small RNA)
    #[arg(long)]
    pub max_length: Option<usize>,

    /// Maximum number of Ns allowed (or fraction if 0 < value < 1)
    #[arg(long = "max_n")]
    pub max_n: Option<f64>,

    /// Remove Ns from both ends of the read
    #[arg(long = "trim-n")]
    pub trim_n: bool,

    /// Compress output with GZIP
    #[arg(long)]
    pub gzip: bool,

    /// Don't compress output with GZIP
    #[arg(long)]
    pub dont_gzip: bool,

    /// Output directory
    #[arg(short = 'o', long = "output_dir")]
    pub output_dir: Option<PathBuf>,

    /// Path to cutadapt executable
    #[arg(long)]
    pub path_to_cutadapt: Option<String>,

    /// Number of cores for trimming
    #[arg(short = 'j', long = "cores", default_value = "1")]
    pub cores: usize,

    /// Clip N bp from 5' end of read 1
    #[arg(long = "clip_R1")]
    pub clip_r1: Option<usize>,

    /// Clip N bp from 5' end of read 2
    #[arg(long = "clip_R2")]
    pub clip_r2: Option<usize>,

    /// Clip N bp from 3' end of read 1
    #[arg(long = "three_prime_clip_R1")]
    pub three_prime_clip_r1: Option<usize>,

    /// Clip N bp from 3' end of read 2
    #[arg(long = "three_prime_clip_R2")]
    pub three_prime_clip_r2: Option<usize>,

    /// RRBS mode (MspI-digested)
    #[arg(long = "rrbs")]
    pub rrbs: bool,

    /// Non-directional RRBS
    #[arg(long)]
    pub non_directional: bool,

    /// Keep quality-trimmed intermediate file (RRBS only)
    #[arg(long)]
    pub keep: bool,

    /// Retain unpaired reads when one mate is too short
    #[arg(long)]
    pub retain_unpaired: bool,

    /// Unpaired read 1 length cutoff (default: 35)
    #[arg(long = "length_1", alias = "r1", default_value = "35")]
    pub length_1: usize,

    /// Unpaired read 2 length cutoff (default: 35)
    #[arg(long = "length_2", alias = "r2", default_value = "35")]
    pub length_2: usize,

    /// Don't generate a report file
    #[arg(long)]
    pub no_report_file: bool,

    /// Suppress all warnings
    #[arg(long)]
    pub suppress_warn: bool,

    /// NextSeq/NovaSeq 2-colour quality trimming cutoff
    #[arg(long = "nextseq", alias = "2colour")]
    pub nextseq: Option<u32>,

    /// Use this basename for output files instead of deriving from input
    #[arg(long)]
    pub basename: Option<String>,

    /// Consider file already trimmed if adapter count <= threshold
    #[arg(long)]
    pub consider_already_trimmed: Option<usize>,

    /// Hard-trim sequences to N bp from the 5' end
    #[arg(long)]
    pub hardtrim5: Option<usize>,

    /// Hard-trim sequences to N bp from the 3' end
    #[arg(long)]
    pub hardtrim3: Option<usize>,

    /// Poly-A tail trimming mode (experimental)
    #[arg(long = "polyA")]
    pub poly_a: bool,

    /// Add clipped sequences to read IDs
    #[arg(long)]
    pub rename: bool,

    /// Extra arguments passed directly to cutadapt
    #[arg(long)]
    pub cutadapt_args: Option<String>,
}

impl TrimGaloreConfig {
    pub fn validate(&self) -> Result<(), String> {
        if self.files.is_empty() {
            return Err("No input files specified".into());
        }
        for f in &self.files {
            if !f.exists() {
                return Err(format!("Input file '{}' does not exist", f.display()));
            }
        }

        if self.paired && self.files.len() % 2 != 0 {
            return Err(
                "Paired-end mode requires an even number of input files".into(),
            );
        }

        if self.phred33 && self.phred64 {
            return Err(
                "Specify only one quality encoding (--phred33 or --phred64)".into(),
            );
        }

        let presets = [
            self.illumina,
            self.nextera,
            self.small_rna,
            self.stranded_illumina,
            self.bgiseq,
        ];
        if presets.iter().filter(|&&x| x).count() > 1 {
            return Err(
                "Cannot use multiple adapter presets simultaneously".into(),
            );
        }
        if self.adapter.is_some() && presets.iter().any(|&x| x) {
            return Err(
                "Cannot specify both a custom adapter and an adapter preset".into(),
            );
        }

        if !(0.0..=1.0).contains(&self.error_rate) {
            return Err("Error rate must be between 0 and 1".into());
        }

        if self.nextseq.is_some() && self.quality != 20 {
            return Err("--quality and --nextseq are mutually exclusive".into());
        }

        if self.non_directional && !self.rrbs {
            return Err("--non_directional requires --rrbs".into());
        }

        if self.adapter2.is_some() && !self.paired {
            return Err("--adapter2 requires --paired".into());
        }

        if let Some(n) = self.hardtrim5 {
            if n == 0 || n >= 1000 {
                return Err("--hardtrim5 must be between 1 and 999".into());
            }
        }
        if let Some(n) = self.hardtrim3 {
            if n == 0 || n >= 1000 {
                return Err("--hardtrim3 must be between 1 and 999".into());
            }
        }

        if let Some(ref dir) = self.output_dir {
            if dir.exists() && !dir.is_dir() {
                return Err(format!("'{}' is not a directory", dir.display()));
            }
        }

        Ok(())
    }

    pub fn phred_encoding(&self) -> u32 {
        if self.phred64 {
            64
        } else {
            33
        }
    }

    pub fn cutadapt_path(&self) -> &str {
        self.path_to_cutadapt.as_deref().unwrap_or("cutadapt")
    }

    pub fn output_dir_str(&self) -> String {
        match &self.output_dir {
            Some(dir) => {
                let s = dir.to_string_lossy().to_string();
                if s.ends_with('/') {
                    s
                } else {
                    format!("{}/", s)
                }
            }
            None => String::new(),
        }
    }

    pub fn should_gzip(&self, input_is_gz: bool) -> bool {
        if self.dont_gzip {
            return false;
        }
        self.gzip || input_is_gz
    }

    /// Build the quality cutoff argument(s) for cutadapt
    pub fn quality_cutoff_args(&self) -> Vec<String> {
        if let Some(nextseq) = self.nextseq {
            vec![format!("--nextseq-trim={}", nextseq)]
        } else {
            let cutoff = if self.phred64 {
                self.quality + 31
            } else {
                self.quality
            };
            vec!["-q".to_string(), cutoff.to_string()]
        }
    }
}
