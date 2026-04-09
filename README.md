# FastQC-RS

A Rust implementation of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), a quality control tool for high throughput sequence data. This is a 1:1 rewrite of FastQC v0.12.1 with identical output format and analysis algorithms.

## Features

- **12 analysis modules** with identical algorithms and pass/warn/fail thresholds:
  1. Basic Statistics
  2. Per Base Sequence Quality
  3. Per Tile Sequence Quality
  4. Per Sequence Quality Scores
  5. Per Base Sequence Content
  6. Per Sequence GC Content
  7. Per Base N Content
  8. Sequence Length Distribution
  9. Sequence Duplication Levels
  10. Overrepresented Sequences
  11. Adapter Content
  12. Kmer Content

- **Input formats**: FASTQ (plain, gzip, bzip2), BAM, SAM
- **Output**: HTML report with SVG graphs, ZIP archive, `fastqc_data.txt`, `summary.txt`
- **Output compatibility**: Text reports (`fastqc_data.txt`) are byte-identical to Java FastQC output
- **Multi-file parallel processing** via rayon
- **Single binary** with embedded configuration files

## Installation

```bash
cargo install --path .
```

Or build from source:

```bash
cargo build --release
```

## Usage

```bash
# Basic usage
fastqc-rs input.fastq

# Multiple files with parallel processing
fastqc-rs -t 4 sample1.fastq.gz sample2.fastq.gz

# Specify output directory
fastqc-rs -o results/ input.fastq

# BAM/SAM files
fastqc-rs input.bam
fastqc-rs -f sam_mapped input.sam   # Only mapped reads, with soft-clip removal

# Extract results from ZIP
fastqc-rs --extract input.fastq

# Quiet mode (suppress progress)
fastqc-rs -q input.fastq
```

## CLI Options

| Option | Description |
|--------|-------------|
| `-o, --outdir <DIR>` | Output directory (must exist) |
| `--extract` | Unzip output after creation |
| `--noextract` | Don't unzip output (default) |
| `--delete` | Delete ZIP after extraction |
| `-f, --format <FMT>` | Force format: `fastq`, `bam`, `sam`, `bam_mapped`, `sam_mapped` |
| `-c, --contaminants <FILE>` | Custom contaminant list |
| `-a, --adapters <FILE>` | Custom adapter list |
| `-l, --limits <FILE>` | Custom pass/warn/fail thresholds |
| `-t, --threads <N>` | Number of files to process simultaneously (default: 1) |
| `-k, --kmers <SIZE>` | Kmer length (default: 7) |
| `-q, --quiet` | Suppress progress messages |
| `--casava` | CASAVA mode (filter flagged reads) |
| `--nogroup` | Disable base position grouping |
| `--expgroup` | Use exponential base grouping |
| `--min-length <BP>` | Minimum sequence length for grouping |
| `--dup-length <BP>` | Truncation length for duplication analysis |
| `--svg` | Output SVG graphs |

## Output

For each input file `sample.fastq.gz`, produces:
- `sample_fastqc.html` — Interactive HTML report
- `sample_fastqc.zip` — Archive containing:
  - `fastqc_report.html`
  - `fastqc_data.txt` — Tab-delimited analysis data
  - `summary.txt` — PASS/WARN/FAIL per module
  - `Icons/` — Status icons
  - `Images/` — SVG graphs

## Compatibility

Output format is fully compatible with Java FastQC v0.12.1:
- `fastqc_data.txt` is identical (verified against golden test files)
- `summary.txt` uses the same PASS/WARN/FAIL format
- Same module ordering and threshold logic
- Same CLI flags (with minor naming convention differences for multi-word flags)

## License

GPL-3.0 (same as the original FastQC)

## Acknowledgements

Based on [FastQC](https://github.com/s-andrews/FastQC) by Simon Andrews at the Babraham Institute.
