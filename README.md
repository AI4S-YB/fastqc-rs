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
- **Output compatibility**: Text reports match Java FastQC output (identical PASS/WARN/FAIL, near-identical data values)
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

## Performance

Benchmarked on a paired-end Illumina dataset (~1.15 GB / ~1.20 GB gzipped FASTQ, ~9.9M reads x 150bp):

### Baseline (direct Rust rewrite, no optimization)

| File | FastQC v0.12.1 (Java) | fastqc-rs (Rust) | Speedup |
|------|----------------------|------------------|---------|
| SPL1E1_raw_1.fastq.gz (1.15 GB) | 48.6s | 46.8s | 1.04x |
| SPL1E1_raw_2.fastq.gz (1.20 GB) | 47.6s | 45.7s | 1.04x |

### Optimized v1 (zlib-rs + 2-thread pipeline + ahash + LTO)

Optimizations applied: zlib-rs decompression backend, reader/processor pipeline (overlapping I/O with compute), AHashMap for hot-path modules, in-place ASCII uppercase, 256KB I/O buffer, LTO + codegen-units=1.

| File | FastQC v0.12.1 (Java) | fastqc-rs (optimized) | Speedup vs Java | Speedup vs baseline |
|------|----------------------|----------------------|-----------------|---------------------|
| SPL1E1_raw_1.fastq.gz (1.15 GB) | 48.6s | 38.8s | 1.25x | 1.21x |
| SPL1E1_raw_2.fastq.gz (1.20 GB) | 47.6s | 39.1s | 1.22x | 1.17x |

> Tested on Linux 6.6 (WSL2). Pipeline uses 2 threads (reader + processor). CPU utilization ~117%. Bottleneck is module processing (~85% of wall time).

### Optimized v2 (data-parallel multi-threaded processing)

Added data-parallel architecture: 1 reader thread + 4 worker threads, each with independent module copies. Workers process sequence subsets in parallel, results merged after completion.

| File | FastQC v0.12.1 (Java) | fastqc-rs (multi-threaded) | Speedup vs Java | Speedup vs baseline |
|------|----------------------|---------------------------|-----------------|---------------------|
| SPL1E1_raw_1.fastq.gz (1.15 GB) | 48.6s | 9.4s | **5.2x** | **5.0x** |
| SPL1E1_raw_2.fastq.gz (1.20 GB) | 47.6s | 9.4s | **5.1x** | **4.9x** |

> Tested on Linux 6.6 (WSL2). Uses 5 threads total (1 reader + 4 workers). CPU utilization ~500%. All PASS/WARN/FAIL assessments remain correct.

## Compatibility

Output format is compatible with Java FastQC v0.12.1:
- `summary.txt` — identical (PASS/WARN/FAIL per module match exactly)
- `fastqc_data.txt` — nearly identical with minor differences:
  - **Per sequence quality scores**: Rust uses proper rounding (`round()`) for per-read mean quality, while Java uses integer division (truncation). This shifts some reads to adjacent quality bins but does not affect the overall PASS/WARN/FAIL assessment.
  - **Number formatting**: Rust outputs small values in decimal notation (e.g., `0.000123`) while Java uses scientific notation (e.g., `1.23E-4`). Values are numerically equivalent.
  - **Floating-point precision**: Last 1-2 digits may differ due to f64 vs Java double rounding.
- Modules with identical data: Basic Statistics, Per Base Sequence Quality, Per Sequence GC Content, Sequence Length Distribution, Overrepresented Sequences
- Same module ordering and threshold logic
- Same CLI flags (with minor naming convention differences for multi-word flags)

## License

GPL-3.0 (same as the original FastQC)

## Acknowledgements

Based on [FastQC](https://github.com/s-andrews/FastQC) by Simon Andrews at the Babraham Institute.
