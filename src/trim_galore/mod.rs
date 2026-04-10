pub mod adapter;
pub mod config;

pub use config::TrimGaloreConfig;

use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use config::TRIM_GALORE_VERSION;

/// Main entry point for the trim_galore module.
pub fn run(config: &TrimGaloreConfig) -> Result<(), Box<dyn std::error::Error>> {
    config
        .validate()
        .map_err(|e| format!("Configuration error: {}", e))?;

    // Create output directory if needed
    if let Some(ref dir) = config.output_dir {
        if !dir.exists() {
            fs::create_dir_all(dir)?;
            eprintln!("Created output directory: {}", dir.display());
        }
    }

    // Check cutadapt
    let cutadapt_path = config.cutadapt_path();
    let cutadapt_version = check_cutadapt(cutadapt_path)?;

    // --- Specialty trimming modes (exit after completion) ---

    if let Some(n) = config.hardtrim5 {
        eprintln!(
            "Hard-trimming from the 3'-end selected. Files will be trimmed to leave \
             the leftmost {} bp on the 5'-end.\n",
            n
        );
        for file in &config.files {
            hardtrim_5prime(file, n, config)?;
        }
        return Ok(());
    }

    if let Some(n) = config.hardtrim3 {
        eprintln!(
            "Hard-trimming from 5'-end selected. Files will be trimmed to leave \
             the rightmost {} bp on the 3'-end.\n",
            n
        );
        for file in &config.files {
            hardtrim_3prime(file, n, config)?;
        }
        return Ok(());
    }

    // --- Adapter detection ---

    let (adapter, adapter_name, report_msg) = determine_adapter(config)?;
    let adapter2 = determine_adapter2(config, &adapter);
    let cores_flag = build_cores_flag(&cutadapt_version, config.cores);

    // Effective clip_r2 for RRBS
    let mut effective_clip_r2 = config.clip_r2;
    if config.rrbs && config.paired && !config.non_directional && config.clip_r2.is_none() {
        eprintln!("Setting --clip_R2 2 (to remove methylation bias from the start of Read 2)");
        effective_clip_r2 = Some(2);
    }

    // Small RNA length cutoff adjustment
    let length_cutoff = if adapter == "TGGAATTCTCGG" && config.length == 20 {
        eprintln!(
            "Reducing length cutoff to 18bp for small RNA-Seq reads"
        );
        18
    } else {
        config.length
    };

    // --- Process files ---

    if config.paired {
        for pair in config.files.chunks(2) {
            trim_paired_end(
                &pair[0],
                &pair[1],
                &adapter,
                &adapter_name,
                &adapter2,
                &cutadapt_version,
                &cores_flag,
                report_msg.as_deref(),
                length_cutoff,
                effective_clip_r2,
                config,
            )?;
        }
    } else {
        for file in &config.files {
            trim_single_end(
                file,
                &adapter,
                &adapter_name,
                &cutadapt_version,
                &cores_flag,
                report_msg.as_deref(),
                length_cutoff,
                config,
            )?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Cutadapt helpers
// ---------------------------------------------------------------------------

fn check_cutadapt(path: &str) -> Result<String, Box<dyn std::error::Error>> {
    eprintln!(
        "Path to Cutadapt set as: '{}' ({})",
        path,
        if path == "cutadapt" {
            "default"
        } else {
            "user defined"
        }
    );

    let output = Command::new(path)
        .arg("--version")
        .output()
        .map_err(|_| {
            format!(
                "Failed to execute Cutadapt. Please install Cutadapt or specify \
                 --path_to_cutadapt /path/to/cutadapt"
            )
        })?;

    if !output.status.success() {
        return Err("Cutadapt --version returned non-zero exit code".into());
    }

    let version = String::from_utf8_lossy(&output.stdout).trim().to_string();
    eprintln!("Cutadapt seems to be working fine (tested command '{} --version')", path);
    eprintln!("Cutadapt version: {}", version);
    Ok(version)
}

fn build_cores_flag(cutadapt_version: &str, cores: usize) -> String {
    if cores <= 1 {
        return String::new();
    }

    // Parse major.minor from version string
    let nums: Vec<u32> = cutadapt_version
        .split('.')
        .take(2)
        .filter_map(|s| s.parse().ok())
        .collect();

    let (major, minor) = match nums.as_slice() {
        [m, s, ..] => (*m, *s),
        [m] => (*m, 0),
        _ => (1, 0),
    };

    if major == 1 && minor < 15 {
        eprintln!(
            "Warning: Cutadapt v{} does not support multi-core. Using single core.",
            cutadapt_version
        );
        return String::new();
    }

    eprintln!(
        "Cutadapt version {} supports multi-core. Setting -j {}",
        cutadapt_version, cores
    );
    format!("-j {}", cores)
}

// ---------------------------------------------------------------------------
// Adapter determination
// ---------------------------------------------------------------------------

fn determine_adapter(
    config: &TrimGaloreConfig,
) -> Result<(String, String, Option<String>), Box<dyn std::error::Error>> {
    if let Some(ref a) = config.adapter {
        return Ok((a.to_uppercase(), "user defined".to_string(), None));
    }
    if config.illumina {
        return Ok((
            "AGATCGGAAGAGC".into(),
            "Illumina TruSeq, Sanger iPCR; user defined".into(),
            None,
        ));
    }
    if config.nextera {
        return Ok((
            "CTGTCTCTTATA".into(),
            "Nextera Transposase sequence; user defined".into(),
            None,
        ));
    }
    if config.small_rna {
        return Ok((
            "TGGAATTCTCGG".into(),
            "Illumina small RNA adapter; user defined".into(),
            None,
        ));
    }
    if config.stranded_illumina {
        return Ok((
            "ACTGTCTCTTATA".into(),
            "Illumina stranded mRNA/Total RNA; user defined".into(),
            None,
        ));
    }
    if config.bgiseq {
        return Ok((
            "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA".into(),
            "BGISEQ/DNBSEQ/MGISEQ; user defined".into(),
            None,
        ));
    }

    // Auto-detect
    let r =
        adapter::autodetect_adapter(&config.files[0], config.consider_already_trimmed);
    Ok((r.sequence, r.name, r.report_message))
}

fn determine_adapter2(config: &TrimGaloreConfig, adapter: &str) -> String {
    if let Some(ref a2) = config.adapter2 {
        return a2.to_uppercase();
    }
    if config.paired {
        if adapter == "TGGAATTCTCGG" {
            eprintln!("Setting Illumina smallRNA 5' adapter as adapter 2: 'GATCGTCGGACT'");
            return "GATCGTCGGACT".into();
        }
        if adapter == "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" {
            let a2 = "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG";
            eprintln!("Setting BGISEQ Read 2 adapter: '{}'", a2);
            return a2.into();
        }
    }
    String::new()
}

// ---------------------------------------------------------------------------
// FASTQ I/O helpers
// ---------------------------------------------------------------------------

struct FastqRecord {
    id: String,
    seq: String,
    plus: String,
    qual: String,
}

fn read_record(reader: &mut dyn BufRead) -> Option<FastqRecord> {
    let mut id = String::new();
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();

    if reader.read_line(&mut id).ok()? == 0 {
        return None;
    }
    if reader.read_line(&mut seq).ok()? == 0 {
        return None;
    }
    if reader.read_line(&mut plus).ok()? == 0 {
        return None;
    }
    if reader.read_line(&mut qual).ok()? == 0 {
        return None;
    }

    fn chomp(s: &mut String) {
        while s.ends_with('\n') || s.ends_with('\r') {
            s.pop();
        }
    }
    chomp(&mut id);
    chomp(&mut seq);
    chomp(&mut plus);
    chomp(&mut qual);

    Some(FastqRecord { id, seq, plus, qual })
}

fn open_reader(path: &Path) -> Box<dyn BufRead> {
    let file =
        File::open(path).unwrap_or_else(|e| panic!("Cannot open '{}': {}", path.display(), e));
    if path.extension().map_or(false, |e| e == "gz") {
        Box::new(BufReader::with_capacity(
            256 * 1024,
            GzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(256 * 1024, file))
    }
}

fn create_writer(path: &Path, gzip: bool) -> Box<dyn Write> {
    let file =
        File::create(path).unwrap_or_else(|e| panic!("Cannot create '{}': {}", path.display(), e));
    if gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::fast())))
    } else {
        Box::new(BufWriter::new(file))
    }
}

fn count_ns(seq: &str) -> usize {
    seq.bytes().filter(|&b| b == b'N' || b == b'n').count()
}

fn is_gz(path: &Path) -> bool {
    path.extension().map_or(false, |e| e == "gz")
}

fn output_path(dir: &str, name: &str) -> PathBuf {
    if dir.is_empty() {
        PathBuf::from(name)
    } else {
        PathBuf::from(dir).join(name)
    }
}

/// Derive the trimmed output filename from the input filename.
fn trimmed_name(input: &Path) -> String {
    let name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    if name.ends_with(".fastq.gz") {
        name.replace(".fastq.gz", "_trimmed.fq")
    } else if name.ends_with(".fq.gz") {
        name.replace(".fq.gz", "_trimmed.fq")
    } else if name.ends_with(".fastq") {
        name.replace(".fastq", "_trimmed.fq")
    } else if name.ends_with(".fq") {
        name.replace(".fq", "_trimmed.fq")
    } else {
        format!("{}_trimmed.fq", name)
    }
}

/// Derive a hardtrim output filename.
fn hardtrim_name(input: &Path, n: usize, end: &str) -> String {
    let name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let suffix = format!(".{}bp_{}.fq", n, end);
    for ext in &[".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if name.ends_with(ext) {
            return name.replace(ext, &suffix);
        }
    }
    format!("{}{}", name, suffix)
}

fn write_record(w: &mut dyn Write, r: &FastqRecord) -> std::io::Result<()> {
    writeln!(w, "{}\n{}\n{}\n{}", r.id, r.seq, r.plus, r.qual)
}

/// Build cutadapt arguments for a single file.
fn build_cutadapt_args(
    adapter: &str,
    filename: &Path,
    cores_flag: &str,
    config: &TrimGaloreConfig,
) -> Vec<String> {
    let mut args: Vec<String> = Vec::new();

    // Extra cutadapt args
    if let Some(ref extra) = config.cutadapt_args {
        args.extend(extra.split_whitespace().map(String::from));
    }

    // Cores
    if !cores_flag.is_empty() {
        args.extend(cores_flag.split_whitespace().map(String::from));
    }

    args.push("-e".into());
    args.push(config.error_rate.to_string());

    // Quality / nextseq
    args.extend(config.quality_cutoff_args());

    args.push("-O".into());
    args.push(config.stringency.to_string());

    if config.trim_n {
        args.push("--trim-n".into());
    }

    args.push("-a".into());
    args.push(adapter.into());

    args.push(filename.to_string_lossy().into_owned());
    args
}

/// Spawn cutadapt, stream FASTQ from stdout.
/// Returns (BufReader over stdout, thread handle collecting stderr, Child).
fn spawn_cutadapt(
    cutadapt: &str,
    args: &[String],
) -> Result<
    (
        BufReader<std::process::ChildStdout>,
        std::thread::JoinHandle<String>,
        std::process::Child,
    ),
    Box<dyn std::error::Error>,
> {
    let mut child = Command::new(cutadapt)
        .args(args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to launch Cutadapt: {}", e))?;

    let stdout = child.stdout.take().unwrap();
    let stderr = child.stderr.take().unwrap();

    let stderr_handle = std::thread::spawn(move || {
        let reader = BufReader::new(stderr);
        let mut out = String::new();
        for line in reader.lines().map_while(Result::ok) {
            out.push_str(&line);
            out.push('\n');
        }
        out
    });

    Ok((BufReader::with_capacity(256 * 1024, stdout), stderr_handle, child))
}

/// Wait for cutadapt child and collect stderr.
fn finish_cutadapt(
    mut child: std::process::Child,
    stderr_handle: std::thread::JoinHandle<String>,
) -> Result<String, Box<dyn std::error::Error>> {
    let status = child.wait()?;
    let stderr_output = stderr_handle.join().unwrap_or_default();
    if !status.success() {
        return Err(format!(
            "Cutadapt terminated with exit signal: '{}'. Check error messages above.",
            status
        )
        .into());
    }
    Ok(stderr_output)
}

/// Write the report header common to both SE and PE modes.
fn report_header(
    filename: &Path,
    mode: &str,
    adapter: &str,
    adapter_name: &str,
    cutadapt_version: &str,
    length_cutoff: usize,
    config: &TrimGaloreConfig,
) -> String {
    format!(
        "\nSUMMARISING RUN PARAMETERS\n\
         ==========================\n\
         Input filename: {}\n\
         Trimming mode: {}\n\
         Trim Galore version: {}\n\
         Cutadapt version: {}\n\
         Quality Phred score cutoff: {}\n\
         Quality encoding type selected: ASCII+{}\n\
         Adapter sequence: '{}' ({})\n\
         Maximum trimming error rate: {}{}\n\
         Minimum required adapter overlap (stringency): {} bp\n\
         Minimum required sequence length before a sequence gets removed: {} bp\n",
        filename.display(),
        mode,
        TRIM_GALORE_VERSION,
        cutadapt_version,
        config.quality,
        config.phred_encoding(),
        adapter,
        adapter_name,
        config.error_rate,
        if (config.error_rate - 0.1).abs() < f64::EPSILON {
            " (default)"
        } else {
            ""
        },
        config.stringency,
        length_cutoff,
    )
}

fn write_report_file(dir: &str, name: &str, content: &str) {
    let path = output_path(dir, name);
    if let Err(e) = fs::write(&path, content) {
        eprintln!("Warning: could not write report '{}': {}", path.display(), e);
    } else {
        eprintln!("Writing report to '{}'", path.display());
    }
}

// ---------------------------------------------------------------------------
// Single-end trimming
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn trim_single_end(
    filename: &Path,
    adapter: &str,
    adapter_name: &str,
    cutadapt_version: &str,
    cores_flag: &str,
    report_msg: Option<&str>,
    length_cutoff: usize,
    config: &TrimGaloreConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let dir = config.output_dir_str();
    let gz = config.should_gzip(is_gz(filename));

    let mut out_name = trimmed_name(filename);
    let report_name = format!(
        "{}_trimming_report.txt",
        filename.file_name().unwrap().to_string_lossy()
    );

    let mut report =
        report_header(filename, "single-end", adapter, adapter_name, cutadapt_version, length_cutoff, config);
    if let Some(msg) = report_msg {
        report.push_str(msg);
    }
    eprintln!("{}", report);

    // Build and launch cutadapt
    let args = build_cutadapt_args(adapter, filename, cores_flag, config);
    eprintln!(
        "\n  >>> Now performing quality and adapter trimming for adapter '{}' from file {} <<<\n",
        adapter,
        filename.display()
    );

    let (mut reader, stderr_handle, child) = spawn_cutadapt(config.cutadapt_path(), &args)?;

    if gz {
        out_name.push_str(".gz");
    }
    let out_path = output_path(&dir, &out_name);
    eprintln!("Writing final adapter and quality trimmed output to {}\n", out_name);

    let mut writer = create_writer(&out_path, gz);

    let mut count = 0usize;
    let mut too_short = 0usize;
    let mut too_long = 0usize;
    let mut too_many_n = 0usize;

    while let Some(mut rec) = read_record(&mut reader) {
        count += 1;
        if count % 10_000_000 == 0 {
            eprintln!("{} sequences processed", count);
        }

        // 5' clip
        if let Some(clip) = config.clip_r1 {
            if rec.seq.len() > clip {
                if config.rename {
                    rec.id = format!("{}:clip5:{}", rec.id, &rec.seq[..clip]);
                }
                rec.seq = rec.seq[clip..].to_string();
                rec.qual = rec.qual[clip..].to_string();
            }
        }

        // 3' clip
        if let Some(clip) = config.three_prime_clip_r1 {
            if rec.seq.len() > clip {
                if config.rename {
                    rec.id = format!("{}:clip3:{}", rec.id, &rec.seq[rec.seq.len() - clip..]);
                }
                rec.seq.truncate(rec.seq.len() - clip);
                rec.qual.truncate(rec.qual.len() - clip);
            }
        }

        // max_n filter
        if let Some(max_n) = config.max_n {
            let nc = count_ns(&rec.seq);
            let exceeds = if max_n > 0.0 && max_n < 1.0 {
                !rec.seq.is_empty() && (nc as f64 / rec.seq.len() as f64) > max_n
            } else {
                nc > max_n as usize
            };
            if exceeds {
                too_many_n += 1;
                continue;
            }
        }

        // Length filters
        if rec.seq.len() < length_cutoff {
            too_short += 1;
            continue;
        }
        if let Some(ml) = config.max_length {
            if rec.seq.len() > ml {
                too_long += 1;
                continue;
            }
        }

        write_record(&mut writer, &rec)?;
    }

    writer.flush()?;
    drop(writer);

    let stderr_output = finish_cutadapt(child, stderr_handle)?;
    eprint!("{}", stderr_output);
    report.push('\n');
    report.push_str(&stderr_output);

    // Statistics
    let stats = format!(
        "\nRUN STATISTICS FOR INPUT FILE: {}\n{}\n{} sequences processed in total\n",
        filename.display(),
        "=".repeat(45),
        count,
    );
    eprintln!("{}", stats);
    report.push_str(&stats);

    let pct = |n: usize| -> String {
        if count > 0 {
            format!("{:.1}%", n as f64 / count as f64 * 100.0)
        } else {
            "N/A".into()
        }
    };

    let msg = format!(
        "Sequences removed because they became shorter than the length cutoff of {} bp:\t{} ({})\n",
        length_cutoff,
        too_short,
        pct(too_short)
    );
    eprintln!("{}", msg);
    report.push_str(&msg);

    if config.max_n.is_some() {
        let msg = format!(
            "Sequences removed because of N content:\t{} ({})\n",
            too_many_n,
            pct(too_many_n)
        );
        eprintln!("{}", msg);
        report.push_str(&msg);
    }
    if let Some(ml) = config.max_length {
        let msg = format!(
            "Sequences removed because they were longer than {} bp:\t{} ({})\n",
            ml,
            too_long,
            pct(too_long)
        );
        eprintln!("{}", msg);
        report.push_str(&msg);
    }

    // Basename rename
    if let Some(ref bn) = config.basename {
        if let Some(pos) = out_name.find("_trimmed") {
            let new_name = format!("{}{}", bn, &out_name[pos..]);
            let old = output_path(&dir, &out_name);
            let new = output_path(&dir, &new_name);
            fs::rename(&old, &new)?;
            out_name = new_name;
        }
    }

    if !config.no_report_file {
        write_report_file(&dir, &report_name, &report);
    }

    // FastQC
    if config.fastqc {
        run_fastqc(&output_path(&dir, &out_name), config);
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Paired-end trimming
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn trim_paired_end(
    file1: &Path,
    file2: &Path,
    adapter: &str,
    adapter_name: &str,
    adapter2: &str,
    cutadapt_version: &str,
    cores_flag: &str,
    report_msg: Option<&str>,
    length_cutoff: usize,
    effective_clip_r2: Option<usize>,
    config: &TrimGaloreConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let r2_adapter = if adapter2.is_empty() {
        adapter
    } else {
        adapter2
    };
    let r2_name = if adapter2.is_empty() {
        adapter_name
    } else {
        "user defined (read 2)"
    };

    // Trim R1
    let trimmed_1 = trim_file_for_paired(
        file1,
        adapter,
        adapter_name,
        cutadapt_version,
        cores_flag,
        report_msg,
        length_cutoff,
        config,
    )?;

    // Trim R2
    let trimmed_2 = trim_file_for_paired(
        file2,
        r2_adapter,
        r2_name,
        cutadapt_version,
        cores_flag,
        report_msg,
        length_cutoff,
        config,
    )?;

    // Validate
    eprintln!(
        "\n>>>>> Now validating the length of the 2 paired-end infiles: {} and {} <<<<<\n",
        trimmed_1.display(),
        trimmed_2.display()
    );

    let gz = config.should_gzip(is_gz(file1));
    validate_paired_end(
        &trimmed_1,
        &trimmed_2,
        gz,
        length_cutoff,
        effective_clip_r2,
        config,
    )?;

    // Clean up intermediates
    eprintln!(
        "Deleting both intermediate output files {} and {}",
        trimmed_1.display(),
        trimmed_2.display()
    );
    let _ = fs::remove_file(&trimmed_1);
    let _ = fs::remove_file(&trimmed_2);

    eprintln!("\n{}\n", "=".repeat(100));
    Ok(())
}

/// Trim a single file with cutadapt for paired-end use (no length filtering).
#[allow(clippy::too_many_arguments)]
fn trim_file_for_paired(
    filename: &Path,
    adapter: &str,
    adapter_name: &str,
    cutadapt_version: &str,
    cores_flag: &str,
    report_msg: Option<&str>,
    length_cutoff: usize,
    config: &TrimGaloreConfig,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let dir = config.output_dir_str();
    let gz = config.should_gzip(is_gz(filename));

    let mut out_name = trimmed_name(filename);
    let report_name = format!(
        "{}_trimming_report.txt",
        filename.file_name().unwrap().to_string_lossy()
    );

    let mut report =
        report_header(filename, "paired-end", adapter, adapter_name, cutadapt_version, length_cutoff, config);
    if let Some(msg) = report_msg {
        report.push_str(msg);
    }
    eprintln!("{}", report);

    let args = build_cutadapt_args(adapter, filename, cores_flag, config);
    eprintln!(
        "\n  >>> Now performing quality and adapter trimming for adapter '{}' from file {} <<<\n",
        adapter,
        filename.display()
    );

    let (mut reader, stderr_handle, child) = spawn_cutadapt(config.cutadapt_path(), &args)?;

    if gz {
        out_name.push_str(".gz");
    }
    let out_path = output_path(&dir, &out_name);
    eprintln!("Writing trimmed output to {}\n", out_name);

    let mut writer = create_writer(&out_path, gz);
    let mut count = 0usize;

    // In paired-end mode we skip per-file length filtering (done in validation)
    while let Some(rec) = read_record(&mut reader) {
        count += 1;
        if count % 10_000_000 == 0 {
            eprintln!("{} sequences processed", count);
        }
        write_record(&mut writer, &rec)?;
    }

    writer.flush()?;
    drop(writer);

    let stderr_output = finish_cutadapt(child, stderr_handle)?;
    eprint!("{}", stderr_output);
    report.push('\n');
    report.push_str(&stderr_output);

    let stats = format!(
        "\nRUN STATISTICS FOR INPUT FILE: {}\n{}\n{} sequences processed in total\n\
         The length threshold of paired-end sequences gets evaluated later on (in the validation step)\n",
        filename.display(),
        "=".repeat(45),
        count,
    );
    eprintln!("{}", stats);
    report.push_str(&stats);

    if !config.no_report_file {
        write_report_file(&dir, &report_name, &report);
    }

    Ok(out_path)
}

/// Validate that both reads of each pair meet length requirements.
fn validate_paired_end(
    file1: &Path,
    file2: &Path,
    gz: bool,
    length_cutoff: usize,
    effective_clip_r2: Option<usize>,
    config: &TrimGaloreConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let dir = config.output_dir_str();

    let mut r1 = open_reader(file1);
    let mut r2 = open_reader(file2);

    let name1 = file1
        .file_name()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let name2 = file2
        .file_name()
        .unwrap()
        .to_string_lossy()
        .to_string();

    let strip_trimmed = |s: &str| -> String {
        s.replace("_trimmed.fq.gz", "")
            .replace("_trimmed.fq", "")
    };

    let mut o1 = format!("{}_val_1.fq", strip_trimmed(&name1));
    let mut o2 = format!("{}_val_2.fq", strip_trimmed(&name2));

    if gz {
        o1.push_str(".gz");
        o2.push_str(".gz");
    }

    eprintln!("Writing validated paired-end Read 1 reads to {}", o1);
    eprintln!("Writing validated paired-end Read 2 reads to {}\n", o2);

    let mut w1 = create_writer(&output_path(&dir, &o1), gz);
    let mut w2 = create_writer(&output_path(&dir, &o2), gz);

    let (mut uw1, mut uw2) = if config.retain_unpaired {
        let u1 = format!("{}_unpaired_1.fq{}", strip_trimmed(&name1), if gz { ".gz" } else { "" });
        let u2 = format!("{}_unpaired_2.fq{}", strip_trimmed(&name2), if gz { ".gz" } else { "" });
        eprintln!("Writing unpaired read 1 reads to {}", u1);
        eprintln!("Writing unpaired read 2 reads to {}\n", u2);
        (
            Some(create_writer(&output_path(&dir, &u1), gz)),
            Some(create_writer(&output_path(&dir, &u2), gz)),
        )
    } else {
        (None, None)
    };

    let mut count = 0usize;
    let mut pairs_removed = 0usize;
    let mut n_pairs = 0usize;
    let mut long_pairs = 0usize;
    let mut unpaired_r1 = 0usize;
    let mut unpaired_r2 = 0usize;

    loop {
        let rec1 = read_record(&mut r1);
        let rec2 = read_record(&mut r2);

        match (rec1, rec2) {
            (Some(mut a), Some(mut b)) => {
                count += 1;

                // Clip R1 5'
                if let Some(c) = config.clip_r1 {
                    clip_5prime(&mut a, c, config.rename);
                }
                // Clip R2 5'
                if let Some(c) = effective_clip_r2 {
                    clip_5prime(&mut b, c, config.rename);
                }
                // Clip R1 3'
                if let Some(c) = config.three_prime_clip_r1 {
                    clip_3prime(&mut a, c, config.rename);
                }
                // Clip R2 3'
                if let Some(c) = config.three_prime_clip_r2 {
                    clip_3prime(&mut b, c, config.rename);
                }

                // N filter
                if let Some(max_n) = config.max_n {
                    if exceeds_n(&a.seq, max_n) || exceeds_n(&b.seq, max_n) {
                        n_pairs += 1;
                        pairs_removed += 1;
                        continue;
                    }
                }

                // Length filter
                if a.seq.len() < length_cutoff || b.seq.len() < length_cutoff {
                    pairs_removed += 1;
                    if config.retain_unpaired {
                        if a.seq.len() >= config.length_1 {
                            if let Some(ref mut w) = uw1 {
                                write_record(w, &a)?;
                            }
                            unpaired_r1 += 1;
                        }
                        if b.seq.len() >= config.length_2 {
                            if let Some(ref mut w) = uw2 {
                                write_record(w, &b)?;
                            }
                            unpaired_r2 += 1;
                        }
                    }
                    continue;
                }

                // Max length filter
                if let Some(ml) = config.max_length {
                    if a.seq.len() > ml || b.seq.len() > ml {
                        long_pairs += 1;
                        continue;
                    }
                }

                write_record(&mut w1, &a)?;
                write_record(&mut w2, &b)?;
            }
            (None, None) => break,
            (Some(_), None) => {
                return Err(format!(
                    "Read 2 output is truncated at sequence {}",
                    count
                )
                .into());
            }
            (None, Some(_)) => {
                return Err(format!(
                    "Read 1 output is truncated at sequence {}",
                    count
                )
                .into());
            }
        }
    }

    w1.flush()?;
    w2.flush()?;
    if let Some(ref mut w) = uw1 {
        w.flush()?;
    }
    if let Some(ref mut w) = uw2 {
        w.flush()?;
    }

    let pct = |n: usize| -> String {
        if count > 0 {
            format!("{:.2}%", n as f64 / count as f64 * 100.0)
        } else {
            "N/A".into()
        }
    };

    eprintln!("Total number of sequences analysed: {}\n", count);
    eprintln!(
        "Number of sequence pairs removed because at least one read was shorter than the length cutoff ({} bp): {} ({})",
        length_cutoff, pairs_removed, pct(pairs_removed)
    );
    if config.max_n.is_some() {
        eprintln!(
            "Number of sequence pairs removed because of N content: {} ({})",
            n_pairs,
            pct(n_pairs)
        );
    }
    if config.max_length.is_some() {
        eprintln!(
            "Sequence pairs removed because of max length: {} ({})",
            long_pairs,
            pct(long_pairs)
        );
    }
    if config.retain_unpaired {
        eprintln!("Number of unpaired read 1 reads printed: {}", unpaired_r1);
        eprintln!("Number of unpaired read 2 reads printed: {}", unpaired_r2);
    }

    // FastQC on validated output
    if config.fastqc {
        eprintln!("\n  >>> Now running FastQC on the validated data <<<\n");
        run_fastqc(&output_path(&dir, &o1), config);
        run_fastqc(&output_path(&dir, &o2), config);
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Hardtrim
// ---------------------------------------------------------------------------

fn hardtrim_5prime(
    filename: &Path,
    n: usize,
    config: &TrimGaloreConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    eprintln!("Input file: {}", filename.display());
    let dir = config.output_dir_str();
    let gz = config.should_gzip(is_gz(filename));

    let mut out_name = hardtrim_name(filename, n, "5prime");
    if gz {
        out_name.push_str(".gz");
    }

    let out = output_path(&dir, &out_name);
    eprintln!(
        "Writing trimmed version (first {} bp) to '{}'",
        n, out_name
    );

    let mut reader = open_reader(filename);
    let mut writer = create_writer(&out, gz);
    let mut count = 0usize;

    while let Some(mut rec) = read_record(&mut reader) {
        count += 1;
        if rec.seq.len() > n {
            if config.rename {
                rec.id = format!("{}:clip5:{}", rec.id, &rec.seq[n..]);
            }
            rec.seq.truncate(n);
            rec.qual.truncate(n);
        }
        write_record(&mut writer, &rec)?;
    }

    writer.flush()?;
    eprintln!(
        "\nFinished writing {} sequences from {}\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
        count,
        filename.display()
    );
    Ok(())
}

fn hardtrim_3prime(
    filename: &Path,
    n: usize,
    config: &TrimGaloreConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    eprintln!("Input file: {}", filename.display());
    let dir = config.output_dir_str();
    let gz = config.should_gzip(is_gz(filename));

    let mut out_name = hardtrim_name(filename, n, "3prime");
    if gz {
        out_name.push_str(".gz");
    }

    let out = output_path(&dir, &out_name);
    eprintln!(
        "Writing trimmed version (last {} bp) to '{}'",
        n, out_name
    );

    let mut reader = open_reader(filename);
    let mut writer = create_writer(&out, gz);
    let mut count = 0usize;

    while let Some(mut rec) = read_record(&mut reader) {
        count += 1;
        if rec.seq.len() > n {
            let start = rec.seq.len() - n;
            if config.rename {
                rec.id = format!("{}:clip3:{}", rec.id, &rec.seq[..start]);
            }
            rec.seq = rec.seq[start..].to_string();
            rec.qual = rec.qual[start..].to_string();
        }
        write_record(&mut writer, &rec)?;
    }

    writer.flush()?;
    eprintln!(
        "\nFinished writing {} sequences from {}\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
        count,
        filename.display()
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

fn clip_5prime(rec: &mut FastqRecord, n: usize, rename: bool) {
    if rec.seq.len() > n {
        if rename {
            rec.id = format!("{}:clip5:{}", rec.id, &rec.seq[..n]);
        }
        rec.seq = rec.seq[n..].to_string();
        rec.qual = rec.qual[n..].to_string();
    }
}

fn clip_3prime(rec: &mut FastqRecord, n: usize, rename: bool) {
    if rec.seq.len() > n {
        let end = rec.seq.len() - n;
        if rename {
            rec.id = format!("{}:clip3:{}", rec.id, &rec.seq[end..]);
        }
        rec.seq.truncate(end);
        rec.qual.truncate(end);
    }
}

fn exceeds_n(seq: &str, max_n: f64) -> bool {
    let nc = count_ns(seq);
    if max_n > 0.0 && max_n < 1.0 {
        !seq.is_empty() && (nc as f64 / seq.len() as f64) > max_n
    } else {
        nc > max_n as usize
    }
}

fn run_fastqc(path: &Path, config: &TrimGaloreConfig) {
    eprintln!("  >>> Now running FastQC on {} <<<", path.display());
    let mut cmd = Command::new("fastqc");
    if let Some(ref args) = config.fastqc_args {
        for a in args.split_whitespace() {
            cmd.arg(a);
        }
    }
    cmd.arg(path.to_string_lossy().as_ref());
    match cmd.status() {
        Ok(s) if !s.success() => {
            eprintln!("Warning: FastQC exited with status {}", s);
        }
        Err(e) => {
            eprintln!("Warning: could not run FastQC: {}", e);
        }
        _ => {}
    }
}
