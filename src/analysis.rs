use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use zip::write::SimpleFileOptions;
use zip::ZipWriter;

use crate::config::Config;
use crate::error::Result;
use crate::modules::{self, ModuleConfig, QcModule};
use crate::report::ReportArchive;
use crate::report::{html, json_report, summary, text_report};
use crate::sequence::bam::BamReader;
use crate::sequence::fastq::FastqReader;
use crate::sequence::Sequence;
use crate::sequence::SequenceFile;

/// Batch size for the pipeline channel.
/// 4096 sequences per batch balances throughput vs memory (~1.6 MB per batch for 150bp reads).
const BATCH_SIZE: usize = 4096;

/// Channel buffer size (number of batches). 16 batches ≈ 6 MB of buffered sequences.
const CHANNEL_BUFFER: usize = 16;

/// Process a single file (or file group) and generate reports.
pub fn process_file(path: &Path, config: &Config) -> Result<()> {
    let module_config = ModuleConfig::load(config);
    let mut modules = modules::create_module_list(config, &module_config);

    // Determine file type and create reader
    let format = detect_format(path, config);

    let file_name: String;

    match format.as_str() {
        "bam" | "sam" | "bam_mapped" | "sam_mapped" => {
            let only_mapped = format == "bam_mapped" || format == "sam_mapped";
            file_name = process_bam(path, config, only_mapped, &mut modules)?;
        }
        _ => {
            file_name = process_fastq_pipeline(path, config, &mut modules)?;
        }
    }

    if !config.quiet {
        eprintln!("Analysis complete for {}", file_name);
    }

    // Generate reports
    generate_reports(path, &file_name, &format, &mut modules, config)?;

    Ok(())
}

/// Process BAM/SAM files (sequential, as before).
fn process_bam(
    path: &Path,
    config: &Config,
    only_mapped: bool,
    modules: &mut Vec<Box<dyn QcModule>>,
) -> Result<String> {
    let mut reader = BamReader::new(path, config, only_mapped)?;
    let file_name = reader.name().to_string();
    let mut seq_count: u64 = 0;
    let mut last_percent: u8 = 0;

    loop {
        let pct = reader.percent_complete();
        match reader.next() {
            Some(Ok(seq)) => {
                seq_count += 1;
                for module in modules.iter_mut() {
                    if seq.is_filtered && module.ignore_filtered_sequences() {
                        continue;
                    }
                    module.process_sequence(&seq);
                }
                report_progress(config, &file_name, seq_count, pct, &mut last_percent);
            }
            Some(Err(e)) => return Err(e),
            None => break,
        }
    }

    Ok(file_name)
}

/// Number of worker threads for data-parallel module processing.
/// Each worker gets its own module set and processes a subset of sequences.
const NUM_WORKERS: usize = 4;

/// Process FASTQ files using data-parallel pipeline:
/// - 1 reader thread decompresses + parses FASTQ + tracks deduplication in file order
/// - N worker threads each process sequences through their own module copies
/// - After completion, merge all worker modules into the primary set
fn process_fastq_pipeline(
    path: &Path,
    config: &Config,
    modules: &mut Vec<Box<dyn QcModule>>,
) -> Result<String> {
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    if file_name == "stdin" || path.to_string_lossy() == "stdin" {
        return process_fastq_sequential(path, config, modules);
    }

    let module_config = ModuleConfig::load(config);
    let duplication_idx = modules
        .iter()
        .position(|m| m.name() == "Sequence Duplication Levels")
        .expect("duplication module must exist");
    let overrep_idx = modules
        .iter()
        .position(|m| m.name() == "Overrepresented sequences")
        .expect("overrepresented module must exist");

    let (reader_result, worker_results) = std::thread::scope(|s| {
        // Create per-worker channels
        let mut worker_txs = Vec::with_capacity(NUM_WORKERS);
        let mut worker_handles = Vec::with_capacity(NUM_WORKERS);

        for _ in 0..NUM_WORKERS {
            let (tx, rx) = mpsc::sync_channel::<Vec<Sequence>>(CHANNEL_BUFFER);
            worker_txs.push(tx);

            // Each worker creates its own module set and processes independently
            let worker_modules = modules::create_module_list(config, &module_config);
            let handle = s.spawn(move || -> Vec<Box<dyn QcModule>> {
                let mut mods = worker_modules;
                for batch in rx.iter() {
                    for seq in &batch {
                        for (module_idx, module) in mods.iter_mut().enumerate() {
                            // Duplication-related modules must observe the original
                            // file order, so the reader thread owns them.
                            if module_idx == duplication_idx || module_idx == overrep_idx {
                                continue;
                            }
                            if seq.is_filtered && module.ignore_filtered_sequences() {
                                continue;
                            }
                            module.process_sequence(seq);
                        }
                    }
                }
                mods
            });
            worker_handles.push(handle);
        }

        // Reader thread: decompress + parse + track deduplication sequentially.
        let overrep_module = &mut modules[overrep_idx];
        let overrep_ignores_filtered = overrep_module.ignore_filtered_sequences();
        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader = FastqReader::new(path, config)?;
            let mut batch = Vec::with_capacity(BATCH_SIZE);
            let mut worker_idx: usize = 0;

            loop {
                match reader.read_sequence() {
                    Ok(Some(seq)) => {
                        if !(seq.is_filtered && overrep_ignores_filtered) {
                            overrep_module.process_sequence(&seq);
                        }
                        batch.push(seq);
                        if batch.len() >= BATCH_SIZE {
                            if worker_txs[worker_idx].send(batch).is_err() {
                                break;
                            }
                            worker_idx = (worker_idx + 1) % NUM_WORKERS;
                            batch = Vec::with_capacity(BATCH_SIZE);
                        }
                    }
                    Ok(None) => {
                        if !batch.is_empty() {
                            let _ = worker_txs[worker_idx].send(batch);
                        }
                        break;
                    }
                    Err(e) => return Err(e),
                }
            }
            // Drop all senders to signal workers to finish
            drop(worker_txs);
            Ok(())
        });

        // Wait for reader
        let read_result = reader_handle.join().unwrap();

        // Collect worker results after the reader finishes and closes the channels.
        let mut worker_results = Vec::with_capacity(NUM_WORKERS);
        for handle in worker_handles {
            worker_results.push(handle.join().unwrap());
        }

        (read_result, worker_results)
    });

    reader_result?;

    // Merge worker-local modules back into the primary module set.
    let mut first = true;
    for worker_mods in worker_results {
        if first {
            for (i, wmod) in worker_mods.into_iter().enumerate() {
                if i == duplication_idx || i == overrep_idx {
                    continue;
                }
                modules[i] = wmod;
            }
            first = false;
        } else {
            for (i, wmod) in worker_mods.into_iter().enumerate() {
                if i == duplication_idx || i == overrep_idx {
                    continue;
                }
                modules[i].merge(wmod.into_any());
            }
        }
    }

    if !config.quiet {
        eprintln!("{}  processing complete", file_name);
    }

    Ok(file_name)
}

/// Sequential FASTQ processing fallback (for stdin).
fn process_fastq_sequential(
    path: &Path,
    config: &Config,
    modules: &mut Vec<Box<dyn QcModule>>,
) -> Result<String> {
    let mut reader = FastqReader::new(path, config)?;
    let file_name = reader.name().to_string();
    let mut seq_count: u64 = 0;
    let mut last_percent: u8 = 0;

    loop {
        let pct = reader.percent_complete();
        match reader.next() {
            Some(Ok(seq)) => {
                seq_count += 1;
                for module in modules.iter_mut() {
                    if seq.is_filtered && module.ignore_filtered_sequences() {
                        continue;
                    }
                    module.process_sequence(&seq);
                }
                report_progress(config, &file_name, seq_count, pct, &mut last_percent);
            }
            Some(Err(e)) => return Err(e),
            None => break,
        }
    }

    Ok(file_name)
}

fn detect_format(path: &Path, config: &Config) -> String {
    if let Some(ref fmt) = config.sequence_format {
        return fmt.clone();
    }

    let path_lower = path.to_string_lossy().to_lowercase();
    if path_lower.ends_with(".bam") || path_lower.ends_with(".ubam") {
        "bam".to_string()
    } else if path_lower.ends_with(".sam") {
        "sam".to_string()
    } else {
        "fastq".to_string()
    }
}

fn report_progress(config: &Config, name: &str, count: u64, percent: u8, last: &mut u8) {
    if config.quiet {
        return;
    }
    if count % 1000 != 0 {
        return;
    }
    let rounded = (percent / 5) * 5;
    if rounded > *last {
        *last = rounded;
        eprintln!("{}  {} sequences ({}%)", name, count, rounded);
    }
}

fn generate_reports(
    path: &Path,
    file_name: &str,
    detected_format: &str,
    modules: &mut Vec<Box<dyn QcModule>>,
    config: &Config,
) -> Result<()> {
    // Determine output paths
    let output_base = get_output_base(path, file_name, config);
    let html_path = PathBuf::from(format!("{}_fastqc.html", output_base));
    let zip_path = PathBuf::from(format!("{}_fastqc.zip", output_base));
    let json_path = PathBuf::from(format!("{}_fastqc.json", output_base));
    let folder_name = format!(
        "{}_fastqc",
        PathBuf::from(&output_base)
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
    );

    // Generate module reports
    let mut module_data: Vec<(String, String, Vec<(String, String)>)> = Vec::new();

    for module in modules.iter_mut() {
        let mut report = ReportArchive::new(file_name, config.clone());

        if !module.ignore_in_report() {
            module.make_report(&mut report)?;
        }

        module_data.push((report.data, report.html_body, report.images));
    }

    // Generate text report
    let data_strings: Vec<String> = module_data.iter().map(|(d, _, _)| d.clone()).collect();
    let text_report = text_report::generate_text_report(modules, &data_strings);

    // Generate summary
    let summary_text = summary::generate_summary(modules, file_name);

    // Generate HTML
    let html_content = html::generate_html(file_name, modules, &module_data);

    // Write HTML file
    let mut html_file = File::create(&html_path)?;
    html_file.write_all(html_content.as_bytes())?;

    // Write ZIP archive
    let zip_file = File::create(&zip_path)?;
    let mut zip = ZipWriter::new(zip_file);
    let options = SimpleFileOptions::default();

    // Create directories
    zip.add_directory(format!("{}/", folder_name), options)?;
    zip.add_directory(format!("{}/Icons/", folder_name), options)?;
    zip.add_directory(format!("{}/Images/", folder_name), options)?;

    // Write icons
    let icons = [
        (
            "fastqc_icon.png",
            include_bytes!("data/icons/fastqc_icon.png").as_slice(),
        ),
        (
            "warning.png",
            include_bytes!("data/icons/warning.png").as_slice(),
        ),
        (
            "error.png",
            include_bytes!("data/icons/error.png").as_slice(),
        ),
        ("tick.png", include_bytes!("data/icons/tick.png").as_slice()),
    ];
    for (name, data) in &icons {
        zip.start_file(format!("{}/Icons/{}", folder_name, name), options)?;
        zip.write_all(data)?;
    }

    // Write images
    for (_, _, images) in &module_data {
        for (img_name, svg_content) in images {
            let svg_name = img_name.replace(".png", ".svg");
            zip.start_file(format!("{}/Images/{}", folder_name, svg_name), options)?;
            zip.write_all(svg_content.as_bytes())?;
        }
    }

    // Write HTML report
    zip.start_file(format!("{}/fastqc_report.html", folder_name), options)?;
    zip.write_all(html_content.as_bytes())?;

    // Write data file
    zip.start_file(format!("{}/fastqc_data.txt", folder_name), options)?;
    zip.write_all(text_report.as_bytes())?;

    // Write summary
    zip.start_file(format!("{}/summary.txt", folder_name), options)?;
    zip.write_all(summary_text.as_bytes())?;

    zip.finish()?;

    // Optional extraction
    if config.do_unzip() {
        // Extract zip to same directory
        let zip_file = File::open(&zip_path)?;
        let mut archive = zip::ZipArchive::new(zip_file)?;
        let extract_dir = zip_path.parent().unwrap_or(Path::new("."));
        archive.extract(extract_dir)?;

        if config.delete {
            std::fs::remove_file(&zip_path)?;
        }
    }

    if !config.quiet {
        eprintln!("Report written: {}", html_path.display());
    }

    if config.json {
        let report = json_report::build_report(path, file_name, detected_format, modules, config);
        let mut json_file = File::create(&json_path)?;
        let pretty = serde_json::to_string_pretty(&report)?;
        json_file.write_all(pretty.as_bytes())?;
        if !config.quiet {
            eprintln!("JSON report written: {}", json_path.display());
        }
    }

    Ok(())
}

fn get_output_base(path: &Path, file_name: &str, config: &Config) -> String {
    let mut base = file_name.to_string();

    // Strip known extensions
    for ext in &[
        ".gz", ".bz2", ".txt", ".fastq", ".fq", ".csfastq", ".sam", ".bam", ".ubam",
    ] {
        if let Some(stripped) = base.strip_suffix(ext) {
            base = stripped.to_string();
        }
    }

    if let Some(ref dir) = config.output_dir {
        format!("{}/{}", dir.display(), base)
    } else if let Some(parent) = path.parent() {
        if parent.to_string_lossy().is_empty() {
            base
        } else {
            format!("{}/{}", parent.display(), base)
        }
    } else {
        base
    }
}
