#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::process;

use clap::Parser;
use rayon::prelude::*;

use fastqc_rs::analysis;
use fastqc_rs::config::Config;
use fastqc_rs::trim_galore::TrimGaloreConfig;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    match args.get(1).map(String::as_str) {
        Some("trim-galore") => run_trim_galore(),
        _ => run_qc(),
    }
}

fn run_trim_galore() {
    let config = TrimGaloreConfig::parse_from(
        std::iter::once("fastqc-rs trim-galore".to_string()).chain(std::env::args().skip(2)),
    );

    if let Err(e) = fastqc_rs::trim_galore::run(&config) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}

fn run_qc() {
    let config = Config::parse();

    if let Err(e) = config.validate() {
        eprintln!("Error: {}", e);
        process::exit(1);
    }

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .unwrap_or_else(|e| {
            eprintln!("Warning: could not configure thread pool: {}", e);
        });

    let files = config.files.clone();
    let had_errors = std::sync::atomic::AtomicBool::new(false);

    if config.threads > 1 && files.len() > 1 {
        // Process files in parallel
        files.par_iter().for_each(|file| {
            if !file.exists() {
                eprintln!("Error: file '{}' does not exist", file.display());
                had_errors.store(true, std::sync::atomic::Ordering::Relaxed);
                return;
            }
            if let Err(e) = analysis::process_file(file, &config) {
                eprintln!("Error processing '{}': {}", file.display(), e);
                had_errors.store(true, std::sync::atomic::Ordering::Relaxed);
            }
        });
    } else {
        // Process files sequentially
        for file in &files {
            if !file.exists() {
                eprintln!("Error: file '{}' does not exist", file.display());
                had_errors.store(true, std::sync::atomic::Ordering::Relaxed);
                continue;
            }
            if let Err(e) = analysis::process_file(file, &config) {
                eprintln!("Error processing '{}': {}", file.display(), e);
                had_errors.store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }
    }

    if had_errors.load(std::sync::atomic::Ordering::Relaxed) {
        process::exit(1);
    }
}
