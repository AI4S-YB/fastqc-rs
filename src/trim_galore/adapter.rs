use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

pub struct AdapterResult {
    pub sequence: String,
    pub name: String,
    pub report_message: Option<String>,
}

struct AdapterCandidate {
    name: &'static str,
    display_name: &'static str,
    sequence: &'static str,
    count: usize,
}

/// Scan the first 1 million reads to auto-detect the most common adapter type.
pub fn autodetect_adapter(
    file: &Path,
    consider_already_trimmed: Option<usize>,
) -> AdapterResult {
    eprintln!("\n\nAUTO-DETECTING ADAPTER TYPE\n===========================");
    eprintln!(
        "Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> {} <<)\n",
        file.display()
    );

    let mut candidates = vec![
        AdapterCandidate {
            name: "Illumina",
            display_name: "Illumina TruSeq, Sanger iPCR; auto-detected",
            sequence: "AGATCGGAAGAGC",
            count: 0,
        },
        AdapterCandidate {
            name: "Nextera",
            display_name: "Nextera Transposase sequence; auto-detected",
            sequence: "CTGTCTCTTATA",
            count: 0,
        },
        AdapterCandidate {
            name: "smallRNA",
            display_name: "Illumina small RNA adapter; auto-detected",
            sequence: "TGGAATTCTCGG",
            count: 0,
        },
    ];

    let reader: Box<dyn BufRead> =
        if file.extension().map_or(false, |e| e == "gz") {
            let f = File::open(file)
                .unwrap_or_else(|e| panic!("Failed to open '{}': {}", file.display(), e));
            Box::new(BufReader::new(GzDecoder::new(f)))
        } else {
            let f = File::open(file)
                .unwrap_or_else(|e| panic!("Failed to open '{}': {}", file.display(), e));
            Box::new(BufReader::new(f))
        };

    let mut total = 0usize;
    let mut lines = reader.lines();

    loop {
        // Read 4 lines (one FASTQ record)
        let _id = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        if lines.next().is_none() {
            break;
        }
        if lines.next().is_none() {
            break;
        }

        total += 1;
        if total >= 1_000_000 {
            break;
        }

        for c in &mut candidates {
            if seq.contains(c.sequence) {
                c.count += 1;
            }
        }
    }

    // Sort by count descending
    candidates.sort_by(|a, b| b.count.cmp(&a.count));

    eprintln!(
        "Found perfect matches for the following adapter sequences:\n\
         Adapter type\tCount\tSequence\tSequences analysed\tPercentage"
    );
    for c in &candidates {
        let pct = if total > 0 {
            c.count as f64 / total as f64 * 100.0
        } else {
            0.0
        };
        eprintln!(
            "{}\t{}\t{}\t{}\t{:.2}",
            c.name, c.count, c.sequence, total, pct
        );
    }

    let highest = &candidates[0];
    let second = &candidates[1];
    let mut report = String::new();

    // Check consider_already_trimmed threshold
    if let Some(threshold) = consider_already_trimmed {
        if highest.count <= threshold {
            let msg = format!(
                "No auto-detected adapter exceeded the 'already adapter-trimmed' limit ({}). Setting -a X\n",
                threshold
            );
            eprintln!("{}", msg);
            report.push_str(&msg);
            return AdapterResult {
                sequence: "X".to_string(),
                name: "No adapter trimming [suppressed by user]".to_string(),
                report_message: Some(report),
            };
        }
    }

    if highest.count == second.count {
        // Tie - default to Illumina
        let msg = format!(
            "Unable to auto-detect most prominent adapter (count {}: {}, count {}: {})\n\
             Defaulting to Illumina universal adapter (AGATCGGAAGAGC).\n\n",
            highest.name, highest.count, second.name, second.count
        );
        eprintln!("{}", msg);
        report.push_str(&msg);
        AdapterResult {
            sequence: "AGATCGGAAGAGC".to_string(),
            name: "Illumina TruSeq, Sanger iPCR; default (inconclusive auto-detection)"
                .to_string(),
            report_message: Some(report),
        }
    } else {
        let msg = format!(
            "Using {} adapter for trimming (count: {}). Second best hit was {} (count: {})\n\n",
            highest.name, highest.count, second.name, second.count
        );
        eprintln!("{}", msg);
        report.push_str(&msg);
        AdapterResult {
            sequence: highest.sequence.to_string(),
            name: highest.display_name.to_string(),
            report_message: Some(report),
        }
    }
}
