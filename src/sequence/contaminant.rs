use crate::config::Config;

pub struct ContaminantFinder {
    contaminants: Vec<(String, String)>, // (name, sequence)
}

impl ContaminantFinder {
    pub fn new(config: &Config) -> Self {
        let text = if let Some(ref path) = config.contaminant_file {
            std::fs::read_to_string(path).unwrap_or_default()
        } else {
            include_str!("../data/contaminant_list.txt").to_string()
        };

        let mut contaminants = Vec::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.splitn(2, '\t').collect();
            if parts.len() == 2 {
                contaminants.push((
                    parts[0].trim().to_string(),
                    parts[1].trim().to_string(),
                ));
            } else {
                // Try splitting on whitespace
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let seq = parts.last().unwrap().to_string();
                    let name = parts[..parts.len() - 1].join(" ");
                    contaminants.push((name, seq));
                }
            }
        }

        Self { contaminants }
    }

    /// Find a contaminant match for the given sequence.
    /// Returns the contaminant name and match details, or None.
    pub fn find_hit(&self, query: &str) -> Option<String> {
        for (name, contaminant_seq) in &self.contaminants {
            // Check forward match
            if let Some(hit) = Self::check_match(query, contaminant_seq) {
                return Some(format!("{} ({})", name, hit));
            }
            // Check reverse complement match
            let rc = reverse_complement(contaminant_seq);
            if let Some(hit) = Self::check_match(query, &rc) {
                return Some(format!("{} ({})", name, hit));
            }
        }
        None
    }

    fn check_match(query: &str, contaminant: &str) -> Option<String> {
        let query_upper = query.to_uppercase();
        let cont_upper = contaminant.to_uppercase();

        // Exact substring match (for short sequences)
        if query_upper.contains(&cont_upper) {
            return Some("100%".to_string());
        }
        if cont_upper.contains(&query_upper) {
            return Some("100%".to_string());
        }

        // For longer sequences, try partial matching with 1 mismatch tolerance
        if contaminant.len() >= 20 && query.len() >= 20 {
            let min_len = 20.min(query.len().min(contaminant.len()));
            let query_bytes = query_upper.as_bytes();
            let cont_bytes = cont_upper.as_bytes();

            // Try aligning the contaminant at each position in the query
            for offset in 0..=(query_bytes.len().saturating_sub(min_len)) {
                let check_len = min_len.min(query_bytes.len() - offset).min(cont_bytes.len());
                let mut mismatches = 0;
                for i in 0..check_len {
                    if query_bytes[offset + i] != cont_bytes[i] {
                        mismatches += 1;
                        if mismatches > 1 {
                            break;
                        }
                    }
                }
                if mismatches <= 1 && check_len >= 20 {
                    let pct = ((check_len - mismatches) * 100) / check_len;
                    return Some(format!("{}%", pct));
                }
            }
        }

        None
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            other => other,
        })
        .collect()
}
