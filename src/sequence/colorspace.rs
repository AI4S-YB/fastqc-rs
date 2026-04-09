/// Detect if a sequence is in colorspace format (SOLiD platform).
/// Pattern: first char is DNA base, remaining chars are color codes (0-6 or .)
pub fn is_colorspace(sequence: &str) -> bool {
    if sequence.len() < 2 {
        return false;
    }
    let chars: Vec<char> = sequence.chars().collect();
    // First char must be a DNA base
    if !"GATCNgatcn".contains(chars[0]) {
        return false;
    }
    // Remaining chars must be color codes
    for &c in &chars[1..] {
        if !"0123456.".contains(c) {
            return false;
        }
    }
    true
}

/// Convert a colorspace sequence to DNA bases.
/// Returns the converted sequence (length is original - 1).
pub fn colorspace_to_bases(colorspace: &str) -> String {
    let chars: Vec<char> = colorspace.chars().collect();
    if chars.len() < 2 {
        return String::new();
    }

    let mut result = Vec::with_capacity(chars.len() - 1);
    let mut last_base = chars[0].to_ascii_uppercase();

    for &color in &chars[1..] {
        let next_base = convert_base(last_base, color);
        result.push(next_base);
        last_base = next_base;
    }

    result.iter().collect()
}

fn convert_base(reference: char, color: char) -> char {
    match color {
        '0' => reference, // Same base
        '1' => match reference {
            'A' => 'C',
            'C' => 'A',
            'G' => 'T',
            'T' => 'G',
            _ => 'N',
        },
        '2' => match reference {
            'A' => 'G',
            'G' => 'A',
            'C' => 'T',
            'T' => 'C',
            _ => 'N',
        },
        '3' => match reference {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        },
        _ => 'N', // '.', '4', '5', '6' all produce N
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_colorspace() {
        assert!(is_colorspace("G0123012"));
        assert!(is_colorspace("A..0123"));
        assert!(!is_colorspace("ATCGATCG"));
        assert!(!is_colorspace(""));
        assert!(!is_colorspace("X0123"));
    }

    #[test]
    fn test_conversion() {
        // G followed by 0=G, 1=T, 2=A, 3=C
        // G->0->G, G->1->T, T->2->C, C->3->G
        let result = colorspace_to_bases("G0123");
        assert_eq!(result, "GTCG");
    }
}
