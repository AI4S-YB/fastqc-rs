use std::fmt;

/// Phred encoding schemes, matching Java's PhredEncoding class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhredEncoding {
    Sanger,
    Illumina1_3,
    Illumina1_5,
}

impl PhredEncoding {
    /// Determine the encoding scheme from the lowest quality character seen.
    pub fn from_lowest_char(lowest_char: u8) -> Self {
        if lowest_char < 33 {
            // Something's wrong, default to Sanger
            PhredEncoding::Sanger
        } else if lowest_char < 64 {
            PhredEncoding::Sanger
        } else if lowest_char >= 66 {
            PhredEncoding::Illumina1_5
        } else {
            PhredEncoding::Illumina1_3
        }
    }

    /// Get the ASCII offset for this encoding.
    pub fn offset(&self) -> usize {
        match self {
            PhredEncoding::Sanger => 33,
            PhredEncoding::Illumina1_3 => 64,
            PhredEncoding::Illumina1_5 => 64,
        }
    }
}

impl fmt::Display for PhredEncoding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PhredEncoding::Sanger => write!(f, "Sanger / Illumina 1.9"),
            PhredEncoding::Illumina1_3 => write!(f, "Illumina 1.3"),
            PhredEncoding::Illumina1_5 => write!(f, "Illumina 1.5"),
        }
    }
}
