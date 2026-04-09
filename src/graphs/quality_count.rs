/// Tracks quality score character frequencies at a single base position.
/// Uses a fixed 150-slot array (ASCII values 0-149), matching Java's QualityCount.
#[derive(Debug, Clone)]
pub struct QualityCount {
    actual_counts: [u64; 150],
    total_counts: u64,
}

impl QualityCount {
    pub fn new() -> Self {
        Self {
            actual_counts: [0u64; 150],
            total_counts: 0,
        }
    }

    /// Add a quality character value.
    pub fn add_value(&mut self, c: u8) {
        self.total_counts += 1;
        let idx = c as usize;
        if idx >= self.actual_counts.len() {
            panic!(
                "Got character {} as a quality value which has ASCII {} which is higher than we can cope with",
                c as char, idx
            );
        }
        self.actual_counts[idx] += 1;
    }

    pub fn total_count(&self) -> u64 {
        self.total_counts
    }

    /// Find the lowest quality character with non-zero count.
    pub fn min_char(&self) -> u8 {
        for i in 0..self.actual_counts.len() {
            if self.actual_counts[i] > 0 {
                return i as u8;
            }
        }
        // Return a high value if nothing found (matches Java's (char)1000 -> wraps)
        255
    }

    /// Find the highest quality character with non-zero count.
    pub fn max_char(&self) -> u8 {
        for i in (0..self.actual_counts.len()).rev() {
            if self.actual_counts[i] > 0 {
                return i as u8;
            }
        }
        255
    }

    /// Calculate mean quality score, adjusted by the given offset.
    /// Matches Java: sum(count[i] * (i - offset)) / total_count
    pub fn mean(&self, offset: usize) -> f64 {
        let mut total: i64 = 0;
        let mut count: u64 = 0;

        for i in offset..self.actual_counts.len() {
            total += self.actual_counts[i] as i64 * (i - offset) as i64;
            count += self.actual_counts[i];
        }

        if count > 0 {
            total as f64 / count as f64
        } else {
            0.0
        }
    }

    /// Calculate a percentile quality score, adjusted by the given offset.
    /// Matches Java's integer-division approach for the threshold.
    pub fn percentile(&self, offset: usize, percentile: u32) -> f64 {
        // Java: total *= percentile; total /= 100; (integer operations on long)
        let mut target = self.total_counts;
        target *= percentile as u64;
        target /= 100;

        let mut count: u64 = 0;
        for i in offset..self.actual_counts.len() {
            count += self.actual_counts[i];
            if count >= target {
                return (i - offset) as f64;
            }
        }

        -1.0
    }
}

impl Default for QualityCount {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_value_and_counts() {
        let mut qc = QualityCount::new();
        qc.add_value(b'I'); // ASCII 73
        qc.add_value(b'I');
        qc.add_value(b'5'); // ASCII 53
        assert_eq!(qc.total_count(), 3);
        assert_eq!(qc.min_char(), b'5');
        assert_eq!(qc.max_char(), b'I');
    }

    #[test]
    fn test_mean() {
        let mut qc = QualityCount::new();
        // Add 10 values all at quality 40 with Sanger offset (33)
        for _ in 0..10 {
            qc.add_value(33 + 40); // ASCII 73
        }
        let mean = qc.mean(33);
        assert!((mean - 40.0).abs() < 0.001);
    }

    #[test]
    fn test_percentile() {
        let mut qc = QualityCount::new();
        // Add values from 33 to 72 (quality 0-39), one each
        for i in 33..=72 {
            qc.add_value(i);
        }
        let median = qc.percentile(33, 50);
        // 40 values, 50th percentile -> target = 40*50/100 = 20
        // Counting from 33: position 20 means quality 19
        assert!((median - 19.0).abs() < 1.0);
    }
}
