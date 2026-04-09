/// GC content model for mapping GC counts to percentage bins.
/// Matches Java's GCModel class.
#[derive(Clone)]
pub struct GcModel {
    /// models[gc_count] = Vec<(percentage, increment)>
    models: Vec<Vec<(usize, f64)>>,
}

impl GcModel {
    pub fn new(read_length: usize) -> Self {
        // First pass: count how many GC counts claim each percentage bucket
        let mut claiming_counts = vec![0usize; 101];

        for pos in 0..=read_length {
            let low_count = if pos == 0 { 0.0 } else { pos as f64 - 0.5 };
            let high_count = (pos as f64 + 0.5).min(read_length as f64);

            let low_pct = ((low_count * 100.0) / read_length as f64).round() as usize;
            let high_pct = ((high_count * 100.0) / read_length as f64).round() as usize;

            for p in low_pct..=high_pct.min(100) {
                claiming_counts[p] += 1;
            }
        }

        // Second pass: build the model with weighted increments
        let mut models = Vec::with_capacity(read_length + 1);

        for pos in 0..=read_length {
            let low_count = if pos == 0 { 0.0 } else { pos as f64 - 0.5 };
            let high_count = (pos as f64 + 0.5).min(read_length as f64);

            let low_pct = ((low_count * 100.0) / read_length as f64).round() as usize;
            let high_pct = ((high_count * 100.0) / read_length as f64).round() as usize;

            let mut values = Vec::new();
            for p in low_pct..=high_pct.min(100) {
                values.push((p, 1.0 / claiming_counts[p] as f64));
            }
            models.push(values);
        }

        Self { models }
    }

    pub fn get_model_values(&self, gc_count: usize) -> &[(usize, f64)] {
        if gc_count < self.models.len() {
            &self.models[gc_count]
        } else {
            &[]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_model_basic() {
        let model = GcModel::new(10);
        // GC count 0 should map to low percentages
        let values = model.get_model_values(0);
        assert!(!values.is_empty());
        assert!(values[0].0 <= 5); // Should be near 0%

        // GC count 10 should map to high percentages
        let values = model.get_model_values(10);
        assert!(!values.is_empty());
        assert!(values.last().unwrap().0 >= 95); // Should be near 100%
    }

    #[test]
    fn test_gc_model_sum() {
        let model = GcModel::new(100);
        // For each percentage bucket, the sum of increments across all GC counts should be ~1
        let mut bucket_sums = vec![0.0f64; 101];
        for gc in 0..=100 {
            for &(pct, inc) in model.get_model_values(gc) {
                bucket_sums[pct] += inc;
            }
        }
        for (i, &sum) in bucket_sums.iter().enumerate() {
            assert!(
                (sum - 1.0).abs() < 0.01,
                "Bucket {} sum = {}, expected ~1.0",
                i,
                sum
            );
        }
    }
}
