use std::f64::consts::PI;

/// Normal distribution, matching Java's NormalDistribution used in PerSequenceGCContent.
/// The getZScoreForValue method in Java actually computes the PDF value, not a z-score.
pub struct NormalDistribution {
    mean: f64,
    stdev: f64,
}

impl NormalDistribution {
    pub fn new(mean: f64, stdev: f64) -> Self {
        Self { mean, stdev }
    }

    /// Compute the probability density function value at x.
    /// (Java calls this getZScoreForValue, but it's actually the PDF)
    pub fn pdf(&self, x: f64) -> f64 {
        if self.stdev <= 0.0 {
            return 0.0;
        }
        let exponent = -((x - self.mean).powi(2)) / (2.0 * self.stdev.powi(2));
        (1.0 / (self.stdev * (2.0 * PI).sqrt())) * exponent.exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_pdf() {
        let nd = NormalDistribution::new(50.0, 10.0);
        let peak = nd.pdf(50.0);
        let off_center = nd.pdf(60.0);
        assert!(peak > off_center);
        assert!(peak > 0.0);
    }
}
