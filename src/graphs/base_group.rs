use crate::config::Config;

/// Represents a group of base positions for visualization.
/// Matches Java's BaseGroup class.
#[derive(Debug, Clone)]
pub struct BaseGroup {
    lower: usize,
    upper: usize,
}

impl BaseGroup {
    pub fn new(lower: usize, upper: usize) -> Self {
        Self { lower, upper }
    }

    /// 1-based lower bound (inclusive).
    pub fn lower_count(&self) -> usize {
        self.lower
    }

    /// 1-based upper bound (inclusive).
    pub fn upper_count(&self) -> usize {
        self.upper
    }

    /// Check if a 1-based position falls within this group.
    pub fn contains_value(&self, value: usize) -> bool {
        value >= self.lower && value <= self.upper
    }
}

impl std::fmt::Display for BaseGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.lower == self.upper {
            write!(f, "{}", self.lower)
        } else {
            write!(f, "{}-{}", self.lower, self.upper)
        }
    }
}

/// Create base groups using the configured strategy.
pub fn make_base_groups(max_length: usize, config: &Config) -> Vec<BaseGroup> {
    let effective_max = if config.min_length > max_length {
        config.min_length
    } else {
        max_length
    };

    if config.nogroup {
        make_ungrouped(effective_max)
    } else if config.expgroup {
        make_exponential_base_groups(effective_max)
    } else {
        make_linear_base_groups(effective_max)
    }
}

/// Every position gets its own group.
pub fn make_ungrouped(max_length: usize) -> Vec<BaseGroup> {
    (1..=max_length).map(|i| BaseGroup::new(i, i)).collect()
}

/// Linear grouping: positions 1-9 individual, then fixed-width intervals.
/// Targets fewer than 75 total groups.
/// Matches Java's BaseGroup.makeLinearBaseGroups().
pub fn make_linear_base_groups(max_length: usize) -> Vec<BaseGroup> {
    if max_length <= 75 {
        return make_ungrouped(max_length);
    }

    // Find the right interval
    let base_values = [2usize, 5, 10];
    let mut interval = 1;

    'outer: for &multiplier in &[1, 10, 100, 1000, 10000, 100000] {
        for &base in &base_values {
            let candidate = base * multiplier;
            // Calculate how many groups this would produce
            let groups_after_9 = if max_length > 9 {
                (max_length - 9 + candidate - 1) / candidate
            } else {
                0
            };
            let total_groups = 9.min(max_length) + groups_after_9;
            if total_groups <= 75 {
                interval = candidate;
                break 'outer;
            }
        }
    }

    let mut groups = Vec::new();

    // First 9 positions are individual
    let individual_end = 9.min(max_length);
    for i in 1..=individual_end {
        groups.push(BaseGroup::new(i, i));
    }

    if max_length <= 9 {
        return groups;
    }

    // Remaining positions grouped by interval
    let mut start = 10;
    while start <= max_length {
        let end = (start + interval - 1).min(max_length);
        groups.push(BaseGroup::new(start, end));
        start = end + 1;
    }

    groups
}

/// Exponential grouping: interval doubles at certain thresholds.
/// Matches Java's BaseGroup.makeExponentialBaseGroups().
pub fn make_exponential_base_groups(max_length: usize) -> Vec<BaseGroup> {
    let mut groups = Vec::new();
    let mut current_position = 1;
    let mut current_interval = 1;

    while current_position <= max_length {
        // Check if we need to increase the interval
        if current_position == 10 && max_length > 75 {
            current_interval = 2;
        }
        if current_position == 50 && max_length > 200 {
            current_interval = 5;
        }
        if current_position == 100 && max_length > 300 {
            current_interval = 10;
        }
        if current_position == 500 && max_length > 1000 {
            current_interval = 50;
        }
        if current_position == 1000 && max_length > 2000 {
            current_interval = 100;
        }

        let end = (current_position + current_interval - 1).min(max_length);
        groups.push(BaseGroup::new(current_position, end));
        current_position = end + 1;
    }

    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ungrouped() {
        let groups = make_ungrouped(5);
        assert_eq!(groups.len(), 5);
        assert_eq!(groups[0].to_string(), "1");
        assert_eq!(groups[4].to_string(), "5");
    }

    #[test]
    fn test_linear_short_reads() {
        // Reads <= 75bp should not be grouped
        let groups = make_linear_base_groups(50);
        assert_eq!(groups.len(), 50);
    }

    #[test]
    fn test_linear_long_reads() {
        let groups = make_linear_base_groups(150);
        // Should have 9 individual + some grouped
        assert!(groups.len() <= 75);
        assert_eq!(groups[0].lower_count(), 1);
        assert_eq!(groups[0].upper_count(), 1);
        assert_eq!(groups[8].lower_count(), 9);
        assert_eq!(groups[8].upper_count(), 9);
        // Last group should end at max_length
        assert_eq!(groups.last().unwrap().upper_count(), 150);
    }

    #[test]
    fn test_exponential() {
        let groups = make_exponential_base_groups(200);
        // Exponential grouping reduces count but may still be > 75 for moderate lengths
        assert!(groups.len() < 200);
        assert_eq!(groups[0].lower_count(), 1);
        assert_eq!(groups.last().unwrap().upper_count(), 200);
    }

    #[test]
    fn test_contains_value() {
        let group = BaseGroup::new(10, 20);
        assert!(group.contains_value(10));
        assert!(group.contains_value(15));
        assert!(group.contains_value(20));
        assert!(!group.contains_value(9));
        assert!(!group.contains_value(21));
    }
}
