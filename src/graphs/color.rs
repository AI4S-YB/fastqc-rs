/// Tol colorblind-safe palette for line graphs (matches Java's LineGraph colors).
pub const LINE_COLORS: &[(u8, u8, u8)] = &[
    (136, 34, 85),   // purple-red
    (51, 34, 136),   // dark blue
    (17, 119, 51),   // green
    (221, 204, 119), // yellow
    (68, 170, 153),  // teal
    (170, 68, 153),  // magenta
    (204, 102, 119), // pink
    (136, 204, 238), // light blue
];

/// Quality box plot background colors.
pub const QUALITY_GREEN: (u8, u8, u8) = (195, 230, 195);
pub const QUALITY_AMBER: (u8, u8, u8) = (230, 220, 195);
pub const QUALITY_RED: (u8, u8, u8) = (230, 195, 195);

/// Hot-Cold colour gradient for tile graphs.
/// Generates a 100-color palette from blue through green to red.
pub struct HotColdColourGradient {
    colors: Vec<(u8, u8, u8)>,
}

impl HotColdColourGradient {
    pub fn new() -> Self {
        let mut colors = Vec::with_capacity(100);
        for i in 0..100 {
            let fraction = i as f64 / 99.0;
            let r;
            let g;
            let b;

            if fraction < 0.5 {
                // Blue to Green
                let local = fraction * 2.0;
                r = 0;
                g = (255.0 * local) as u8;
                b = (255.0 * (1.0 - local)) as u8;
            } else {
                // Green to Red
                let local = (fraction - 0.5) * 2.0;
                r = (255.0 * local) as u8;
                g = (255.0 * (1.0 - local)) as u8;
                b = 0;
            }
            colors.push((r, g, b));
        }
        Self { colors }
    }

    /// Get color for a value between min and max, using square-root scaling.
    pub fn get_color(&self, value: f64, min: f64, max: f64) -> (u8, u8, u8) {
        if (max - min).abs() < f64::EPSILON {
            return self.colors[50];
        }

        let proportion = (value - min) / (max - min);
        let proportion = proportion.clamp(0.0, 1.0);

        // Square root scaling for better visual differentiation
        let scaled = if proportion >= 0.5 {
            0.5 + ((proportion - 0.5) * 2.0).sqrt() / 2.0
        } else {
            0.5 - ((0.5 - proportion) * 2.0).sqrt() / 2.0
        };

        let index = (scaled * 99.0) as usize;
        self.colors[index.min(99)]
    }
}

pub fn rgb_to_hex(r: u8, g: u8, b: u8) -> String {
    format!("#{:02x}{:02x}{:02x}", r, g, b)
}
