use super::color::HotColdColourGradient;
use super::color::rgb_to_hex;

/// Generate an SVG heatmap for per-tile quality scores.
pub fn render_tile_graph(
    x_labels: &[String],
    tiles: &[i32],
    means: &[Vec<f64>],
    width: usize,
    height: usize,
) -> String {
    let margin_left = 80;
    let margin_right = 20;
    let margin_top = 40;
    let margin_bottom = 100;

    let plot_width = width - margin_left - margin_right;
    let plot_height = height - margin_top - margin_bottom;

    let n_cols = x_labels.len();
    let n_rows = tiles.len();

    if n_cols == 0 || n_rows == 0 {
        return String::from(r##"<svg xmlns="http://www.w3.org/2000/svg" width="800" height="600"><text x="400" y="300" text-anchor="middle">No tile data</text></svg>"##);
    }

    // Find min/max deviation for color scaling
    let mut max_dev = 0.0f64;
    for row in means {
        for &val in row {
            if val.abs() > max_dev {
                max_dev = val.abs();
            }
        }
    }
    if max_dev < f64::EPSILON {
        max_dev = 1.0;
    }

    let gradient = HotColdColourGradient::new();

    let cell_width = plot_width as f64 / n_cols as f64;
    let cell_height = plot_height as f64 / n_rows as f64;

    let mut svg = String::new();
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">"##,
        width, height, width, height
    ));

    svg.push_str(&format!(
        r##"<rect x="0" y="0" width="{}" height="{}" fill="white"/>"##,
        width, height
    ));

    svg.push_str(&format!(
        r##"<text x="{}" y="20" text-anchor="middle" font-size="14" font-weight="bold">Quality per tile</text>"##,
        width / 2
    ));

    // Heatmap cells
    for (t, row) in means.iter().enumerate() {
        for (p, &val) in row.iter().enumerate() {
            let x = margin_left as f64 + p as f64 * cell_width;
            let y = margin_top as f64 + t as f64 * cell_height;
            let (r, g, b) = gradient.get_color(val, -max_dev, max_dev);
            svg.push_str(&format!(
                r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}"/>"##,
                x,
                y,
                cell_width.ceil(),
                cell_height.ceil(),
                rgb_to_hex(r, g, b)
            ));
        }
    }

    // Y-axis tile labels (show a subset if too many)
    let y_step = if n_rows > 50 { n_rows / 25 } else { 1 };
    for i in (0..n_rows).step_by(y_step.max(1)) {
        let y = margin_top as f64 + (i as f64 + 0.5) * cell_height;
        svg.push_str(&format!(
            r##"<text x="{}" y="{:.0}" text-anchor="end" font-size="9" fill="#333">{}</text>"##,
            margin_left - 5,
            y + 3.0,
            tiles[i]
        ));
    }

    // X-axis labels
    let x_step = if n_cols > 50 { n_cols / 25 } else { 1 };
    for i in (0..n_cols).step_by(x_step.max(1)) {
        let x = margin_left as f64 + (i as f64 + 0.5) * cell_width;
        let y = (margin_top + plot_height + 15) as f64;
        svg.push_str(&format!(
            r##"<text x="{:.0}" y="{:.0}" text-anchor="end" font-size="9" fill="#333" transform="rotate(-90, {:.0}, {:.0})">{}</text>"##,
            x, y, x, y,
            escape_xml(&x_labels[i])
        ));
    }

    svg.push_str("</svg>");
    svg
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
