use super::color::{rgb_to_hex, LINE_COLORS};

/// Generate an SVG line graph with multiple data series.
pub fn render_line_graph(
    data: &[&[f64]],
    y_min: f64,
    y_max: f64,
    x_label: &str,
    series_labels: &[&str],
    x_categories: &[String],
    title: &str,
    width: usize,
    height: usize,
) -> String {
    let margin_left = 80;
    let margin_right = 20;
    let margin_top = 40;
    let margin_bottom = if x_categories.len() > 20 { 100 } else { 60 };
    let legend_height = if series_labels.len() > 1 { 30 } else { 0 };

    let plot_width = width - margin_left - margin_right;
    let plot_height = height - margin_top - margin_bottom - legend_height;

    let y_range = if (y_max - y_min).abs() < f64::EPSILON {
        1.0
    } else {
        y_max - y_min
    };

    let mut svg = String::new();
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">"##,
        width, height, width, height
    ));

    // Background
    svg.push_str(&format!(
        r##"<rect x="0" y="0" width="{}" height="{}" fill="white"/>"##,
        width, height
    ));

    // Title
    svg.push_str(&format!(
        r##"<text x="{}" y="20" text-anchor="middle" font-size="14" font-weight="bold">{}</text>"##,
        width / 2,
        escape_xml(title)
    ));

    // Plot area background
    svg.push_str(&format!(
        r##"<rect x="{}" y="{}" width="{}" height="{}" fill="#f8f8f8" stroke="#cccccc"/>"##,
        margin_left, margin_top, plot_width, plot_height
    ));

    // Y-axis gridlines and labels
    let y_ticks = 5;
    for i in 0..=y_ticks {
        let y_val = y_min + y_range * i as f64 / y_ticks as f64;
        let y_pos = margin_top + plot_height - (plot_height as f64 * i as f64 / y_ticks as f64) as usize;
        svg.push_str(&format!(
            r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#e0e0e0" stroke-width="1"/>"##,
            margin_left, y_pos, margin_left + plot_width, y_pos
        ));
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" text-anchor="end" font-size="10" fill="#333">{:.1}</text>"##,
            margin_left - 5,
            y_pos + 4,
            y_val
        ));
    }

    // X-axis labels
    let n_points = x_categories.len();
    if n_points > 0 {
        let step = if n_points > 50 { n_points / 25 } else { 1 };
        for i in (0..n_points).step_by(step.max(1)) {
            let x_pos = margin_left + (plot_width as f64 * (i as f64 + 0.5) / n_points as f64) as usize;
            let y_pos = margin_top + plot_height + 15;

            if x_categories.len() > 20 {
                svg.push_str(&format!(
                    r##"<text x="{}" y="{}" text-anchor="end" font-size="9" fill="#333" transform="rotate(-90, {}, {})">{}</text>"##,
                    x_pos, y_pos, x_pos, y_pos, escape_xml(&x_categories[i])
                ));
            } else {
                svg.push_str(&format!(
                    r##"<text x="{}" y="{}" text-anchor="middle" font-size="10" fill="#333">{}</text>"##,
                    x_pos, y_pos, escape_xml(&x_categories[i])
                ));
            }
        }
    }

    // X-axis label
    svg.push_str(&format!(
        r##"<text x="{}" y="{}" text-anchor="middle" font-size="12" fill="#333">{}</text>"##,
        margin_left + plot_width / 2,
        height - 5,
        escape_xml(x_label)
    ));

    // Data lines
    for (series_idx, series) in data.iter().enumerate() {
        let color_idx = series_idx % LINE_COLORS.len();
        let (r, g, b) = LINE_COLORS[color_idx];
        let color = rgb_to_hex(r, g, b);

        let mut path = String::new();
        for (i, &val) in series.iter().enumerate() {
            let x = margin_left as f64 + plot_width as f64 * (i as f64 + 0.5) / n_points as f64;
            let y = margin_top as f64 + plot_height as f64 * (1.0 - (val - y_min) / y_range);
            let y = y.clamp(margin_top as f64, (margin_top + plot_height) as f64);

            if i == 0 {
                path.push_str(&format!("M{:.1},{:.1}", x, y));
            } else {
                path.push_str(&format!(" L{:.1},{:.1}", x, y));
            }
        }
        svg.push_str(&format!(
            r##"<path d="{}" fill="none" stroke="{}" stroke-width="2"/>"##,
            path, color
        ));
    }

    // Legend
    if series_labels.len() > 1 {
        let legend_y = margin_top + plot_height + margin_bottom - 10;
        let legend_x_start = margin_left;
        let mut x_offset = legend_x_start;

        for (i, label) in series_labels.iter().enumerate() {
            let color_idx = i % LINE_COLORS.len();
            let (r, g, b) = LINE_COLORS[color_idx];
            let color = rgb_to_hex(r, g, b);

            svg.push_str(&format!(
                r##"<rect x="{}" y="{}" width="12" height="12" fill="{}"/>"##,
                x_offset,
                legend_y - 10,
                color
            ));
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" font-size="10" fill="#333">{}</text>"##,
                x_offset + 16,
                legend_y,
                escape_xml(label)
            ));
            x_offset += 16 + label.len() * 7 + 20;
        }
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
