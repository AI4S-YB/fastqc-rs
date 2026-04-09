use super::color::{rgb_to_hex, QUALITY_AMBER, QUALITY_GREEN, QUALITY_RED};

/// Generate an SVG quality box plot.
pub fn render_quality_box_plot(
    means: &[f64],
    medians: &[f64],
    lowest: &[f64],
    highest: &[f64],
    lower_quartile: &[f64],
    upper_quartile: &[f64],
    y_min: f64,
    y_max: f64,
    x_labels: &[String],
    title: &str,
    width: usize,
    height: usize,
) -> String {
    let margin_left = 60;
    let margin_right = 20;
    let margin_top = 40;
    let margin_bottom = if x_labels.len() > 20 { 100 } else { 60 };

    let plot_width = width - margin_left - margin_right;
    let plot_height = height - margin_top - margin_bottom;

    let y_range = if (y_max - y_min).abs() < f64::EPSILON {
        1.0
    } else {
        y_max - y_min
    };

    let n_boxes = means.len();
    let box_width = if n_boxes > 0 {
        (plot_width as f64 / n_boxes as f64 * 0.8).max(3.0)
    } else {
        10.0
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

    // Quality zone backgrounds
    let quality_zones: [(f64, f64, (u8, u8, u8)); 3] = [
        (0.0, 20.0, QUALITY_RED),
        (20.0, 28.0, QUALITY_AMBER),
        (28.0, y_max, QUALITY_GREEN),
    ];

    for &(zone_min, zone_max, (r, g, b)) in quality_zones.iter() {
        let zone_min = zone_min.max(y_min);
        let zone_max = zone_max.min(y_max);
        if zone_max <= zone_min {
            continue;
        }
        let y_top = margin_top as f64 + plot_height as f64 * (1.0 - (zone_max - y_min) / y_range);
        let y_bottom = margin_top as f64 + plot_height as f64 * (1.0 - (zone_min - y_min) / y_range);
        svg.push_str(&format!(
            r##"<rect x="{}" y="{:.0}" width="{}" height="{:.0}" fill="{}"/>"##,
            margin_left,
            y_top,
            plot_width,
            y_bottom - y_top,
            rgb_to_hex(r, g, b)
        ));
    }

    // Plot border
    svg.push_str(&format!(
        r##"<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="#cccccc"/>"##,
        margin_left, margin_top, plot_width, plot_height
    ));

    // Y-axis gridlines
    let y_ticks = 5;
    for i in 0..=y_ticks {
        let y_val = y_min + y_range * i as f64 / y_ticks as f64;
        let y_pos = margin_top as f64 + plot_height as f64 * (1.0 - i as f64 / y_ticks as f64);
        svg.push_str(&format!(
            r##"<text x="{}" y="{:.0}" text-anchor="end" font-size="10" fill="#333">{:.0}</text>"##,
            margin_left - 5,
            y_pos + 4.0,
            y_val
        ));
    }

    // Box plots
    for i in 0..n_boxes {
        let x_center =
            margin_left as f64 + plot_width as f64 * (i as f64 + 0.5) / n_boxes as f64;

        let val_to_y = |val: f64| -> f64 {
            let v = val.clamp(y_min, y_max);
            margin_top as f64 + plot_height as f64 * (1.0 - (v - y_min) / y_range)
        };

        if medians[i].is_nan() {
            continue;
        }

        let y_low = val_to_y(lowest[i]);
        let y_high = val_to_y(highest[i]);
        let y_lq = val_to_y(lower_quartile[i]);
        let y_uq = val_to_y(upper_quartile[i]);
        let y_med = val_to_y(medians[i]);
        let y_mean = val_to_y(means[i]);

        let half_box = box_width / 2.0;

        // Whisker lines
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#333" stroke-width="1"/>"##,
            x_center, y_low, x_center, y_lq
        ));
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#333" stroke-width="1"/>"##,
            x_center, y_high, x_center, y_uq
        ));

        // Whisker caps
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#333" stroke-width="1"/>"##,
            x_center - half_box * 0.5, y_low, x_center + half_box * 0.5, y_low
        ));
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#333" stroke-width="1"/>"##,
            x_center - half_box * 0.5, y_high, x_center + half_box * 0.5, y_high
        ));

        // Box (IQR)
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="#ffdf00" stroke="#333" stroke-width="1"/>"##,
            x_center - half_box,
            y_uq,
            box_width,
            (y_lq - y_uq).abs().max(1.0)
        ));

        // Median line (red)
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="red" stroke-width="2"/>"##,
            x_center - half_box, y_med, x_center + half_box, y_med
        ));

        // Mean point
        svg.push_str(&format!(
            r##"<circle cx="{:.1}" cy="{:.1}" r="2" fill="blue"/>"##,
            x_center, y_mean
        ));
    }

    // X-axis labels
    let step = if n_boxes > 50 { n_boxes / 25 } else { 1 };
    for i in (0..n_boxes).step_by(step.max(1)) {
        let x_pos = margin_left as f64 + plot_width as f64 * (i as f64 + 0.5) / n_boxes as f64;
        let y_pos = (margin_top + plot_height + 15) as f64;

        if n_boxes > 20 {
            svg.push_str(&format!(
                r##"<text x="{:.0}" y="{:.0}" text-anchor="end" font-size="9" fill="#333" transform="rotate(-90, {:.0}, {:.0})">{}</text>"##,
                x_pos, y_pos, x_pos, y_pos, escape_xml(&x_labels[i])
            ));
        } else {
            svg.push_str(&format!(
                r##"<text x="{:.0}" y="{:.0}" text-anchor="middle" font-size="10" fill="#333">{}</text>"##,
                x_pos, y_pos, escape_xml(&x_labels[i])
            ));
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
