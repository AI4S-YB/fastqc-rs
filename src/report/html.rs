use base64::Engine;
use base64::engine::general_purpose::STANDARD as BASE64;

use crate::modules::QcModule;

const HEADER_CSS: &str = include_str!("../data/header_template.html");
const ICON_FASTQC: &[u8] = include_bytes!("../data/icons/fastqc_icon.png");
const ICON_TICK: &[u8] = include_bytes!("../data/icons/tick.png");
const ICON_WARNING: &[u8] = include_bytes!("../data/icons/warning.png");
const ICON_ERROR: &[u8] = include_bytes!("../data/icons/error.png");

fn png_to_base64(data: &[u8]) -> String {
    format!("data:image/png;base64,{}", BASE64.encode(data))
}

fn svg_to_base64(svg: &str) -> String {
    format!("data:image/svg+xml;base64,{}", BASE64.encode(svg.as_bytes()))
}

/// Generate a complete HTML report document.
pub fn generate_html(
    file_name: &str,
    modules: &mut [Box<dyn QcModule>],
    module_data: &[(String, String, Vec<(String, String)>)], // (data, html_body, images) per module
) -> String {
    let mut html = String::new();

    html.push_str("<!DOCTYPE html><html><head>");
    html.push_str(&format!("<title>{} FastQC Report</title>", escape_xml(file_name)));
    html.push_str("<style type=\"text/css\">");
    html.push_str(HEADER_CSS);
    html.push_str("</style></head><body>");

    // Header
    html.push_str("<div class=\"header\">");
    html.push_str("<div id=\"header_title\">");
    html.push_str(&format!(
        "<img src=\"{}\" alt=\"FastQC\"/>FastQC Report",
        png_to_base64(ICON_FASTQC)
    ));
    html.push_str("</div>");
    html.push_str("<div id=\"header_filename\">");
    html.push_str(&chrono_date());
    html.push_str("<br/>");
    html.push_str(&escape_xml(file_name));
    html.push_str("</div></div>");

    // Summary sidebar
    html.push_str("<div class=\"summary\"><h2>Summary</h2><ul>");
    for (m, module) in modules.iter_mut().enumerate() {
        if module.ignore_in_report() {
            continue;
        }
        let (icon, alt) = if module.raises_error() {
            (png_to_base64(ICON_ERROR), "[FAIL]")
        } else if module.raises_warning() {
            (png_to_base64(ICON_WARNING), "[WARNING]")
        } else {
            (png_to_base64(ICON_TICK), "[PASS]")
        };
        html.push_str(&format!(
            "<li><img src=\"{}\" alt=\"{}\"/><a href=\"#M{}\">{}</a></li>",
            icon,
            alt,
            m,
            escape_xml(module.name())
        ));
    }
    html.push_str("</ul></div>");

    // Main content
    html.push_str("<div class=\"main\">");
    for (m, module) in modules.iter_mut().enumerate() {
        if module.ignore_in_report() {
            continue;
        }

        let (icon, alt) = if module.raises_error() {
            (png_to_base64(ICON_ERROR), "[FAIL]")
        } else if module.raises_warning() {
            (png_to_base64(ICON_WARNING), "[WARN]")
        } else {
            (png_to_base64(ICON_TICK), "[OK]")
        };

        html.push_str("<div class=\"module\">");
        html.push_str(&format!(
            "<h2 id=\"M{}\"><img src=\"{}\" alt=\"{}\"/>{}</h2>",
            m,
            icon,
            alt,
            escape_xml(module.name())
        ));

        // Module content (images + html body)
        if m < module_data.len() {
            let (_, ref body, ref images) = module_data[m];
            // Embed images
            for (_, svg_content) in images {
                html.push_str(&format!(
                    "<p><img src=\"{}\" alt=\"Graph\"/></p>",
                    svg_to_base64(svg_content)
                ));
            }
            html.push_str(body);
        }

        html.push_str("</div>");
    }
    html.push_str("</div>");

    // Footer
    html.push_str("<div class=\"footer\">");
    html.push_str(&format!(
        "Produced by <a href=\"http://www.bioinformatics.babraham.ac.uk/projects/fastqc/\">FastQC</a>  (version {})",
        crate::VERSION
    ));
    html.push_str("</div>");

    html.push_str("</body></html>");
    html
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn chrono_date() -> String {
    // Simple date without external dependency
    let now = std::time::SystemTime::now();
    let since_epoch = now
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = since_epoch.as_secs();
    // Simple formatting: just show the timestamp
    // For proper formatting, would need chrono crate
    let days = secs / 86400;
    let year = 1970 + (days / 365); // Approximate
    format!("{}", year)
}
