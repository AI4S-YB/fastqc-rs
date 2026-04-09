pub mod html;
pub mod text_report;
pub mod summary;

use crate::config::Config;

/// Report archive that collects module output data.
/// Modules write their data/html during make_report().
pub struct ReportArchive {
    /// Text data for fastqc_data.txt
    pub data: String,
    /// HTML body content for modules
    pub html_body: String,
    /// SVG images: (filename, svg_content)
    pub images: Vec<(String, String)>,
    /// Reference to config for grouping calculations
    pub config: Config,
    /// Sequence file name
    pub file_name: String,
}

impl ReportArchive {
    pub fn new(file_name: &str, config: Config) -> Self {
        Self {
            data: String::new(),
            html_body: String::new(),
            images: Vec::new(),
            config,
            file_name: file_name.to_string(),
        }
    }

    pub fn add_image(&mut self, filename: &str, svg_content: &str) {
        self.images.push((filename.to_string(), svg_content.to_string()));
    }
}
