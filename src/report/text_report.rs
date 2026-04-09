use crate::modules::QcModule;

/// Generate fastqc_data.txt content.
pub fn generate_text_report(
    modules: &mut [Box<dyn QcModule>],
    module_data: &[String], // data string per module
) -> String {
    let mut text = String::new();

    // Header
    text.push_str(&format!("##FastQC\t{}\n", crate::VERSION));

    // Module sections
    for (m, module) in modules.iter_mut().enumerate() {
        if module.ignore_in_report() {
            continue;
        }

        text.push_str(">>");
        text.push_str(module.name());
        text.push('\t');
        if module.raises_error() {
            text.push_str("fail");
        } else if module.raises_warning() {
            text.push_str("warn");
        } else {
            text.push_str("pass");
        }
        text.push('\n');

        if m < module_data.len() {
            text.push_str(&module_data[m]);
        }

        text.push_str(">>END_MODULE\n");
    }

    text
}
