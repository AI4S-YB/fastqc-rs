use crate::modules::QcModule;

/// Generate summary.txt content.
pub fn generate_summary(
    modules: &mut [Box<dyn QcModule>],
    file_name: &str,
) -> String {
    let mut text = String::new();

    for module in modules.iter_mut() {
        if module.ignore_in_report() {
            continue;
        }

        if module.raises_error() {
            text.push_str("FAIL");
        } else if module.raises_warning() {
            text.push_str("WARN");
        } else {
            text.push_str("PASS");
        }
        text.push('\t');
        text.push_str(module.name());
        text.push('\t');
        text.push_str(file_name);
        text.push('\n');
    }

    text
}
