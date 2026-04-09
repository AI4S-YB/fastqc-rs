use std::path::{Path, PathBuf};

/// Group files by CASAVA sample name.
pub fn get_casava_groups(files: &[PathBuf]) -> Vec<Vec<PathBuf>> {
    // Simplified: just return each file individually for now
    // Full CASAVA grouping requires parsing sample names from filenames
    files.iter().map(|f| vec![f.clone()]).collect()
}

/// Check if a filename matches CASAVA naming convention.
pub fn is_casava_file(path: &Path) -> bool {
    if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
        // CASAVA files match: SampleName_S1_L001_R1_001.fastq.gz
        name.contains("_S") && name.contains("_L") && name.contains("_R")
    } else {
        false
    }
}
