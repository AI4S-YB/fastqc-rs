use std::path::PathBuf;

#[derive(thiserror::Error, Debug)]
pub enum FastqcError {
    #[error("Sequence format error: {0}")]
    SequenceFormat(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Configuration error: {0}")]
    Config(String),

    #[error("File not found: {0}")]
    FileNotFound(PathBuf),

    #[error("ZIP error: {0}")]
    Zip(#[from] zip::result::ZipError),

    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
}

pub type Result<T> = std::result::Result<T, FastqcError>;
