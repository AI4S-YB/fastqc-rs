pub mod fastq;
pub mod phred;
pub mod colorspace;
pub mod contaminant;
pub mod bam;
pub mod file_group;

/// Core sequence data structure, equivalent to Java's Sequence class.
#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub sequence: String,
    pub quality: String,
    pub colorspace: Option<String>,
    pub is_filtered: bool,
    pub file_name: String,
}

impl Sequence {
    pub fn new(id: String, sequence: String, quality: String, file_name: String) -> Self {
        Self {
            id,
            sequence: sequence.to_uppercase(),
            quality,
            colorspace: None,
            is_filtered: false,
            file_name,
        }
    }

    pub fn with_colorspace(
        id: String,
        sequence: String,
        colorspace: String,
        quality: String,
        file_name: String,
    ) -> Self {
        Self {
            id,
            sequence: sequence.to_uppercase(),
            quality,
            colorspace: Some(colorspace),
            is_filtered: false,
            file_name,
        }
    }
}

/// Trait for sequence file readers (FASTQ, BAM, etc.)
pub trait SequenceFile: Iterator<Item = crate::error::Result<Sequence>> {
    fn name(&self) -> &str;
    fn percent_complete(&self) -> u8;
    fn is_colorspace(&self) -> bool;
}
