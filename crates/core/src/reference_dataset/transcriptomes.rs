mod flex_refdata_gex_grch38_2020_a;
mod flex_refdata_gex_mm10_2020_a;
mod refdata_gex_grch38_2020_a;
mod refdata_gex_mm10_2020_a;

#[derive(Debug, Clone, Copy)]
pub struct Transcriptome {
    pub inner: TranscriptomeInner,
    pub flex: bool,
}

#[derive(Debug, Clone, Copy, strum::EnumString)]
pub enum TranscriptomeInner {
    #[strum(serialize = "GRCh38-2020-A", serialize = "h", ascii_case_insensitive)]
    Grch382020A,
    #[strum(serialize = "mm10-2020-A", serialize = "m", ascii_case_insensitive)]
    Mm102020A,
}
