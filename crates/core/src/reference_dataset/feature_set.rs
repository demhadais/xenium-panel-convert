use crate::reference_dataset::feature_set::{
    grch38_2020_a::GRCH38_2020_A, grch38_2020_a_flex::GRCH38_2020_A_FLEX, mm10_2020_a::MM10,
    mm10_2020_a_flex::MM10_FLEX,
};

mod grch38_2020_a;
mod grch38_2020_a_flex;
mod mm10_2020_a;
mod mm10_2020_a_flex;

#[derive(Debug, Clone, Copy)]
pub struct FeatureSet {
    genes: &'static phf::Map<&'static str, &'static str>,
}

impl FeatureSet {
    pub fn new(transcriptome: Transcriptome, flex: bool) -> Self {
        match (transcriptome, flex) {
            (Transcriptome::Grch382020A, false) => Self {
                genes: &GRCH38_2020_A,
            },
            (Transcriptome::Grch382020A, true) => Self {
                genes: &GRCH38_2020_A_FLEX,
            },
            (Transcriptome::Mm102020A, false) => Self { genes: &MM10 },
            (Transcriptome::Mm102020A, true) => Self { genes: &MM10_FLEX },
        }
    }

    pub fn reference_transcriptome(&self) -> &'static phf::Map<&'static str, &'static str> {
        self.genes
    }
}

#[derive(Debug, Clone, Copy, strum::EnumString)]
pub enum Transcriptome {
    #[strum(serialize = "GRCh38-2020-A", serialize = "h", ascii_case_insensitive)]
    Grch382020A,
    #[strum(serialize = "mm10-2020-A", serialize = "m", ascii_case_insensitive)]
    Mm102020A,
}
