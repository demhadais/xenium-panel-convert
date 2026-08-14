use crate::reference_dataset::feature_set::{
    grch38_2020_a::GRCH38_2020_A, grch38_2020_a_flex::GRCH38_2020_A_FLEX,
    grch38_2024_a::GRCH38_2024_A, grch38_2024_a_flex_v1_1::GRCH38_2024_A_FLEX_V1_1,
    grch38_2024_a_flex_v2_0::GRCH38_2024_A_FLEX_V2_0, grcm39_2024_a::GRCM39_2024_A,
    grcm39_2024_a_flex_v1_1::GRCM39_2024_A_FLEX_V1_1,
    grcm39_2024_a_flex_v2_0::GRCM39_2024_A_FLEX_V2_0, mm10_2020_a::MM10_2020_A,
    mm10_2020_a_flex::MM10_2020_A_FLEX,
};

mod grch38_2020_a;
mod grch38_2020_a_flex;
mod grch38_2024_a;
mod grch38_2024_a_flex_v1_1;
mod grch38_2024_a_flex_v2_0;
mod grcm39_2024_a;
mod grcm39_2024_a_flex_v1_1;
mod grcm39_2024_a_flex_v2_0;
mod mm10_2020_a;
mod mm10_2020_a_flex;

#[derive(Debug, Clone, Copy)]
pub enum FeatureSet {
    ThreePrime(&'static phf::Map<&'static str, &'static str>),
    Flex2020A(&'static phf::Map<&'static str, &'static str>),
    Flex2024A {
        v1: &'static phf::Map<&'static str, &'static str>,
        v2: &'static phf::Map<&'static str, &'static str>,
    },
}

impl FeatureSet {
    pub fn new(transcriptome: Transcriptome, flex: bool) -> Self {
        match (transcriptome, flex) {
            (Transcriptome::Grch382020A, false) => Self::ThreePrime(&GRCH38_2020_A),
            (Transcriptome::Grch382020A, true) => Self::Flex2020A(&GRCH38_2020_A_FLEX),
            (Transcriptome::Grch382024A, false) => Self::ThreePrime(&GRCH38_2024_A),
            (Transcriptome::Grch382024A, true) => Self::Flex2024A {
                v1: &GRCH38_2024_A_FLEX_V1_1,
                v2: &GRCH38_2024_A_FLEX_V2_0,
            },
            (Transcriptome::Mm102020A, false) => Self::ThreePrime(&MM10_2020_A),
            (Transcriptome::Mm102020A, true) => Self::Flex2020A(&MM10_2020_A_FLEX),
            (Transcriptome::Grcm392024A, false) => Self::ThreePrime(&GRCM39_2024_A),
            (Transcriptome::Grcm392024A, true) => Self::Flex2024A {
                v1: &GRCM39_2024_A_FLEX_V1_1,
                v2: &GRCM39_2024_A_FLEX_V2_0,
            },
        }
    }

    pub fn genes(self, n_genes: usize) -> Option<&'static phf::Map<&'static str, &'static str>> {
        match self {
            Self::ThreePrime(genes) | Self::Flex2020A(genes) if n_genes == genes.len() => {
                Some(genes)
            }
            Self::Flex2024A { v1, v2 } => {
                let is_v1 = v1.len() == n_genes;
                let is_v2 = v2.len() == n_genes;

                let v1 = is_v1.then_some(v1);
                let v2 = is_v2.then_some(v2);

                // The only feature-sets with the same number of genes are probe-sets
                // GRCm39-2024-A v1.1.1 and probe-sets GRCm39-2024-A v2.0.0. Luckily, they are
                // exactly the same, so just default to v2 since that's more recent
                v2.or(v1)
            }
            Self::ThreePrime(_) | Self::Flex2020A(_) => None,
        }
    }
}

#[derive(Debug, Clone, Copy, strum::EnumString)]
pub enum Transcriptome {
    #[strum(
        serialize = "GRCh38-2020-A",
        serialize = "h2020",
        ascii_case_insensitive
    )]
    Grch382020A,
    #[strum(
        serialize = "GRCh38-2024-A",
        serialize = "h2024",
        ascii_case_insensitive
    )]
    Grch382024A,
    #[strum(serialize = "mm10-2020-A", serialize = "m2020", ascii_case_insensitive)]
    Mm102020A,
    #[strum(
        serialize = "GRCm39-2024-A",
        serialize = "m2024",
        ascii_case_insensitive
    )]
    Grcm392024A,
}
