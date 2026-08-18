use serde::Serialize;

use crate::target_list::{
    ValidTarget,
    chemistry::{EnsemblId, GeneName},
    target::{self, ValidGene},
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(transparent)]
pub struct XeniumPanelDesignerGeneList(Vec<XeniumPanelDesignerGene>);

impl XeniumPanelDesignerGeneList {
    pub fn from_valid_targets(mut valid_targets: Vec<ValidTarget>) -> Self {
        valid_targets.sort_by_key(|target| target.priority);

        Self(
            valid_targets
                .iter()
                .map(XeniumPanelDesignerGene::from_valid_target)
                .collect(),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
struct XeniumPanelDesignerGene {
    #[serde(rename = "Gene")]
    gene: GeneName,
    #[serde(rename = "Ensembl ID")]
    ensembl_id: EnsemblId,
    #[serde(rename = "Probe sets")]
    probe_sets: Option<u16>,
    #[serde(rename = "Force")]
    force: Option<Force>,
}

impl XeniumPanelDesignerGene {
    fn from_valid_target(
        ValidTarget {
            gene: ValidGene {
                ensembl_id,
                gene_name,
            },
            group: _,
            priority,
        }: &ValidTarget,
    ) -> Self {
        Self {
            gene: *gene_name,
            ensembl_id: *ensembl_id,
            probe_sets: None,
            force: (*priority == target::Priority::MustHave).then_some(Force::Forced),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum Force {
    Forced,
}
