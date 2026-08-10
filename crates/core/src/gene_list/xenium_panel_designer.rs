use serde::Serialize;

use crate::gene_list::{
    TargetPriority, ValidGene, ValidTarget,
    chemistry::{EnsemblId, GeneName},
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(transparent)]
pub struct XeniumPanelDesignerGeneList(Vec<XeniumPanelDesignerTarget>);

impl XeniumPanelDesignerGeneList {
    pub fn from_valid_target_list(mut valid_targets: Vec<ValidTarget>) -> Self {
        valid_targets.sort_by_key(|target| target.priority);

        Self(
            valid_targets
                .iter()
                .map(XeniumPanelDesignerTarget::from_valid_target)
                .collect(),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
struct XeniumPanelDesignerTarget {
    #[serde(rename = "Gene")]
    gene: GeneName,
    #[serde(rename = "Ensembl ID")]
    ensembl_id: EnsemblId,
    #[serde(rename = "Probe sets")]
    probe_sets: Option<u16>,
    #[serde(rename = "Force")]
    force: Option<Force>,
}

impl XeniumPanelDesignerTarget {
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
            force: (*priority == TargetPriority::MustHave).then_some(Force::Forced),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
enum Force {
    Forced,
}
