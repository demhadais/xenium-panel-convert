use crate::{
    reference_dataset::{
        feature_set::{FeatureSet, Transcriptome},
        pseudo_anndata::PseudoAnndata,
    },
    target_list::{
        chemistry::Species, target::ValidGene, xenium_panel_designer::XeniumPanelDesignerGeneList,
    },
};

#[must_use]
pub fn validate_target_list_and_reference_dataset_compatibility(
    target_list: &XeniumPanelDesignerGeneList,
    target_list_species: Species,
    reference_dataset: &PseudoAnndata,
    reference_dataset_transcriptome: Transcriptome,
    reference_dataset_is_flex: bool,
) -> Vec<TargetListReferenceDatasetCompatibilityWarning> {
    // Return early because this warning implies the other two warning types
    if !species_and_transcriptome_match(target_list_species, reference_dataset_transcriptome) {
        return vec![
            TargetListReferenceDatasetCompatibilityWarning::SpeciesTranscriptomeMismatch {
                target_list_species,
                reference_dataset_transcriptome,
            },
        ];
    }

    let mut warnings = Vec::with_capacity(target_list.len());

    for target in target_list.as_slice() {
        match validate_gene_is_in_feature_set_with_correct_name(
            target.gene(),
            reference_dataset,
            reference_dataset_transcriptome,
            reference_dataset_is_flex,
        ) {
            Ok(()) => (),
            Err(w) => {
                warnings.push(w);
            }
        }
    }

    warnings
}

fn species_and_transcriptome_match(species: Species, transcriptome: Transcriptome) -> bool {
    matches!(
        (species, transcriptome),
        (
            Species::HomoSapiens,
            Transcriptome::Grch382020A | Transcriptome::Grch382024A
        ) | (
            Species::MusMusculus,
            Transcriptome::Mm102020A | Transcriptome::Grcm392024A
        )
    )
}

fn validate_gene_is_in_feature_set_with_correct_name(
    target: ValidGene,
    reference_dataset: &PseudoAnndata,
    reference_dataset_transcriptome: Transcriptome,
    reference_dataset_is_flex: bool,
) -> Result<(), TargetListReferenceDatasetCompatibilityWarning> {
    let feature_set = FeatureSet::new(reference_dataset_transcriptome, reference_dataset_is_flex);
    let feature_set = feature_set
        .genes(reference_dataset.features().len())
        .expect("if we have a PseudoAnndata, we know its features are exactly the transcriptome");

    let gene_name_from_transcriptome = feature_set.get(target.ensembl_id.as_str()).ok_or(
        TargetListReferenceDatasetCompatibilityWarning::TargetNotInReferenceDataset {
            gene: target,
            transcriptome: reference_dataset_transcriptome,
        },
    )?;

    if target.gene_name != *gene_name_from_transcriptome {
        return Err(
            TargetListReferenceDatasetCompatibilityWarning::GeneNameMismatch {
                gene_in_target_list: target,
                gene_name_in_reference_dataset: gene_name_from_transcriptome,
            },
        );
    }

    Ok(())
}

#[derive(Clone, Copy, Debug, thiserror::Error, serde::Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum TargetListReferenceDatasetCompatibilityWarning {
    #[error(
        "target-list and reference dataset transcriptome do not have the same species \
         ({target_list_species} and {reference_dataset_transcriptome})"
    )]
    SpeciesTranscriptomeMismatch {
        target_list_species: Species,
        reference_dataset_transcriptome: Transcriptome,
    },
    #[error("{} ({}) not in reference dataset, whose transcriptome is {transcriptome}", gene.gene_name, gene.ensembl_id)]
    TargetNotInReferenceDataset {
        gene: ValidGene,
        transcriptome: Transcriptome,
    },
    #[error("{} is called {} in the target-list, but it is called {gene_name_in_reference_dataset} in the reference dataset - this is likely because the reference-dataset was aligned against GRCh38-2024-A or GRCm39-2024-A - add a new column to .var in the reference dataset with gene names from the 2020-A version of the transcriptome", gene_in_target_list.ensembl_id, gene_in_target_list.gene_name)]
    GeneNameMismatch {
        gene_in_target_list: ValidGene,
        gene_name_in_reference_dataset: &'static str,
    },
}
