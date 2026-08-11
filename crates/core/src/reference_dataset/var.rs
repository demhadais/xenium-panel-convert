use hdf5_metno::{File, types::FixedAscii};
use ndarray::Array1;
use serde::Serialize;

use crate::reference_dataset::h5_util::{read_1d_string_dataset, to_ascii};

pub(crate) fn read_features_from_h5ad(
    file: &File,
    ensembl_id_col: &str,
    gene_name_col: &str,
) -> Result<Features, Error> {
    let features = [ensembl_id_col, gene_name_col, "feature_types"]
        .map(|s| format!("var/{s}"))
        .map(|path| read_1d_string_dataset(file, &path));

    let [Ok(ensembl_ids), Ok(gene_names), Ok(feature_type)] = features else {
        return Err(Error::IncompleteFeatures {
            reason: format!(
                "one of the columns '{ensembl_id_col}', '{gene_name_col}', or 'feature_types' \
                 were not found"
            ),
        });
    };

    Ok(Features {
        id: ensembl_ids.mapv_into_any(|s| to_ascii(&s)),
        name: gene_names.mapv_into_any(|s| to_ascii(&s)),
        feature_type: feature_type.mapv_into_any(|s| to_ascii(&s)),
    })
}

// Placeholder for the planned wrong-reference-genome check (see docs/GOALS.md)
#[allow(dead_code)]
fn check_genes_are_from_correct_genome() {}

// Human Ensembl IDs are 15 characters while mouse Ensembl IDs are 18
pub(crate) type EnsemblId = FixedAscii<18>;

pub(crate) type EnsemblIds = Array1<EnsemblId>;

// No gene name is likely to exceed 32 characters
pub(crate) type GeneName = FixedAscii<32>;

pub(crate) type GeneNames = Array1<GeneName>;

pub(crate) type FeatureType = FixedAscii<32>;

pub(crate) type FeatureTypes = Array1<FeatureType>;

#[derive(Debug, PartialEq)]
pub(crate) struct Features {
    pub(crate) id: EnsemblIds,
    pub(crate) name: GeneNames,
    pub(crate) feature_type: FeatureTypes,
}

#[derive(Clone, thiserror::Error, Serialize, Debug)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Error {
    #[error("incomplete features: {reason}")]
    IncompleteFeatures { reason: String },
}
