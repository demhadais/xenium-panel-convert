use std::collections::HashSet;

use hdf5_metno::{
    File,
    types::{FixedAscii, VarLenUnicode},
};
use ndarray::Array1;
use serde::Serialize;

use crate::reference_dataset::{
    columns::{EnsemblIdCol, GeneNameCol},
    feature_set::FeatureSet,
    h5_util::{read_1d_string_dataset, to_ascii},
};

pub(super) fn read_features_from_h5ad(
    file: &File,
    ensembl_id_col: &EnsemblIdCol,
    gene_name_col: &GeneNameCol,
    expected_feature_set: FeatureSet,
) -> Result<Features, Error> {
    let features = [
        ensembl_id_col.as_str(),
        gene_name_col.as_str(),
        "feature_types",
    ]
    .map(|s| format!("var/{s}"))
    .map(|path| read_1d_string_dataset(file, &path));

    let [Ok(ensembl_ids), Ok(gene_names), Ok(feature_types)] = features else {
        return Err(Error::IncompleteFeatures {
            reason: format!(
                "one of the columns '{ensembl_id_col}', '{gene_name_col}', or 'feature_types' \
                 were not found"
            ),
        });
    };

    check_feature_array_lens(&ensembl_ids, &gene_names, &feature_types)?;

    check_feature_set(
        &ensembl_ids,
        &gene_names,
        expected_feature_set.genes(ensembl_ids.len()),
    )?;

    Ok(Features {
        ensembl_ids: ensembl_ids.mapv_into_any(|s| to_ascii(&s)),
        gene_names: gene_names.mapv_into_any(|s| to_ascii(&s)),
        feature_types: feature_types.mapv_into_any(|s| to_ascii(&s)),
    })
}

fn check_feature_array_lens(
    ensembl_ids: &Array1<VarLenUnicode>,
    gene_names: &Array1<VarLenUnicode>,
    feature_types: &Array1<VarLenUnicode>,
) -> Result<(), Error> {
    if ensembl_ids.len() != gene_names.len() || ensembl_ids.len() != feature_types.len() {
        return Err(Error::InvalidShapes {
            ensembl_ids_len: ensembl_ids.len(),
            gene_names_len: gene_names.len(),
            feature_types_len: feature_types.len(),
        });
    }

    Ok(())
}

fn check_feature_set(
    ensembl_ids: &Array1<VarLenUnicode>,
    gene_names: &Array1<VarLenUnicode>,
    expected_features: Option<&phf::Map<&str, &str>>,
) -> Result<(), Error> {
    let err = Err(Error::UnexpectedFeatureSet {
        detail: "the set of features in the dataset must be the exact same as the unfiltered \
                 features of a cellranger output",
    });

    let Some(expected_features) = expected_features else {
        return err;
    };

    let mut seen = HashSet::with_capacity(ensembl_ids.len());

    for (id, name) in ensembl_ids.iter().zip(gene_names) {
        if !seen.insert(id) {
            return err;
        }

        let Some(expected_name) = expected_features.get(id) else {
            return err;
        };

        if name != *expected_name {
            return err;
        }
    }

    Ok(())
}

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
    ensembl_ids: EnsemblIds,
    gene_names: GeneNames,
    feature_types: FeatureTypes,
}

impl Features {
    pub(super) fn ensembl_ids(&self) -> &EnsemblIds {
        &self.ensembl_ids
    }

    pub(super) fn gene_names(&self) -> &GeneNames {
        &self.gene_names
    }

    pub(super) fn feature_types(&self) -> &FeatureTypes {
        &self.feature_types
    }

    pub(super) fn len(&self) -> usize {
        self.ensembl_ids.len()
    }
}

#[derive(Clone, thiserror::Error, Serialize, Debug)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Error {
    #[error("incomplete features: {reason}")]
    IncompleteFeatures { reason: String },
    #[error("unexpected feature-set")]
    UnexpectedFeatureSet { detail: &'static str },
    #[error(
        "ensembl_ids, gene_names, and feature_types should have the same length - found \
         ({ensembl_ids_len}, {gene_names_len}, {feature_types_len})"
    )]
    InvalidShapes {
        ensembl_ids_len: usize,
        gene_names_len: usize,
        feature_types_len: usize,
    },
}
