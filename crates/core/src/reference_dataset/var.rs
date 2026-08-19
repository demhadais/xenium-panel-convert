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
    h5_util::{ReadH5FieldError, read_1d_string_dataset, to_ascii},
};

pub(super) fn read_features_from_h5ad(
    file: &File,
    ensembl_id_col: &EnsemblIdCol,
    gene_name_col: &GeneNameCol,
    expected_feature_set: FeatureSet,
) -> Result<Features, VarError> {
    let features = [
        ensembl_id_col.as_str(),
        gene_name_col.as_str(),
        "feature_types",
    ]
    .map(|s| format!("var/{s}"))
    .map(|path| read_1d_string_dataset(file, &path));

    let [Ok(ensembl_ids), Ok(gene_names), Ok(feature_types)] = features else {
        return Err(VarError::InvalidFields {
            errors: features
                .into_iter()
                .filter_map(std::result::Result::err)
                .collect(),
        });
    };

    check_feature_array_lens(&ensembl_ids, &gene_names, &feature_types)?;

    let Some(expected_genes) = expected_feature_set.genes(ensembl_ids.len()) else {
        return Err(VarError::FilteredGenes {
            n_expected_genes: expected_feature_set.n_genes(),
            n_found_genes: ensembl_ids.len(),
        });
    };

    check_feature_set(&ensembl_ids, &gene_names, expected_genes)?;

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
) -> Result<(), VarError> {
    if ensembl_ids.len() != gene_names.len() || ensembl_ids.len() != feature_types.len() {
        return Err(VarError::InvalidShapes {
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
    expected_features: &phf::Map<&str, &str>,
) -> Result<(), VarError> {
    let err = Err(VarError::UnexpectedFeatureSet {
        detail: "the set of features in the dataset must be the exact same as the unfiltered \
                 features of a cellranger output",
    });

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
pub(super) type EnsemblId = FixedAscii<18>;

type EnsemblIds = Array1<EnsemblId>;

// No gene name is likely to exceed 32 characters
pub(super) type GeneName = FixedAscii<32>;

type GeneNames = Array1<GeneName>;

type FeatureType = FixedAscii<32>;

type FeatureTypes = Array1<FeatureType>;

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

    pub(crate) fn len(&self) -> usize {
        self.ensembl_ids.len()
    }
}

#[derive(Clone, Serialize, Debug, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum VarError {
    #[error("one or more fields in var were improperly formatted - was scanpy used correctly?")]
    InvalidFields { errors: Vec<ReadH5FieldError> },
    #[error("unexpected feature set: {detail} - were any genes filtered?")]
    UnexpectedFeatureSet { detail: &'static str },
    #[error("some genes were filtered out of the dataset - expected: {}, found: {n_found_genes}", n_expected_genes.1.map(|n1| format!("{} or {n1}", n_expected_genes.0)).unwrap_or(n_expected_genes.0.to_string()))]
    FilteredGenes {
        n_expected_genes: (usize, Option<usize>),
        n_found_genes: usize,
    },
    #[error(
        "invalid shapes: {ensembl_ids_len} Ensembl IDs, {gene_names_len} gene names, \
         {feature_types_len} feature types - expected them all to be the same"
    )]
    InvalidShapes {
        ensembl_ids_len: usize,
        gene_names_len: usize,
        feature_types_len: usize,
    },
}

#[cfg(test)]
mod tests {
    use hdf5_metno::File;

    use crate::reference_dataset::{
        columns::{EnsemblIdCol, GeneNameCol},
        feature_set::{FeatureSet, Transcriptome},
        var::{VarError, read_features_from_h5ad},
    };

    fn read_generated_features(
        ensembl_id_col: &str,
        gene_name_col: &str,
    ) -> Result<super::Features, VarError> {
        read_features_from_h5ad(
            &File::open("test-data/csr_adata.h5ad").unwrap(),
            &EnsemblIdCol(ensembl_id_col.to_owned()),
            &GeneNameCol(gene_name_col.to_owned()),
            FeatureSet::new(Transcriptome::Grch382020A, false),
        )
    }

    #[test]
    fn filtered_genes_are_rejected() {
        // The generated datasets have 100 genes
        let err = read_generated_features("ensembl_id", "gene_name").unwrap_err();

        std::assert_matches!(
            err,
            VarError::FilteredGenes {
                n_found_genes: 100,
                ..
            }
        );
    }

    #[test]
    fn missing_var_columns_are_collected() {
        let err = read_generated_features("nonexistent", "also_nonexistent").unwrap_err();

        let VarError::InvalidFields { errors } = err else {
            panic!("expected invalid fields, got {err:?}");
        };

        assert_eq!(
            errors.len(),
            2,
            "both missing columns in var should be reported"
        );
    }
}
