use hdf5_metno::{File, types::VarLenUnicode};
use ndarray::Array1;

use crate::reference_dataset::{Error, common::read_string_array_from_file};

pub fn read_features_from_h5ad(
    file: &File,
    ensembl_id_col: &str,
    gene_name_col: &str,
) -> Result<Features, Error> {
    let features = [ensembl_id_col, gene_name_col, "feature_type"]
        .map(|s| format!("var/{s}"))
        .map(|s| read_string_array_from_file(file, &s));

    let [Ok(ensembl_ids), Ok(gene_names), Ok(feature_types)] = features else {
        return Err(Error::IncompleteFeatures {
            reason: "columns containing the Ensembl IDs, gene names, and feature types were not \
                     found in 'var'",
        });
    };

    Ok(Features {
        ensembl_ids,
        gene_names,
        feature_types,
    })
}

fn check_genes_are_from_correct_genome() {}

pub struct Features {
    ensembl_ids: Array1<VarLenUnicode>,
    gene_names: Array1<VarLenUnicode>,
    feature_types: Array1<VarLenUnicode>,
}
