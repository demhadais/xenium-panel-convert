use hdf5_metno::{File, types::VarLenUnicode};
use ndarray::Array1;

use crate::reference_dataset::{Error, common::read_string_array_from_file};

pub fn read_features_from_h5ad(
    file: &File,
    ensembl_id_col: &str,
    gene_name_col: &str,
) -> Result<Features, Error> {
    let features = [ensembl_id_col, gene_name_col]
        .map(|s| format!("var/{s}"))
        .map(|s| read_string_array_from_file(file, &s));

    let [Ok(ensembl_ids), Ok(gene_names)] = features else {
        return Err(Error::IncompleteFeatures {
            reason: format!("columns {ensembl_id_col} and/or {gene_name_col} were not found"),
        });
    };

    let shape = ensembl_ids.raw_dim();

    Ok(Features {
        id: ensembl_ids,
        name: gene_names,
        feature_type: Array1::from_elem(shape, "Gene Expression"),
    })
}

fn check_genes_are_from_correct_genome() {}

pub struct Features {
    id: Array1<VarLenUnicode>,
    name: Array1<VarLenUnicode>,
    feature_type: Array1<&'static str>,
}

impl Features {
    pub fn id(&self) -> &[VarLenUnicode] {
        &self.id.as_slice().unwrap()
    }

    pub fn name(&self) -> &[VarLenUnicode] {
        &self.name.as_slice().unwrap()
    }

    pub fn feature_type(&self) -> &[&'static str] {
        &self.feature_type.as_slice().unwrap()
    }
}
