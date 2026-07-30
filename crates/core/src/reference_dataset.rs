use std::path::Path;

use hdf5_metno::{
    File, H5Type,
    types::{FixedUnicode, VarLenUnicode},
};
use ndarray::Array1;
use serde::Serialize;

use crate::reference_dataset::{
    obs::read_cell_annotations_from_h5ad,
    umi_counts::{RawCscUmiCounts, read_umi_counts_from_h5ad},
    var::{Features, read_features_from_h5ad},
};

mod common;
mod obs;
mod umi_counts;
mod var;

pub fn validate_reference_dataset(
    path: impl AsRef<Path>,
    cell_annotations_col: &str,
    ensembl_id_col: &str,
    gene_name_col: &str,
) -> Result<ReferenceDataset, Vec<Error>> {
    let mut errors = Vec::new();

    let file = hdf5_metno::File::open(path).map_err(|e| vec![e.into()])?;
    let counts = match read_umi_counts_from_h5ad(&file) {
        Ok(c) => Some(c),
        Err(e) => {
            errors.push(e);
            None
        }
    };

    let annotations = match read_cell_annotations_from_h5ad(&file, cell_annotations_col) {
        Ok(a) => Some(a),
        Err(e) => {
            errors.push(e);
            None
        }
    };

    let features = match read_features_from_h5ad(&file, ensembl_id_col, gene_name_col) {
        Ok(f) => Some(f),
        Err(e) => {
            errors.push(e);
            None
        }
    };

    todo!()
}

pub struct ReferenceDataset {
    counts: RawCscUmiCounts,
    cell_annotations: Array1<VarLenUnicode>,
    features: Features,
}

#[derive(Debug, thiserror::Error, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Error {
    #[error("HDF5 error: {reason}")]
    Hdf5 { reason: String },
    #[error("empty UMI counts")]
    EmptyCounts,
    #[error("transformed UMI counts")]
    TransformedCounts,
    #[error("normalized UMI counts")]
    NormalizedCounts,
    #[error("malformed count matrix: {reason}")]
    MalformedMatrix { reason: String },
    #[error("unknown encoding type")]
    UnknownEncodingType,
    #[error("incomplete features: {reason}")]
    IncompleteFeatures { reason: &'static str },
}

impl From<hdf5_metno::Error> for Error {
    fn from(err: hdf5_metno::Error) -> Self {
        Self::Hdf5 {
            reason: err.to_string(),
        }
    }
}

impl From<sprs::errors::StructureError> for Error {
    fn from(err: sprs::errors::StructureError) -> Self {
        Self::MalformedMatrix {
            reason: err.to_string(),
        }
    }
}
