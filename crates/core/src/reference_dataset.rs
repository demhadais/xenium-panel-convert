use std::{
    fs,
    path::{Path, PathBuf},
    str::FromStr,
};

use hdf5_metno::{
    File, H5Type,
    types::{FixedUnicode, VarLenUnicode},
};
use ndarray::{Array1, ArrayView};
use serde::Serialize;

use crate::reference_dataset::{
    self,
    obs::{read_cell_annotations_from_h5ad, read_cell_barcodes_from_h5ad},
    umi_counts::{RawCscUmiCounts, read_umi_counts_from_h5ad},
    var::{Features, read_features_from_h5ad},
};

mod common;
mod feature_sets;
mod obs;
mod umi_counts;
mod var;

pub fn validate_reference_dataset(
    path: impl AsRef<Path>,
    cell_barcode_col: &str,
    cell_annotation_col: &str,
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

    let barcodes = match read_cell_barcodes_from_h5ad(&file, cell_barcode_col) {
        Ok(b) => Some(b),
        Err(e) => {
            errors.push(e);
            None
        }
    };

    let cell_annotations = match read_cell_annotations_from_h5ad(&file, cell_annotation_col) {
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

    let (Some(counts), Some(barcodes), Some(cell_annotations), Some(features)) =
        (counts, barcodes, cell_annotations, features)
    else {
        return Err(errors);
    };

    Ok(ReferenceDataset {
        barcodes,
        counts,
        cell_annotations,
        features,
    })
}

// See https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices for format specifics
pub fn write_reference_dataset(
    path: impl AsRef<Path>,
    ReferenceDataset {
        barcodes,
        counts,
        cell_annotations,
        features,
    }: &ReferenceDataset,
) -> Result<(), Error> {
    fs::create_dir_all(&path).map_err(|e| Error::InvalidOutputPath {
        path: path.as_ref().to_str().map(str::to_owned).unwrap(),
        reason: e.to_string(),
    })?;

    let file = File::create_excl(path)?;

    write_dataset_to_h5(
        &file,
        "matrix/barcodes",
        barcodes.as_slice().ok_or(Error::Other)?,
    )?;
    write_dataset_to_h5(&file, "matrix/data", counts.data())?;
    write_dataset_to_h5(&file, "matrix/indices", counts.indices())?;
    write_dataset_to_h5(&file, "matrix/indptr", counts.indptr().iter().as_slice())?;
    write_dataset_to_h5(&file, "matrix/shape", &counts.shape())?;

    let feature_types: Vec<_> = features
        .feature_type()
        .into_iter()
        .map(|s| VarLenUnicode::from_str(s).unwrap())
        .collect();
    write_dataset_to_h5(&file, "features/feature_type", &feature_types)?;
    write_dataset_to_h5(&file, "features/id", features.id())?;
    write_dataset_to_h5(&file, "features/name", features.name())?;

    Ok(())
}

fn write_dataset_to_h5(file: &File, path: &str, data: &[impl H5Type]) -> Result<(), Error> {
    file.new_dataset_builder().with_data(data).create(path)?;

    Ok(())
}

pub struct ReferenceDataset {
    barcodes: Array1<VarLenUnicode>,
    counts: RawCscUmiCounts,
    cell_annotations: Array1<VarLenUnicode>,
    features: Features,
}

// TODO: don't need thiserror
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
    IncompleteFeatures { reason: String },
    #[error("unable to write files to {path}: {reason}")]
    InvalidOutputPath { path: String, reason: String },
    #[error("other")]
    Other,
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
