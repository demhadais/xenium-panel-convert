use std::path::Path;

use serde::Serialize;

use crate::reference_dataset::umi_counts::{RawCscUmiCounts, read_umi_counts_from_h5ad};

mod obs;
mod umi_counts;
mod var;

pub fn validate_reference_dataset(
    path: impl AsRef<Path>,
    _cell_annotations_col: &str,
    _ensembl_id_col: &str,
    _gene_name_col: &str,
) -> Result<ReferenceDataset, Vec<Error>> {
    let file = hdf5_metno::File::open(path).map_err(|e| vec![e.into()])?;
    let counts = read_umi_counts_from_h5ad(&file).map_err(|e| vec![e.into()])?;

    Ok(ReferenceDataset { counts })
}

pub struct ReferenceDataset {
    counts: RawCscUmiCounts,
}

#[derive(Debug, thiserror::Error, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Error {
    #[error("HDF5 error: {reason}")]
    Hdf5 { reason: String },
    #[error(transparent)]
    #[serde(untagged)]
    Counts(#[from] umi_counts::Error),
}

impl From<hdf5_metno::Error> for Error {
    fn from(err: hdf5_metno::Error) -> Self {
        Self::Hdf5 {
            reason: err.to_string(),
        }
    }
}
