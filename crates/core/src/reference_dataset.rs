use std::path::Path;

use serde::Serialize;

use crate::reference_dataset::umi_counts::{RawUmiCounts, read_umi_counts_from_h5ad};

mod obs;
mod umi_counts;
mod var;

pub fn validate_reference_dataset(
    path: impl AsRef<Path>,
    cell_annotations_col: &str,
    ensembl_id_col: &str,
    gene_name_col: &str,
) -> Result<ReferenceDataset, Vec<Error>> {
    let mut errors = vec![];

    let file = hdf5_metno::File::open(path).map_err(|e| vec![e.into()])?;

    let counts = match read_umi_counts_from_h5ad(&file) {
        Ok(counts) => Some(counts),
        Err(e) => {
            errors.push(e.into());
            None
        }
    };

    let Some(counts) = counts else {
        return Err(errors);
    };

    Ok(ReferenceDataset { counts })
}

pub struct ReferenceDataset {
    counts: RawUmiCounts,
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
