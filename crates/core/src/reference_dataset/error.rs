use camino::Utf8Path;
use serde::Serialize;

use crate::{
    common::ErrorVecExt,
    reference_dataset::{obs::ObsError, umi_counts::UmiCountsError, var::VarError},
};

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ReferenceDatasetError {
    #[error("invalid H5 file - {reason}")]
    InvalidH5File {
        reason: String,
    },
    #[error(transparent)]
    UmiCounts(UmiCountsError),
    #[error(transparent)]
    Obs(ObsError),
    #[error(transparent)]
    Var(VarError),
    #[error(
        "invalid shape - {n_barcodes} barcodes, {n_annotations} cell annotations, {n_features} \
         features, counts shape {counts_shape:?}"
    )]
    Shape {
        n_barcodes: usize,
        n_annotations: usize,
        n_features: usize,
        counts_shape: [i32; 2],
    },
    #[error("invalid output path {path} - {reason}")]
    InvalidOutputPath {
        path: String,
        reason: String,
    },
    #[error("failed to write H5 object {path} - {reason}")]
    WriteH5Object {
        path: String,
        reason: String,
    },
}

impl ReferenceDatasetError {
    pub(super) fn from_path_error(path: &Utf8Path, error: impl std::error::Error) -> Self {
        Self::InvalidOutputPath {
            path: path.to_string(),
            reason: error.to_string(),
        }
    }
}

impl From<UmiCountsError> for ReferenceDatasetError {
    fn from(value: UmiCountsError) -> Self {
        Self::UmiCounts(value)
    }
}

impl From<ObsError> for ReferenceDatasetError {
    fn from(value: ObsError) -> Self {
        Self::Obs(value)
    }
}

impl From<VarError> for ReferenceDatasetError {
    fn from(value: VarError) -> Self {
        Self::Var(value)
    }
}

impl<E> ErrorVecExt<E> for Vec<ReferenceDatasetError>
where
    E: Into<ReferenceDatasetError>,
{
    fn push_err<T>(&mut self, err: E) -> Option<T> {
        self.push(err.into());

        None
    }
}
