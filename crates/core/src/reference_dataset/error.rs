use camino::Utf8PathBuf;
use serde::Serialize;

use crate::{
    common::ErrorVecExt,
    reference_dataset::{
        h5_util::{CreateH5GroupError, WriteH5DatasetError},
        obs::ObsError,
        pseudo_anndata::ShapeMismatchError,
        umi_counts::UmiCountsError,
        var::VarError,
    },
};

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ReadReferenceDatasetError {
    #[error("ensure the H5AD is properly formatted")]
    InvalidH5File { reason: String },
    #[error(transparent)]
    UmiCounts { error: UmiCountsError },
    #[error(transparent)]
    Obs { error: ObsError },
    #[error(transparent)]
    Var { error: VarError },
    #[error(transparent)]
    Shape { error: ShapeMismatchError },
}

impl From<UmiCountsError> for ReadReferenceDatasetError {
    fn from(error: UmiCountsError) -> Self {
        Self::UmiCounts { error }
    }
}

impl From<ObsError> for ReadReferenceDatasetError {
    fn from(error: ObsError) -> Self {
        Self::Obs { error }
    }
}

impl From<VarError> for ReadReferenceDatasetError {
    fn from(error: VarError) -> Self {
        Self::Var { error }
    }
}

impl From<ShapeMismatchError> for ReadReferenceDatasetError {
    fn from(error: ShapeMismatchError) -> Self {
        Self::Shape { error }
    }
}

impl<E> ErrorVecExt<E> for Vec<ReadReferenceDatasetError>
where
    E: Into<ReadReferenceDatasetError>,
{
    fn push_err<T>(&mut self, err: E) -> Option<T> {
        self.push(err.into());

        None
    }
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum WriteReferenceDatasetError {
    #[error("failed to create output directory - {reason}")]
    CreateOutputDir { path: Utf8PathBuf, reason: String },
    #[error("failed to create matrix.h5 - {reason}")]
    CreateMatrixFile { path: Utf8PathBuf, reason: String },
    #[error("{error}")]
    CreateH5Group {
        path: Utf8PathBuf,
        error: CreateH5GroupError,
    },
    #[error("{error}")]
    WriteH5Dataset {
        path: Utf8PathBuf,
        error: WriteH5DatasetError,
    },
    #[error("cannot overwrite {path} - move or delete the existing annotation.csv file")]
    AnnotationsCsvExists { path: Utf8PathBuf },
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[error("{error}\nhint: {hint}")]
pub struct ReadReferenceDatasetErrorWrapper {
    pub error: ReadReferenceDatasetError,
    pub hint: String,
}

impl From<ReadReferenceDatasetError> for ReadReferenceDatasetErrorWrapper {
    fn from(error: ReadReferenceDatasetError) -> Self {
        Self {
            hint: error.to_string(),
            error,
        }
    }
}

#[derive(Clone, Debug, Serialize, thiserror::Error)]
#[error("{error}\nhint: {hint}")]
pub struct WriteReferenceDatasetErrorWrapper {
    pub error: WriteReferenceDatasetError,
    pub hint: String,
}

impl From<WriteReferenceDatasetError> for WriteReferenceDatasetErrorWrapper {
    fn from(error: WriteReferenceDatasetError) -> Self {
        Self {
            hint: error.to_string(),
            error,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ReadReferenceDatasetErrors {
    pub path: Utf8PathBuf,
    pub errors: Vec<ReadReferenceDatasetErrorWrapper>,
}
