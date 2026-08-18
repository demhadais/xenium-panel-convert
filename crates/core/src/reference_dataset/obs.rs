use hdf5_metno::File;
use serde::Serialize;

use crate::reference_dataset::{
    Barcodes, CellAnnotations,
    columns::{CellAnnotationCol, CellBarcodeCol},
    h5_util::{self, ReadH5FieldError, read_1d_string_dataset, to_ascii},
};

// Read these into String because we will write them to a CSV, so we need serde
// support
pub(super) fn read_cell_annotations_from_h5ad(
    file: &File,
    annotation_col: &CellAnnotationCol,
) -> Result<CellAnnotations, ObsError> {
    let strings = read_1d_string_dataset(file, &format!("obs/{annotation_col}"))?;

    Ok(strings.mapv_into_any(|s| s.to_string()))
}

// If the dataset combines multiple smaller datasets, then the barcodes may
// contain sample names, which can be arbitrarily long. As such, we allow a
// conservative 64 bytes per barcode (18 bytes for the barcode itself and 46 for
// the sample name). This is easy to adjust should we find it to be too small
pub(super) fn read_cell_barcodes_from_h5ad(
    file: &File,
    barcode_col: &CellBarcodeCol,
) -> Result<Barcodes, ObsError> {
    let barcodes = h5_util::read_1d_string_dataset(file, &format!("obs/{barcode_col}"))?;

    Ok(barcodes.mapv_into_any(|b| to_ascii(&b)))
}

#[derive(Debug, Clone, Serialize, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum ObsError {
    #[error(transparent)]
    MalformedObs { error: ReadH5FieldError },
}

impl From<ReadH5FieldError> for ObsError {
    fn from(error: ReadH5FieldError) -> Self {
        Self::MalformedObs { error }
    }
}
