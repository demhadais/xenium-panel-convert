use hdf5_metno::File;
use serde::Serialize;

use crate::reference_dataset::{
    Barcodes, CellAnnotations,
    h5::{self, ReadFieldError, read_1d_string_dataset, to_ascii},
};

// Read these into String because we will write them to a CSV, so we need serde support
pub fn read_cell_annotations_from_h5ad(
    file: &File,
    annotations_col: &str,
) -> Result<CellAnnotations, Error> {
    let strings = read_1d_string_dataset(file, &format!("obs/{annotations_col}"))?;

    Ok(strings.mapv_into_any(|s| s.to_string()))
}

// If the dataset combines multiple smaller datasets, then the barcodes may contain sample names, which can be arbitrarily long. As such, we allow a conservative 64 bytes per barcode (18 bytes for the barcode itself and 46 for the sample name). This is easy to adjust should we find it to be too small
pub fn read_cell_barcodes_from_h5ad(file: &File, barcodes_col: &str) -> Result<Barcodes, Error> {
    let barcodes = h5::read_1d_string_dataset(file, &format!("obs/{barcodes_col}"))?;

    Ok(barcodes.mapv_into_any(|b| to_ascii(&b)))
}

#[derive(Debug, Clone, thiserror::Error, Serialize)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum Error {
    #[error(transparent)]
    MalformedObs {
        #[from]
        error: ReadFieldError,
    },
}
