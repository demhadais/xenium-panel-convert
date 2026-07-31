use hdf5_metno::{File, types::VarLenUnicode};
use ndarray::Array1;

use crate::reference_dataset::{Error, common::read_string_array_from_file};

pub fn read_cell_annotations_from_h5ad(
    file: &File,
    annotations_col: &str,
) -> Result<Array1<VarLenUnicode>, Error> {
    let annotations = read_string_array_from_file(file, &format!("obs/{annotations_col}"))?;

    Ok(annotations)
}

pub fn read_cell_barcodes_from_h5ad(
    file: &File,
    barcodes_col: &str,
) -> Result<Array1<VarLenUnicode>, Error> {
    let barcodes = read_string_array_from_file(file, &format!("obs/{barcodes_col}"))?;

    Ok(barcodes)
}
