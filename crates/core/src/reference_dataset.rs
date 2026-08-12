use std::{fs, path::Path};

use hdf5_metno::{File, types::FixedAscii};
use ndarray::Array1;
use serde::Serialize;

use crate::{
    Species,
    common::ErrorVecExt,
    reference_dataset::{
        columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneNameCol},
        h5_util::write_dataset_to_h5_group,
        obs::{read_cell_annotations_from_h5ad, read_cell_barcodes_from_h5ad},
        pseudo_anndata::PseudoAnndata,
        umi_counts::read_umi_counts_from_h5ad,
        var::{Features, read_features_from_h5ad},
    },
};

pub mod columns;
pub mod feature_sets;
pub mod h5_util;
pub mod obs;
mod pseudo_anndata;
pub mod specification;
pub mod umi_counts;
pub mod var;

pub fn read_reference_dataset(
    path: impl AsRef<Path>,
    cell_barcode_col: &CellBarcodeCol,
    cell_annotation_col: &CellAnnotationCol,
    ensembl_id_col: &EnsemblIdCol,
    gene_name_col: &GeneNameCol,
    species: Species,
) -> Result<PseudoAnndata, Vec<Error>> {
    let mut errors = Vec::new();

    let file = hdf5_metno::File::open(path).map_err(|e| {
        vec![Error::InvalidH5File {
            reason: e.to_string(),
        }]
    })?;

    let counts = read_umi_counts_from_h5ad(&file).map_or_else(|err| errors.push_err(err), Some);

    let barcodes = read_cell_barcodes_from_h5ad(&file, cell_barcode_col)
        .map_or_else(|err| errors.push_err(err), Some);

    let cell_annotations = read_cell_annotations_from_h5ad(&file, cell_annotation_col)
        .map_or_else(|err| errors.push_err(err), Some);

    let features = read_features_from_h5ad(&file, ensembl_id_col, gene_name_col, todo!())
        .map_or_else(|err| errors.push_err(err), Some);

    let (Some(counts), Some(barcodes), Some(cell_annotations), Some(features)) =
        (counts, barcodes, cell_annotations, features)
    else {
        return Err(errors);
    };

    let anndata =
        PseudoAnndata::new(counts, barcodes, cell_annotations, features).map_err(|err| {
            errors.push_err::<()>(err);
            errors
        })?;

    Ok(anndata)
}

// See https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices for format specifics
pub fn write_reference_dataset(dir: impl AsRef<Path>, ds: &PseudoAnndata) -> Result<(), Error> {
    let dir = dir.as_ref();

    if !dir.exists() {
        fs::create_dir_all(dir).map_err(|e| Error::from_path_error(dir, e))?;
    }

    let matrix_path = dir.join("matrix.h5");
    let file =
        File::create_excl(&matrix_path).map_err(|e| Error::from_path_error(&matrix_path, e))?;

    let matrix_group = file.create_group("matrix").unwrap();

    let barcodes = ds.barcodes();
    write_dataset_to_h5_group(&matrix_group, "barcodes", barcodes)?;

    let counts = ds.counts();
    write_dataset_to_h5_group(&matrix_group, "data", counts.data())?;
    write_dataset_to_h5_group(&matrix_group, "indices", counts.indices())?;
    write_dataset_to_h5_group(&matrix_group, "indptr", counts.indptr().iter().as_slice())?;
    write_dataset_to_h5_group(&matrix_group, "shape", &counts.shape())?;

    let features = ds.features();
    write_dataset_to_h5_group(
        &matrix_group,
        "features/feature_type",
        features.feature_types(),
    )?;
    write_dataset_to_h5_group(&matrix_group, "features/id", features.ensembl_ids())?;
    write_dataset_to_h5_group(&matrix_group, "features/name", features.gene_names())?;

    write_annotations_csv(dir.join("annotations.csv"), barcodes, ds.cell_annotations())?;

    Ok(())
}

fn write_annotations_csv(
    path: impl AsRef<Path>,
    barcodes: &Barcodes,
    annotations: &CellAnnotations,
) -> Result<(), Error> {
    #[derive(Debug, Serialize)]
    struct CellAnnotation<'a> {
        barcode: &'a str,
        annotation: &'a str,
    }

    let mut writer = csv::Writer::from_path(path).unwrap();
    for (barcode, annotation) in barcodes.iter().zip(annotations) {
        writer
            .serialize(CellAnnotation {
                barcode,
                annotation,
            })
            .unwrap();
    }

    Ok(())
}

type Barcode = FixedAscii<64>;

type Barcodes = Array1<Barcode>;

type CellAnnotations = Array1<String>;

#[derive(Clone, Debug, thiserror::Error, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Error {
    #[error("invalid H5 file: {reason}")]
    InvalidH5File { reason: String },
    #[error(transparent)]
    UmiCounts {
        #[from]
        error: umi_counts::Error,
    },
    #[error(transparent)]
    Obs {
        #[from]
        error: obs::Error,
    },
    #[error(transparent)]
    Var {
        #[from]
        error: var::Error,
    },
    #[error("number of cell barcodes != number of cell annotations")]
    Shape {
        n_barcodes: usize,
        n_annotations: usize,
        n_features: usize,
        counts_shape: [i32; 2],
    },
    #[error("invalid output path: {reason}")]
    InvalidOutputPath { path: String, reason: String },
    #[error("something went wrong writing H5 object: {reason}")]
    WriteH5Object { path: String, reason: String },
}

impl<E> ErrorVecExt<E> for Vec<Error>
where
    E: Into<Error>,
{
    fn push_err<T>(&mut self, err: E) -> Option<T> {
        self.push(err.into());

        None
    }
}

impl Error {
    fn from_path_error(path: impl AsRef<Path>, error: impl std::error::Error) -> Self {
        Self::InvalidOutputPath {
            path: path.as_ref().to_str().map(str::to_owned).unwrap(),
            reason: error.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use anyhow::Context;
    use hdf5_metno::types::FixedAscii;

    use crate::{
        Species,
        reference_dataset::{
            h5_util::read_1d_dataset,
            read_reference_dataset,
            var::{EnsemblId, GeneName},
        },
    };

    #[test]
    fn read_generated_h5ad_files() {
        let paths = [
            "test-data/csr_adata.h5ad",
            "test-data/csc_adata.h5ad",
            "test-data/dense_adata.h5ad",
        ]
        .map(|f| Path::new(f));

        let mut fake_counts = Vec::with_capacity(paths.len());

        for path in paths {
            let filename = path.to_str().unwrap();
            let reference_dataset = read_reference_dataset(
                path,
                &"barcode".into(),
                &"annotation".into(),
                &"ensembl_id".into(),
                &"gene_name".into(),
                Species::HomoSapiens,
            )
            .map_err(|e| e[0].clone())
            .context(format!("failed to read dataset from {filename}"))
            .unwrap();

            assert_eq!(
                reference_dataset.counts().data()[0],
                10,
                "first entry in UMI counts of {filename} != 10"
            );

            fake_counts.push(reference_dataset.counts().to_owned());
        }

        assert_eq!(fake_counts[0], fake_counts[1]);
        assert_eq!(fake_counts[0], fake_counts[2]);
    }

    #[test]
    fn read_real_h5ad() {
        let filename = "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5ad";

        let reconstructed_dataset = read_reference_dataset(
            filename,
            &"barcode".into(),
            &"annotation".into(),
            &"gene_ids".into(),
            &"gene_name".into(),
            Species::HomoSapiens,
        )
        .map_err(|e| e[0].clone())
        .context(format!("failed to validate {filename}"))
        .unwrap();

        // Compare it against the original data
        let original_h5 = hdf5_metno::File::open(
            "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5",
        )
        .unwrap();

        let original_counts = read_1d_dataset::<i32>(&original_h5, "matrix/data").unwrap();
        assert_eq!(
            original_counts.as_slice().unwrap(),
            reconstructed_dataset.counts().data(),
            "UMI counts were not correctly reconstructed"
        );

        let original_indices = read_1d_dataset::<i64>(&original_h5, "matrix/indices").unwrap();
        assert_eq!(
            original_indices.as_slice().unwrap(),
            reconstructed_dataset.counts().indices(),
            "UMI count indices were not correctly reconstructed"
        );

        let original_indptr = read_1d_dataset::<i64>(&original_h5, "matrix/indptr").unwrap();
        assert_eq!(
            original_indptr.as_slice().unwrap(),
            reconstructed_dataset.counts().indptr().iter().as_slice(),
            "UMI counts indptr was not correctly reconstructed"
        );

        let original_barcodes =
            read_1d_dataset::<FixedAscii<64>>(&original_h5, "matrix/barcodes").unwrap();
        assert_eq!(original_barcodes, reconstructed_dataset.barcodes());

        let original_feature_ids =
            read_1d_dataset::<EnsemblId>(&original_h5, "matrix/features/id").unwrap();
        assert_eq!(
            original_feature_ids,
            reconstructed_dataset.features().ensembl_ids()
        );

        let original_feature_names =
            read_1d_dataset::<GeneName>(&original_h5, "matrix/features/name").unwrap();
        assert_eq!(
            original_feature_names,
            reconstructed_dataset.features().gene_names()
        );
    }
}
