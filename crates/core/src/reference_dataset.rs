use std::{fs, path::Path};

use hdf5_metno::{File, Group, H5Type, types::FixedAscii};
use ndarray::{Array1, ArrayView};
use serde::Serialize;

use crate::reference_dataset::{
    obs::{read_cell_annotations_from_h5ad, read_cell_barcodes_from_h5ad},
    umi_counts::{RawCscUmiCounts, read_umi_counts_from_h5ad},
    var::{Features, read_features_from_h5ad},
};

mod feature_sets;
mod h5;
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

    let file = hdf5_metno::File::open(path).map_err(|e| {
        vec![Error::InvalidH5File {
            reason: e.to_string(),
        }]
    })?;

    let counts = match read_umi_counts_from_h5ad(&file) {
        Ok(c) => Some(c),
        Err(error) => {
            errors.push(Error::UmiCounts { error });
            None
        }
    };

    let barcodes = match read_cell_barcodes_from_h5ad(&file, cell_barcode_col) {
        Ok(b) => Some(b),
        Err(error) => {
            errors.push(Error::Obs { error });
            None
        }
    };

    let cell_annotations = match read_cell_annotations_from_h5ad(&file, cell_annotation_col) {
        Ok(a) => Some(a),
        Err(error) => {
            errors.push(Error::Obs { error });
            None
        }
    };

    let features = match read_features_from_h5ad(&file, ensembl_id_col, gene_name_col) {
        Ok(f) => Some(f),
        Err(error) => {
            errors.push(Error::Var { error });
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
    let path = path.as_ref();

    fs::create_dir_all(path).map_err(|e| Error::from_path_error(path, e))?;

    let matrix_path = path.join("matrix.h5");
    let file =
        File::create_excl(&matrix_path).map_err(|e| Error::from_path_error(&matrix_path, e))?;

    let matrix_group = file.create_group("matrix").unwrap();

    write_dataset_to_h5_group(&matrix_group, "barcodes", barcodes)?;
    write_dataset_to_h5_group(&matrix_group, "data", counts.data())?;
    write_dataset_to_h5_group(&matrix_group, "indices", counts.indices())?;
    write_dataset_to_h5_group(&matrix_group, "indptr", counts.indptr().iter().as_slice())?;
    write_dataset_to_h5_group(&matrix_group, "shape", &counts.shape())?;

    write_dataset_to_h5_group(
        &matrix_group,
        "features/feature_type",
        features.feature_type(),
    )?;
    write_dataset_to_h5_group(&matrix_group, "features/id", features.id())?;
    write_dataset_to_h5_group(&matrix_group, "features/name", features.name())?;

    write_annotations_csv(path, barcodes, cell_annotations)?;

    Ok(())
}

fn write_dataset_to_h5_group<'d, A, T, D>(file: &Group, path: &str, data: A) -> Result<(), Error>
where
    A: Into<ArrayView<'d, T, D>>,
    T: H5Type,
    D: ndarray::Dimension,
{
    file.new_dataset_builder()
        .with_data(data)
        .create(path)
        .map_err(|e| Error::WriteH5Object {
            path: path.to_owned(),
            reason: e.to_string(),
        })?;

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

    if barcodes.len() != annotations.len() {
        return Err(Error::BarcodeAnnotationMismatch);
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

#[derive(Debug, PartialEq)]
pub struct ReferenceDataset {
    barcodes: Barcodes,
    counts: RawCscUmiCounts,
    cell_annotations: CellAnnotations,
    features: Features,
}

type Barcodes = Array1<FixedAscii<64>>;

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
    BarcodeAnnotationMismatch,
    #[error("invalid output path: {reason}")]
    InvalidOutputPath { path: String, reason: String },
    #[error("something went wrong writing H5 object: {reason}")]
    WriteH5Object { path: String, reason: String },
    #[error("something went wrong: {reason}")]
    Other { reason: String },
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
    use std::path::{Path, PathBuf};

    use anyhow::Context;

    use crate::reference_dataset::{
        umi_counts::read_umi_counts_from_h5ad, validate_reference_dataset,
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
            let reference_dataset = validate_reference_dataset(
                path,
                "barcode",
                "annotation",
                "ensembl_id",
                "gene_name",
            )
            .map_err(|e| e[0].clone())
            .context(format!("failed to validate {filename}"))
            .unwrap();

            assert_eq!(
                reference_dataset.counts.data()[0],
                10,
                "first entry in UMI counts of {filename} != 10"
            );

            fake_counts.push(reference_dataset.counts);
        }

        assert_eq!(fake_counts[0], fake_counts[1]);
        assert_eq!(fake_counts[0], fake_counts[2]);
    }

    #[test]
    fn read_real_h5ad() {
        let filename = "test-data/10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx_count_sample_filtered_feature_bc_matrix.h5ad";

        let data =
            validate_reference_dataset(filename, "barcode", "annotation", "gene_ids", "gene_name")
                .map_err(|e| e[0].clone())
                .context(format!("failed to validate {filename}"))
                .unwrap();

        assert!(data.features.id().len() > 30_000)
    }
}
