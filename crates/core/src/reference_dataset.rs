use std::fs;

use camino::Utf8Path;
use hdf5_metno::{File, types::FixedAscii};
use ndarray::Array1;
use serde::Serialize;

use crate::{
    common::ErrorVecExt,
    reference_dataset::{
        columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneNameCol},
        error::{
            ReadReferenceDatasetError, ReadReferenceDatasetErrors, WriteReferenceDatasetError,
            WriteReferenceDatasetErrorWrapper,
        },
        feature_set::FeatureSet,
        h5_util::{create_h5_group, write_dataset_to_h5_group},
        obs::{read_cell_annotations_from_h5ad, read_cell_barcodes_from_h5ad},
        pseudo_anndata::PseudoAnndata,
        umi_counts::read_umi_counts_from_h5ad,
        var::read_features_from_h5ad,
    },
};

pub mod columns;
pub mod error;
pub mod feature_set;
pub mod h5_util;
pub mod obs;
pub mod pseudo_anndata;
pub mod umi_counts;
pub mod var;

pub fn read_reference_dataset(
    path: &Utf8Path,
    cell_barcode_col: &CellBarcodeCol,
    cell_annotation_col: &CellAnnotationCol,
    ensembl_id_col: &EnsemblIdCol,
    gene_name_col: &GeneNameCol,
    feature_set: FeatureSet,
) -> Result<PseudoAnndata, ReadReferenceDatasetErrors> {
    let mut errors = Vec::new();

    let collect = |errors: Vec<ReadReferenceDatasetError>| ReadReferenceDatasetErrors {
        path: path.to_owned(),
        errors: errors.into_iter().map(Into::into).collect(),
    };

    let file = hdf5_metno::File::open(path).map_err(|e| {
        collect(vec![ReadReferenceDatasetError::InvalidH5File {
            reason: e.to_string(),
        }])
    })?;

    let counts = read_umi_counts_from_h5ad(&file).map_or_else(|err| errors.push_err(err), Some);

    let barcodes = read_cell_barcodes_from_h5ad(&file, cell_barcode_col)
        .map_or_else(|err| errors.push_err(err), Some);

    let cell_annotations = read_cell_annotations_from_h5ad(&file, cell_annotation_col)
        .map_or_else(|err| errors.push_err(err), Some);

    let features = read_features_from_h5ad(&file, ensembl_id_col, gene_name_col, feature_set)
        .map_or_else(|err| errors.push_err(err), Some);

    let (Some(counts), Some(barcodes), Some(cell_annotations), Some(features)) =
        (counts, barcodes, cell_annotations, features)
    else {
        return Err(collect(errors));
    };

    let anndata =
        PseudoAnndata::new(counts, barcodes, cell_annotations, features).map_err(|err| {
            errors.push_err::<()>(err);
            collect(errors)
        })?;

    Ok(anndata)
}

// See https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices for format specifics
pub fn write_reference_dataset(
    dir: &Utf8Path,
    ds: &PseudoAnndata,
) -> Result<(), WriteReferenceDatasetErrorWrapper> {
    if !dir.exists() {
        fs::create_dir_all(dir).map_err(|e| WriteReferenceDatasetError::CreateOutputDir {
            path: dir.to_owned(),
            reason: e.to_string(),
        })?;
    }

    let annotations_path = dir.join("annotations.csv");
    if annotations_path.exists() {
        return Err(WriteReferenceDatasetError::AnnotationsCsvExists {
            path: annotations_path,
        }
        .into());
    }

    let matrix_path = dir.join("matrix.h5");
    let file = File::create_excl(&matrix_path).map_err(|e| {
        WriteReferenceDatasetError::CreateMatrixFile {
            path: matrix_path.clone(),
            reason: e.to_string(),
        }
    })?;

    let matrix_group = create_h5_group(&file, "matrix").map_err(|error| {
        WriteReferenceDatasetError::CreateH5Group {
            path: matrix_path.clone(),
            error,
        }
    })?;

    let write_err = |error| WriteReferenceDatasetError::WriteH5Dataset {
        path: matrix_path.clone(),
        error,
    };

    let barcodes = ds.barcodes();
    write_dataset_to_h5_group(&matrix_group, "barcodes", barcodes).map_err(write_err)?;

    let counts = ds.counts();
    write_dataset_to_h5_group(&matrix_group, "data", counts.data()).map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "indices", counts.indices()).map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "indptr", counts.indptr().iter().as_slice())
        .map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "shape", &counts.shape()).map_err(write_err)?;

    let features = ds.features();
    write_dataset_to_h5_group(
        &matrix_group,
        "features/feature_type",
        features.feature_types(),
    )
    .map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "features/id", features.ensembl_ids())
        .map_err(write_err)?;
    write_dataset_to_h5_group(&matrix_group, "features/name", features.gene_names())
        .map_err(write_err)?;

    write_annotations_csv(&annotations_path, barcodes, ds.cell_annotations());

    Ok(())
}

fn write_annotations_csv(path: &Utf8Path, barcodes: &Barcodes, annotations: &CellAnnotations) {
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
}

type Barcode = FixedAscii<64>;

type Barcodes = Array1<Barcode>;

type CellAnnotations = Array1<String>;

#[cfg(test)]
mod tests {
    use camino::Utf8Path;
    use hdf5_metno::types::FixedAscii;

    use crate::reference_dataset::{
        feature_set::{FeatureSet, Transcriptome},
        h5_util::read_1d_dataset,
        read_reference_dataset,
        var::{EnsemblId, GeneName},
    };

    #[test]
    fn read_generated_h5ad_files() {
        let paths = [
            "test-data/csr_adata.h5ad",
            "test-data/csc_adata.h5ad",
            "test-data/dense_adata.h5ad",
        ]
        .map(Utf8Path::new);

        let mut fake_counts = Vec::with_capacity(paths.len());

        for path in paths {
            let reference_dataset = read_reference_dataset(
                path,
                &"barcode".into(),
                &"annotation".into(),
                &"ensembl_id".into(),
                &"gene_name".into(),
                FeatureSet::new(Transcriptome::Grch382020A, false),
            )
            .unwrap_or_else(|e| {
                panic!(
                    "failed to read {path}: {}",
                    serde_json::to_string(&e).unwrap()
                )
            });

            assert_eq!(
                reference_dataset.counts().data()[0],
                10,
                "first entry in UMI counts of {path} != 10"
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
            Utf8Path::new(filename),
            &"barcode".into(),
            &"annotation".into(),
            &"gene_ids".into(),
            &"gene_name".into(),
            FeatureSet::new(Transcriptome::Grch382020A, false),
        )
        .unwrap_or_else(|e| {
            panic!(
                "failed to validate {filename}: {}",
                serde_json::to_string(&e).unwrap()
            )
        });

        // Compare it against the original data
        let original_h5 = hdf5_metno::File::open(
            "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5",
        )
        .unwrap();

        let original_counts = read_1d_dataset::<i32>(&original_h5, "matrix/data").unwrap();
        assert_eq!(
            original_counts.as_slice().unwrap(),
            reconstructed_dataset.counts().data(),
            "UMI counts were not correctly? reconstructed"
        );

        let original_indices = read_1d_dataset::<i64>(&original_h5, "matrix/indices").unwrap();
        assert_eq!(
            original_indices.as_slice().unwrap(),
            reconstructed_dataset.counts().indices(),
            "UMI count indices were not correctly? reconstructed"
        );

        let original_indptr = read_1d_dataset::<i64>(&original_h5, "matrix/indptr").unwrap();
        assert_eq!(
            original_indptr.as_slice().unwrap(),
            reconstructed_dataset.counts().indptr().iter().as_slice(),
            "UMI counts indptr was not correctly? reconstructed"
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
