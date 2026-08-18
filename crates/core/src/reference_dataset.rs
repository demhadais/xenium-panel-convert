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
    use std::fs;

    use camino::Utf8Path;
    use hdf5_metno::{
        File,
        types::{FixedAscii, TypeDescriptor},
    };
    use tempfile::TempDir;

    use crate::reference_dataset::{
        Barcode,
        error::{
            ReadReferenceDatasetError, ReadReferenceDatasetErrorWrapper, WriteReferenceDatasetError,
        },
        feature_set::{FeatureSet, Transcriptome},
        h5_util::read_1d_dataset,
        pseudo_anndata::PseudoAnndata,
        read_reference_dataset,
        var::{EnsemblId, GeneName},
        write_reference_dataset,
    };

    const REAL_H5AD: &str = "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5ad";
    const REAL_H5: &str = "test-data/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5";

    fn real_dataset() -> PseudoAnndata {
        read_reference_dataset(
            Utf8Path::new(REAL_H5AD),
            &"barcode".into(),
            &"annotation".into(),
            &"gene_ids".into(),
            &"gene_name".into(),
            FeatureSet::new(Transcriptome::Mm102020A, false),
        )
        .unwrap()
    }

    fn temp_dir() -> TempDir {
        TempDir::new().unwrap()
    }

    fn utf8_path_from_temp_dir(dir: &TempDir) -> &Utf8Path {
        Utf8Path::from_path(dir.path()).unwrap()
    }

    fn dtype(file: &File, path: &str) -> TypeDescriptor {
        file.dataset(path)
            .unwrap()
            .dtype()
            .unwrap()
            .to_descriptor()
            .unwrap()
    }

    #[test]
    fn read_collects_errors_from_every_field() {
        let path = "test-data/csr_adata.h5ad";

        let errors = read_reference_dataset(
            Utf8Path::new(path),
            &"nonexistent".into(),
            &"nonexistent".into(),
            &"nonexistent".into(),
            &"nonexistent".into(),
            FeatureSet::new(Transcriptome::Grch382020A, false),
        )
        .unwrap_err();

        assert_eq!(errors.path.as_str(), path);

        std::assert_matches!(
            errors.errors.as_slice(),
            [
                ReadReferenceDatasetErrorWrapper {
                    error: ReadReferenceDatasetError::Obs { .. },
                    ..
                },
                ReadReferenceDatasetErrorWrapper {
                    error: ReadReferenceDatasetError::Obs { .. },
                    ..
                },
                ReadReferenceDatasetErrorWrapper {
                    error: ReadReferenceDatasetError::Var { .. },
                    ..
                },
            ],
            "every unreadable field should be reported, not just the first"
        );
    }

    #[test]
    fn write_dataset_produces_cellranger_h5() {
        let dataset = real_dataset();
        let dir = temp_dir();
        let output_dir = utf8_path_from_temp_dir(&dir).join("dataset");

        write_reference_dataset(&output_dir, &dataset).unwrap();

        let written = File::open(output_dir.join("matrix.h5")).unwrap();

        // Every array should land at the path cellranger puts it
        let counts = dataset.counts();
        let shape = counts.shape();

        assert_eq!(
            read_1d_dataset::<i32>(&written, "matrix/data")
                .unwrap()
                .as_slice()
                .unwrap(),
            counts.data()
        );
        assert_eq!(
            read_1d_dataset::<i64>(&written, "matrix/indices")
                .unwrap()
                .as_slice()
                .unwrap(),
            counts.indices()
        );
        assert_eq!(
            read_1d_dataset::<i64>(&written, "matrix/indptr")
                .unwrap()
                .as_slice()
                .unwrap(),
            counts.indptr().iter().as_slice()
        );
        assert_eq!(
            read_1d_dataset::<i32>(&written, "matrix/shape")
                .unwrap()
                .as_slice()
                .unwrap(),
            shape.as_slice()
        );

        assert_eq!(
            read_1d_dataset::<Barcode>(&written, "matrix/barcodes").unwrap(),
            dataset.barcodes()
        );

        let features = dataset.features();
        assert_eq!(
            read_1d_dataset::<EnsemblId>(&written, "matrix/features/id").unwrap(),
            features.ensembl_ids()
        );
        assert_eq!(
            read_1d_dataset::<GeneName>(&written, "matrix/features/name").unwrap(),
            features.gene_names()
        );
        assert_eq!(
            read_1d_dataset::<FixedAscii<32>>(&written, "matrix/features/feature_type").unwrap(),
            features.feature_types()
        );

        // The panel designer reads the counts as integers, so they must have the
        // same on-disk types cellranger writes. The string types are allowed to
        // differ, since we write wider fixed-length strings than cellranger does
        let original = File::open(REAL_H5).unwrap();
        for path in [
            "matrix/data",
            "matrix/indices",
            "matrix/indptr",
            "matrix/shape",
        ] {
            assert_eq!(
                dtype(&written, path),
                dtype(&original, path),
                "{path} was not written with the type cellranger writes"
            );
        }

        let annotations = fs::read_to_string(output_dir.join("annotations.csv")).unwrap();
        let mut rows = annotations.lines();

        assert_eq!(rows.next(), Some("barcode,annotation"));

        let rows: Vec<_> = rows.collect();
        assert_eq!(
            rows.len(),
            dataset.barcodes().len(),
            "expected one row per barcode"
        );
        // The annotations of the test-data alternate between two groups
        assert!(rows[0].ends_with(",group_0"));
        assert!(rows[1].ends_with(",group_1"));
    }

    #[test]
    fn does_not_overwrite_existing_outputs() {
        let dataset = real_dataset();
        let dir = temp_dir();

        let existing_annotations = utf8_path_from_temp_dir(&dir).join("existing-annotations");
        fs::create_dir_all(&existing_annotations).unwrap();
        fs::write(existing_annotations.join("annotations.csv"), "").unwrap();

        std::assert_matches!(
            write_reference_dataset(&existing_annotations, &dataset).map_err(|e| e.error),
            Err(WriteReferenceDatasetError::AnnotationsCsvExists { .. })
        );

        let existing_matrix = utf8_path_from_temp_dir(&dir).join("existing-matrix");
        fs::create_dir_all(&existing_matrix).unwrap();
        fs::write(existing_matrix.join("matrix.h5"), "").unwrap();

        std::assert_matches!(
            write_reference_dataset(&existing_matrix, &dataset).map_err(|e| e.error),
            Err(WriteReferenceDatasetError::CreateMatrixFile { .. })
        );
    }

    #[test]
    fn read_real_h5ad() {
        let reconstructed_dataset = real_dataset();

        // Compare it against the original data
        let original_h5 = hdf5_metno::File::open(REAL_H5).unwrap();

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
