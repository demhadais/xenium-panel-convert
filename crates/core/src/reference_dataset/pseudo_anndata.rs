use serde::Serialize;

use crate::reference_dataset::{
    Barcodes, CellAnnotations, umi_counts::RawCscUmiCounts, var::Features,
};

#[derive(Debug)]
pub struct PseudoAnndata {
    counts: RawCscUmiCounts,
    barcodes: Barcodes,
    cell_annotations: CellAnnotations,
    features: Features,
}

impl PseudoAnndata {
    pub(super) fn new(
        counts: RawCscUmiCounts,
        barcodes: Barcodes,
        cell_annotations: CellAnnotations,
        features: Features,
    ) -> Result<Self, ShapeMismatchError> {
        let n_barcodes = barcodes.len();
        let n_annotations = cell_annotations.len();
        let n_features = features.len();
        let counts_shape = counts.shape();

        let err = Err(ShapeMismatchError {
            n_barcodes,
            n_annotations,
            n_features,
            counts_shape,
        });

        let [n_genes, n_cells] = counts_shape.map(i128::from);

        if n_genes != n_features as i128 {
            return err;
        }

        if n_cells != n_barcodes as i128 || n_cells != n_annotations as i128 {
            return err;
        }

        Ok(Self {
            counts,
            barcodes,
            cell_annotations,
            features,
        })
    }

    pub(super) fn counts(&self) -> &RawCscUmiCounts {
        &self.counts
    }

    pub(super) fn barcodes(&self) -> &Barcodes {
        &self.barcodes
    }

    pub(super) fn cell_annotations(&self) -> &CellAnnotations {
        &self.cell_annotations
    }

    pub(super) fn features(&self) -> &Features {
        &self.features
    }
}

#[derive(Clone, Copy, Debug, Serialize, thiserror::Error)]
#[error(
    "invalid shape - {n_barcodes} barcodes, {n_annotations} cell annotations, {n_features} \
     features, counts shape {counts_shape:?}"
)]
pub struct ShapeMismatchError {
    pub n_barcodes: usize,
    pub n_annotations: usize,
    pub n_features: usize,
    pub counts_shape: [i32; 2],
}
