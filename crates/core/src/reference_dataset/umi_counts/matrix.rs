use ndarray::Array2;
use sprs::CsMatI;

use crate::reference_dataset::{
    Error,
    umi_counts::{encoding_type::SparseEncodingType, matrix::csc::CscMatrix},
};

mod csc;

type Matrix = CsMatI<f32, i32>;

#[derive(Clone, Debug, PartialEq)]
pub struct RawCscUmiCounts(CscMatrix);

impl RawCscUmiCounts {
    pub fn from_sparse_matrix(
        shape: (usize, usize),
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
        encoding_type: SparseEncodingType,
    ) -> Result<Self, Error> {
        let mtx = match encoding_type {
            SparseEncodingType::CsrMatrix => Matrix::try_new(shape, indptr, indices, data),
            SparseEncodingType::CscMatrix => Matrix::try_new_csc(shape, indptr, indices, data),
        };

        mtx.map(CscMatrix::new)
            .map_err(|(.., err)| err.into())
            .and_then(Self::from_csc_matrix)
    }

    fn from_csc_matrix(mtx: CscMatrix) -> Result<Self, Error> {
        validate_counts(&mtx)?;
        Ok(Self(mtx))
    }

    pub fn from_dense_matrix(counts: &Array2<f32>) -> Result<Self, Error> {
        Self::from_csc_matrix(CscMatrix::new(Matrix::csc_from_dense(counts.view(), 0.)))
    }
}

fn validate_counts(mtx: &CscMatrix) -> Result<(), Error> {
    let mtx = mtx.get();

    let mut total_counts = vec![0.; mtx.rows()];

    for (cell_idx, count) in mtx.indices().iter().zip(mtx.data()) {
        if !f32_is_nonnegative_integer(*count) {
            return Err(Error::TransformedCounts);
        }

        total_counts[*cell_idx as usize] += *count;
    }

    if all_total_counts_are_equal(&total_counts)? {
        return Err(Error::NormalizedCounts);
    }

    Ok(())
}

fn all_total_counts_are_equal(total_counts: &[f32]) -> Result<bool, Error> {
    let mut nonempty_cells = total_counts.iter().copied().filter(|tc| *tc > 0.);
    let first_nonempty_cell = nonempty_cells.next().ok_or(Error::EmptyCounts)?;

    Ok(nonempty_cells.all(|tc| tc == first_nonempty_cell))
}

#[allow(clippy::float_cmp)]
fn f32_is_nonnegative_integer(f: f32) -> bool {
    f.round() == f && f >= 0.0
}

#[cfg(test)]
mod tests {
    use ndarray::{Array2, array};

    use crate::reference_dataset::umi_counts::{
        Error, RawCscUmiCounts,
        matrix::{Matrix, csc::CscMatrix},
    };

    impl RawCscUmiCounts {
        pub fn data(&self) -> &[f32] {
            self.0.get().data()
        }
    }

    fn counts() -> Array2<f32> {
        array![[0., 1., 2.], [2., 4., 6.]]
    }

    fn csr() -> Matrix {
        Matrix::csr_from_dense(counts().view(), 0.)
    }

    fn csc() -> CscMatrix {
        CscMatrix::new(Matrix::csc_from_dense(counts().view(), 0.))
    }

    #[test]
    fn storage_orders_are_equivalent() {
        let from_csr = RawCscUmiCounts::from_csc_matrix(CscMatrix::new(csr())).unwrap();
        let from_csc = RawCscUmiCounts::from_csc_matrix(csc()).unwrap();

        assert_eq!(
            from_csr, from_csc,
            "the same counts in CSR and CSC did not produce the same matrix"
        );

        let from_dense = RawCscUmiCounts::from_dense_matrix(&counts()).unwrap();
        assert_eq!(
            from_csr, from_dense,
            "dense matrix and sparse matrix provided from the same counts did not produce the \
             same matrix"
        );
    }

    #[test]
    fn equal_cell_sums_return_normalized_error() {
        // Both cells sum to 3
        let normalized_counts = array![[0., 1., 2.], [1., 1., 1.]];

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&normalized_counts),
            Err(Error::NormalizedCounts)
        );
    }

    #[test]
    fn all_zero_counts_are_empty() {
        let all_zero_counts = Array2::zeros((2, 2));

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&all_zero_counts),
            Err(Error::EmptyCounts)
        );
    }

    #[test]
    fn transformed_counts_error() {
        let counts = array![[1., 2.5, 3.], [4., 5., 6.]];

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&counts),
            Err(Error::TransformedCounts)
        );
    }
}
