use ndarray::{Array2, ArrayView2};
use sprs::CsMatI;

use super::f32_to_i32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawCscUmiCounts(CscMatrix);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CscMatrix(Matrix);

impl CsrMatrix {
    pub fn new(shape: (usize, usize), indptr: Vec<u32>, indices: Vec<u32>, data: Vec<i32>) -> Self {
        Self(Matrix::new(shape, indptr, indices, data))
    }

    pub fn into_raw_umi_counts(self) -> Result<RawCscUmiCounts, super::Error> {
        self.validate_cell_sums_differ()?;

        Ok(RawCscUmiCounts(CscMatrix(self.0.into_csc())))
    }

    fn validate_cell_sums_differ(&self) -> Result<(), super::Error> {
        let mtx = &self.0;

        let mut total_counts_by_cell = mtx
            .outer_iterator()
            .map(|cell| cell.data().iter().sum::<i32>());

        let Some(first) = total_counts_by_cell.next() else {
            return Err(super::Error::EmptyCounts);
        };

        // Write this check using Iterator::any so we can short-circuit on success
        if total_counts_by_cell.any(|s| s != first) {
            Ok(())
        } else {
            Err(super::Error::NormalizedCounts)
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CsrMatrix(Matrix);

impl CscMatrix {
    pub fn new(shape: (usize, usize), indptr: Vec<u32>, indices: Vec<u32>, data: Vec<i32>) -> Self {
        Self(Matrix::new_csc(shape, indptr, indices, data))
    }

    pub fn from_dense(counts: &Array2<f32>) -> Result<Self, super::Error> {
        let data: Vec<_> = counts
            .iter()
            .copied()
            .map(f32_to_i32)
            .collect::<Result<_, _>>()?;

        let counts = ArrayView2::from_shape(counts.dim(), &data)
            .expect("it's the same array, so it has the same shape");

        Ok(Self(Matrix::csc_from_dense(counts, 0)))
    }

    pub fn into_raw_umi_counts(self) -> Result<RawCscUmiCounts, super::Error> {
        // This looks shitty because we are converting to a CSR, which then converts
        // back to a CSC (so that we can reuse code we've already written). These
        // conversions are actually kind of expensive, but testing out the utility shows
        // no perceptible difference, so we opt for the simpler implementation
        let mtx = self.0;
        CsrMatrix(mtx.into_csr()).into_raw_umi_counts()
    }
}

type Matrix = CsMatI<i32, u32>;

#[cfg(test)]
mod tests {
    use ndarray::{Array2, array};

    use crate::reference_dataset::umi_counts::{
        Error, RawCscUmiCounts,
        matrix::{CscMatrix, CsrMatrix, Matrix},
    };

    impl RawCscUmiCounts {
        pub fn data(&self) -> &[i32] {
            self.0.0.data()
        }
    }

    fn counts() -> Array2<i32> {
        array![[0, 1, 2], [2, 4, 6]]
    }

    fn csr() -> CsrMatrix {
        CsrMatrix(Matrix::csr_from_dense(counts().view(), 0))
    }

    fn csc() -> CscMatrix {
        CscMatrix(Matrix::csc_from_dense(counts().view(), 0))
    }

    #[test]
    fn storage_orders_are_equivalent() {
        let from_csr = csr().into_raw_umi_counts().unwrap();
        let from_csc = csc().into_raw_umi_counts().unwrap();

        assert_eq!(
            from_csr, from_csc,
            "the same counts in CSR and CSC did not produce the same matrix"
        );
    }

    #[test]
    fn equal_cell_sums_return_normalized_error() {
        // Both cells sum to 3
        let normalized_counts = array![[0, 1, 2], [1, 1, 1]];
        let csr = CsrMatrix(Matrix::csr_from_dense(normalized_counts.view(), 0));
        let csc = CscMatrix(Matrix::csc_from_dense(normalized_counts.view(), 0));

        std::assert_matches!(csr.into_raw_umi_counts(), Err(Error::NormalizedCounts));
        std::assert_matches!(csc.into_raw_umi_counts(), Err(Error::NormalizedCounts));
    }

    #[test]
    fn dense_mtx_equals_compressed() {
        let float_counts = counts().mapv_into_any(|i| i as f32);
        let from_dense = CscMatrix::from_dense(&float_counts)
            .unwrap()
            .into_raw_umi_counts()
            .unwrap();

        assert_eq!(from_dense, csr().into_raw_umi_counts().unwrap());
        assert_eq!(from_dense, csc().into_raw_umi_counts().unwrap());
    }

    #[test]
    fn dense_transformed_counts_are_rejected() {
        let counts = array![[1., 2.5, 3.], [4., 5., 6.]];

        std::assert_matches!(
            CscMatrix::from_dense(&counts),
            Err(Error::TransformedCounts)
        );
    }
}
