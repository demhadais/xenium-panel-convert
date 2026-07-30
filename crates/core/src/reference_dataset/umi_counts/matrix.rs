use ndarray::Array2;
use sprs::CsMatI;

use crate::reference_dataset::umi_counts::matrix::cell_totals::CellTotals;

#[derive(Clone, Debug, PartialEq)]
pub struct RawCscUmiCounts(CscMatrix);

#[derive(Clone, Debug, PartialEq)]
pub struct CscMatrix(Matrix);

#[derive(Clone, Debug, PartialEq)]
pub struct CsrMatrix(Matrix);

impl CsrMatrix {
    pub fn new(
        shape: (usize, usize),
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
    ) -> Result<Self, super::Error> {
        Matrix::try_new(shape, indptr, indices, data)
            .map(Self)
            .map_err(|(.., structure_error)| structure_error.into())
    }

    pub fn into_raw_umi_counts(self) -> Result<RawCscUmiCounts, super::Error> {
        self.validate_counts()?;

        Ok(RawCscUmiCounts(CscMatrix(self.0.into_csc())))
    }

    fn validate_counts(&self) -> Result<(), super::Error> {
        let mtx = &self.0;

        if mtx.nnz() == 0 {
            return Err(super::Error::EmptyCounts);
        }

        let mut cell_totals = CellTotals::default();

        for cell in mtx.outer_iterator() {
            let counts = cell.data();

            if !counts.iter().copied().all(f32_is_nonnegative_integer) {
                return Err(super::Error::TransformedCounts);
            }

            cell_totals.add(counts.iter().copied().sum());
        }

        cell_totals.check_differ()
    }
}

impl CscMatrix {
    pub fn new(
        shape: (usize, usize),
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
    ) -> Result<Self, super::Error> {
        Matrix::try_new_csc(shape, indptr, indices, data)
            .map(Self)
            .map_err(|(.., structure_error)| structure_error.into())
    }

    pub fn from_dense(counts: &Array2<f32>) -> Self {
        Self(Matrix::csc_from_dense(counts.view(), 0.))
    }

    pub fn into_raw_umi_counts(self) -> Result<RawCscUmiCounts, super::Error> {
        self.validate_counts()?;

        Ok(RawCscUmiCounts(self))
    }

    fn validate_counts(&self) -> Result<(), super::Error> {
        let mtx = &self.0;

        if mtx.nnz() == 0 {
            return Err(super::Error::EmptyCounts);
        }

        // Cells are the inner dimension of a CSC matrix, so no cell's total is known
        // until the whole matrix has been traversed
        let mut totals_by_cell = vec![0.; mtx.rows()];
        let cell_indices = mtx.indices();

        for (cell_idx, count) in cell_indices.iter().zip(mtx.data()) {
            if !f32_is_nonnegative_integer(*count) {
                return Err(super::Error::TransformedCounts);
            }

            totals_by_cell[*cell_idx as usize] += *count;
        }

        let mut cell_totals = CellTotals::default();

        for total in totals_by_cell {
            cell_totals.add(total);
        }

        cell_totals.check_differ()
    }
}

#[allow(clippy::float_cmp)]
fn f32_is_nonnegative_integer(f: f32) -> bool {
    f.round() == f && f >= 0.0
}

type Matrix = CsMatI<f32, i32>;

mod cell_totals {
    use crate::reference_dataset::umi_counts::Error;

    #[derive(Default)]
    pub struct CellTotals {
        first: Option<f32>,
        n_cells: usize,
        any_differ: bool,
    }

    impl CellTotals {
        #[allow(clippy::float_cmp)]
        pub fn add(&mut self, cell_total: f32) {
            self.n_cells += 1;

            match self.first {
                None => self.first = Some(cell_total),
                Some(first) => self.any_differ |= cell_total != first,
            }
        }

        pub fn check_differ(&self) -> Result<(), Error> {
            if self.n_cells == 1 || self.any_differ {
                Ok(())
            } else {
                Err(Error::NormalizedCounts)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use ndarray::{Array2, array};

    use crate::reference_dataset::umi_counts::{
        Error, RawCscUmiCounts,
        matrix::{CscMatrix, CsrMatrix, Matrix},
    };

    impl RawCscUmiCounts {
        pub fn data(&self) -> &[f32] {
            self.0.0.data()
        }
    }

    fn counts() -> Array2<f32> {
        array![[0., 1., 2.], [2., 4., 6.]]
    }

    fn csr() -> CsrMatrix {
        CsrMatrix(Matrix::csr_from_dense(counts().view(), 0.))
    }

    fn csc() -> CscMatrix {
        CscMatrix(Matrix::csc_from_dense(counts().view(), 0.))
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
        let normalized_counts = array![[0., 1., 2.], [1., 1., 1.]];
        let csr = CsrMatrix(Matrix::csr_from_dense(normalized_counts.view(), 0.));
        let csc = CscMatrix(Matrix::csc_from_dense(normalized_counts.view(), 0.));

        std::assert_matches!(csr.into_raw_umi_counts(), Err(Error::NormalizedCounts));
        std::assert_matches!(csc.into_raw_umi_counts(), Err(Error::NormalizedCounts));
    }

    #[test]
    fn single_cell_counts_are_not_normalized() {
        let single_cell_counts = array![[1., 2., 3.]];
        let csr = CsrMatrix(Matrix::csr_from_dense(single_cell_counts.view(), 0.));
        let csc = CscMatrix(Matrix::csc_from_dense(single_cell_counts.view(), 0.));

        csr.into_raw_umi_counts().unwrap();
        csc.into_raw_umi_counts().unwrap();
    }

    #[test]
    fn all_zero_counts_are_empty() {
        let all_zero_counts = Array2::zeros((2, 3));
        let csr = CsrMatrix(Matrix::csr_from_dense(all_zero_counts.view(), 0.));
        let csc = CscMatrix(Matrix::csc_from_dense(all_zero_counts.view(), 0.));

        std::assert_matches!(csr.into_raw_umi_counts(), Err(Error::EmptyCounts));
        std::assert_matches!(csc.into_raw_umi_counts(), Err(Error::EmptyCounts));
    }

    #[test]
    fn dense_mtx_equals_compressed() {
        let from_dense = CscMatrix::from_dense(&counts())
            .into_raw_umi_counts()
            .unwrap();

        assert_eq!(from_dense, csr().into_raw_umi_counts().unwrap());
        assert_eq!(from_dense, csc().into_raw_umi_counts().unwrap());
    }

    #[test]
    fn dense_transformed_counts_are_rejected() {
        let counts = array![[1., 2.5, 3.], [4., 5., 6.]];

        std::assert_matches!(
            CscMatrix::from_dense(&counts).into_raw_umi_counts(),
            Err(Error::TransformedCounts)
        );
    }
}
