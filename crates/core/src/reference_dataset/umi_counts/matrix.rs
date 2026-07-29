use ndarray::{Array2, ArrayView2};
use sprs::{CsMatBase, CsMatI};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawCscUmiCounts(CscMatrix);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CscMatrix(Matrix);

impl CsrMatrix {
    pub fn new(shape: (usize, usize), indptr: Vec<u32>, indices: Vec<u32>, data: Vec<i32>) -> Self {
        Self(Matrix::new(shape, indptr, indices, data))
    }

    pub fn into_raw_umi_counts(self) -> Result<RawCscUmiCounts, super::Error> {
        let mtx = self.0;

        let mut total_counts_per_cell = mtx
            .outer_iterator()
            .map(|cell| cell.data().iter().sum::<i32>());

        let Some(first) = total_counts_per_cell.next() else {
            return Err(super::Error::EmptyCounts);
        };

        if total_counts_per_cell.any(|s| s != first) {
            drop(total_counts_per_cell);
            Ok(RawCscUmiCounts(CscMatrix(mtx.into_csc())))
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

    pub fn into_raw_umi_counts(self) -> Result<RawCscUmiCounts, super::Error> {
        let mtx = self.0.into_csr();
        CsrMatrix(mtx).into_raw_umi_counts()
    }

    pub fn from_dense(m: ArrayView2<i32>) -> Self {
        Self(Matrix::csc_from_dense(m, 0))
    }
}

type Matrix = CsMatI<i32, u32>;

#[cfg(test)]
impl RawCscUmiCounts {
    pub fn data(&self) -> &[i32] {
        self.0.0.data()
    }
}
