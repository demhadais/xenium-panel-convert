use std::borrow::Cow;

use ndarray::Array2;
use sprs::CsMatI;

use crate::reference_dataset::umi_counts::{
    Error,
    encoding_type::{DenseEncodingType, SparseEncodingType},
    matrix::csc::CscMatrix,
};

mod csc;

type Matrix<N> = CsMatI<N, u32>;

#[derive(Clone, Debug, PartialEq)]
pub struct RawCscUmiCounts(CscMatrix<u32>);

impl RawCscUmiCounts {
    pub fn from_sparse_matrix(
        shape: (usize, usize),
        indptr: Vec<u32>,
        indices: Vec<u32>,
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

    fn from_csc_matrix(mtx: CscMatrix<f32>) -> Result<Self, Error> {
        let mtx = convert_to_u32_mtx(mtx)?;
        Ok(Self(mtx))
    }

    pub fn from_dense_matrix(
        counts: &Array2<f32>,
        encoding_type: DenseEncodingType,
    ) -> Result<Self, Error> {
        match encoding_type {
            DenseEncodingType::Array => {
                Self::from_csc_matrix(CscMatrix::new(Matrix::csc_from_dense(counts.view(), 0.)))
            }
        }
    }

    pub fn data(&self) -> &[u32] {
        self.0.get().data()
    }

    pub fn indices(&self) -> &[u32] {
        self.0.get().indices()
    }

    pub fn indptr(&self) -> Cow<'_, [u32]> {
        self.0.get().proper_indptr()
    }

    pub fn shape(&self) -> [u64; 2] {
        let (nrows, ncols) = self.0.get().shape();

        [nrows as u64, ncols as u64]
    }
}

fn convert_to_u32_mtx(mtx: CscMatrix<f32>) -> Result<CscMatrix<u32>, Error> {
    let shape = mtx.get().shape();
    let nrows = mtx.get().rows();

    let (indptr, indices, f32_data) = mtx.into_inner().into_raw_storage();

    let mut total_counts = vec![0; nrows];
    let mut i32_data: Vec<_> = (0..f32_data.len()).map(|_| 0).collect();

    for (cell_idx, count) in indices.iter().zip(&f32_data) {
        let cell_idx = *cell_idx as usize;
        let count = f32_to_u32(*count)?;

        i32_data[cell_idx] = count;
        total_counts[cell_idx] += count;
    }

    if all_total_counts_are_equal(&total_counts)? {
        return Err(Error::NormalizedCounts);
    }

    Ok(CscMatrix::new(Matrix::new_csc(
        shape, indptr, indices, i32_data,
    )))
}

fn all_total_counts_are_equal(total_counts: &[u32]) -> Result<bool, Error> {
    let mut nonempty_cells = total_counts.iter().copied().filter(|tc| *tc > 0);
    let first_nonempty_cell = nonempty_cells.next().ok_or(Error::EmptyCounts)?;

    Ok(nonempty_cells.all(|tc| tc == first_nonempty_cell))
}

#[allow(clippy::float_cmp)]
fn f32_to_u32(f: f32) -> Result<u32, Error> {
    let is_nonnegative_integral = f.round() == f && f >= 0.0;

    if is_nonnegative_integral {
        Ok(f as u32)
    } else {
        Err(Error::TransformedCounts)
    }
}

#[cfg(test)]
mod tests {
    use ndarray::{Array2, array};

    use crate::reference_dataset::umi_counts::{
        Error, RawCscUmiCounts,
        encoding_type::DenseEncodingType,
        matrix::{Matrix, csc::CscMatrix},
    };

    fn counts() -> Array2<f32> {
        array![[0., 1., 2.], [2., 4., 6.]]
    }

    fn csr() -> Matrix<f32> {
        Matrix::csr_from_dense(counts().view(), 0.)
    }

    fn csc() -> CscMatrix<f32> {
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

        let from_dense =
            RawCscUmiCounts::from_dense_matrix(&counts(), DenseEncodingType::Array).unwrap();
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
            RawCscUmiCounts::from_dense_matrix(&normalized_counts, DenseEncodingType::Array),
            Err(Error::NormalizedCounts)
        );
    }

    #[test]
    fn all_zero_counts_are_empty() {
        let all_zero_counts = Array2::zeros((2, 2));

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&all_zero_counts, DenseEncodingType::Array),
            Err(Error::EmptyCounts)
        );
    }

    #[test]
    fn transformed_counts_error() {
        let counts = array![[1., 2.5, 3.], [4., 5., 6.]];

        std::assert_matches!(
            RawCscUmiCounts::from_dense_matrix(&counts, DenseEncodingType::Array),
            Err(Error::TransformedCounts)
        );
    }
}
