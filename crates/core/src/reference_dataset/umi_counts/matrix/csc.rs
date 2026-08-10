use crate::reference_dataset::umi_counts::matrix::Matrix;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct CscMatrix<N>(Matrix<N>);

impl<N> CscMatrix<N>
where
    N: Clone + Default,
{
    pub(super) fn new(mtx: Matrix<N>) -> Self {
        Self(mtx.into_csc())
    }

    pub(super) fn as_matrix(&self) -> &Matrix<N> {
        &self.0
    }

    pub(super) fn into_inner(self) -> Matrix<N> {
        self.0
    }
}
