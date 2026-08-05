use crate::reference_dataset::umi_counts::matrix::Matrix;

#[derive(Clone, Debug, PartialEq)]
pub struct CscMatrix<N>(Matrix<N>);

impl<N> CscMatrix<N>
where
    N: Clone + Default,
{
    pub fn new(mtx: Matrix<N>) -> Self {
        Self(mtx.into_csc())
    }

    pub fn get(&self) -> &Matrix<N> {
        &self.0
    }

    pub fn into_inner(self) -> Matrix<N> {
        self.0
    }
}
