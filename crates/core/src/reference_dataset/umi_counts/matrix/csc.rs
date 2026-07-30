use crate::reference_dataset::umi_counts::matrix::Matrix;

#[derive(Clone, Debug, PartialEq)]
pub struct CscMatrix(Matrix);

impl CscMatrix {
    pub fn new(mtx: Matrix) -> Self {
        Self(mtx.into_csc())
    }

    pub fn get(&self) -> &Matrix {
        &self.0
    }
}
