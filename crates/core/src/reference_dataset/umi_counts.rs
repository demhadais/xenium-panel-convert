use hdf5_metno::{Dataset, File, Group};
pub use matrix::RawCscUmiCounts;
use ndarray::{Array2, ArrayView2, Axis};
use serde::Serialize;

use crate::reference_dataset::umi_counts::{
    encoding_type::{DenseEncodingType, SparseEncodingType},
    matrix::{CscMatrix, CsrMatrix},
};

mod encoding_type;
mod matrix;

pub fn read_umi_counts_from_h5ad(file: &File) -> Result<RawCscUmiCounts, Error> {
    const X: &str = "X";
    // The actual counts (stored by scanpy as the highly-descriptive name "X") are usually stored as in compressed-sparse matrix format in a group, but they might also be in a dataset
    match file.group(X) {
        Ok(g) => read_x_group(&g),
        Err(_) => read_x_dataset(&file.dataset(X)?),
    }
}

fn read_x_group(x: &Group) -> Result<RawCscUmiCounts, Error> {
    let encoding_type = SparseEncodingType::from_group(x)?;

    let data = x.dataset("data").and_then(|ds| ds.read_raw())?;
    let data = data.into_iter().map(f32_to_i32).collect::<Result<_, _>>()?;

    let indptr: Vec<usize> = x.dataset("indptr").and_then(|ds| ds.read_raw())?;
    let indices: Vec<usize> = x.dataset("indices").and_then(|ds| ds.read_raw())?;

    // It's very nice that scanpy decides to store the shape as an attribute rather than the following 10x Genomics and storing it as a dataset. It's great when a library built to analyze data ends up changing the format of the data :)
    let shape = x.attr("shape").and_then(|sh| sh.read_1d())?;
    let shape = (shape[0], shape[1]);

    let counts = match encoding_type {
        SparseEncodingType::CsrMatrix => {
            CsrMatrix::new(shape, indptr, indices, data).into_raw_umi_counts()
        }
        SparseEncodingType::CscMatrix => {
            CscMatrix::new(shape, indptr, indices, data).into_raw_umi_counts()
        }
    };

    counts
}

fn read_x_dataset(x: &Dataset) -> Result<RawCscUmiCounts, Error> {
    let encoding_type = DenseEncodingType::from_dataset(x)?;
    let data = match encoding_type {
        DenseEncodingType::Array => x.read_2d()?,
    };

    let shape = data.raw_dim();

    let data: Vec<_> = data.into_iter().map(f32_to_i32).collect::<Result<_, _>>()?;

    let mtx = CscMatrix::from_dense(
        ArrayView2::from_shape(shape, &data)
            .expect("it's the same array, so it has the same shape"),
    );

    mtx.into_raw_umi_counts()
}

fn f32_to_i32(f: f32) -> Result<i32, Error> {
    let is_integer = f.round() == f;
    let is_nonnegative = f >= 0.0;

    if is_integer && is_nonnegative {
        Ok(f as i32)
    } else {
        Err(Error::TransformedCounts)
    }
}

#[derive(Debug, Serialize, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum Error {
    #[error("empty UMI counts")]
    EmptyCounts,
    #[error("transformed UMI counts")]
    TransformedCounts,
    #[error("normalized UMI counts")]
    NormalizedCounts,
    #[error("HDF5 error: {reason}")]
    Hdf5 { reason: String },
    #[error("unknown encoding type")]
    UnknownEncodingType,
}

impl From<hdf5_metno::Error> for Error {
    fn from(err: hdf5_metno::Error) -> Self {
        Self::Hdf5 {
            reason: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use hdf5_metno::File;

    use crate::reference_dataset::umi_counts::read_umi_counts_from_h5ad;

    #[test]
    fn read_h5ad_files() {
        // The last one apparently has some compression applied to it because of scanpy.write's default behavior. It's really nice that scanpy's default behavior differs from anndata's default behavior :)
        let files = ["csr_adata", "csc_adata", "dense_adata", "WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix"]
            .map(|fname| format!("test-data/{fname}.h5ad"))
            .map(|path| File::open(path).unwrap());

        let mut all_counts = Vec::with_capacity(3);

        for f in files {
            let filename = f.filename();
            let counts = read_umi_counts_from_h5ad(&f)
                .expect(&format!("failed to read UMI counts from {filename}"));

            if filename.contains("adata") {
                assert_eq!(
                    counts.data()[0],
                    10,
                    "first entry in UMI counts of {filename} != 10"
                );
            }

            all_counts.push(counts);
        }

        // We know the first 3 files are generated from the same data
        assert_eq!(all_counts[0], all_counts[1]);
        assert_eq!(all_counts[0], all_counts[2]);
    }
}
