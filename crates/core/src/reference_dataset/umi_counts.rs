use hdf5_metno::{Dataset, File, Group};
use ndarray::{Array2, ArrayView2, Axis};
use serde::Serialize;
use sprs::CsMatBase;

use crate::reference_dataset::umi_counts::encoding_type::{DenseEncodingType, SparseEncodingType};

mod encoding_type;

type UmiCounts = CsMatBase<i32, usize, Vec<usize>, Vec<usize>, Vec<i32>>;

pub struct RawUmiCounts(UmiCounts);

pub fn read_umi_counts_from_h5ad(file: &File) -> Result<RawUmiCounts, Error> {
    const X: &str = "X";
    // The actual counts (stored by scanpy as the highly-descriptive name "X") are usually stored as in compressed-sparse matrix format in a group, but they might also be in a dataset
    match file.group(X) {
        Ok(g) => read_x_group(&g),
        Err(_) => read_x_dataset(&file.dataset(X)?),
    }
}

fn read_x_group(x: &Group) -> Result<RawUmiCounts, Error> {
    let encoding_type = SparseEncodingType::from_group(x)?;

    let data = x.dataset("data").and_then(|ds| ds.read_raw())?;
    let data = data.into_iter().map(f32_to_i32).collect::<Result<_, _>>()?;

    let indptr: Vec<usize> = x.dataset("indptr").and_then(|ds| ds.read_raw())?;
    let indices: Vec<usize> = x.dataset("indices").and_then(|ds| ds.read_raw())?;

    // It's very nice that scanpy decides to store the shape as an attribute rather than the following 10x Genomics and storing it as a dataset. It's great when a library built to analyze data ends up changing the format of the data :)
    let shape = x.attr("shape").and_then(|sh| sh.read_1d())?;
    let shape = (shape[0], shape[1]);

    let counts = match encoding_type {
        SparseEncodingType::CsrMatrix => UmiCounts::new(shape, indptr, indices, data),
        SparseEncodingType::CscMatrix => UmiCounts::new_csc(shape, indptr, indices, data),
    };

    validate_total_counts_are_different_sparse(&counts)?;

    Ok(RawUmiCounts(counts))
}

fn read_x_dataset(x: &Dataset) -> Result<RawUmiCounts, Error> {
    let encoding_type = DenseEncodingType::from_dataset(x)?;
    let data = match encoding_type {
        DenseEncodingType::Array => x.read_2d()?,
    };

    validate_total_counts_are_different_dense(&data)?;

    let shape = data.raw_dim();

    let data: Vec<_> = data.into_iter().map(f32_to_i32).collect::<Result<_, _>>()?;

    Ok(ArrayView2::from_shape(shape, &data)
        .map(|a| UmiCounts::csc_from_dense(a, 0))
        .map(RawUmiCounts)
        .expect("it's the same array, so it has the same shape"))
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

fn validate_total_counts_are_different_dense(counts: &Array2<f32>) -> Result<(), Error> {
    let total_counts_per_cell = counts.sum_axis(Axis(1));
    let Some(first_nonzero) = total_counts_per_cell.first().copied() else {
        return Err(Error::EmptyCounts);
    };

    let is_different_from_first = |f: f32| {
        let diff: f32 = first_nonzero - f;
        diff.abs() > 1.
    };

    let mut cell_sums = total_counts_per_cell.into_iter();
    // Advance the iterator because obviously the first element will cause this function to fail
    cell_sums.next();

    // Write this check using Iterator::any so we can short-circuit on success
    if cell_sums.any(is_different_from_first) {
        Ok(())
    } else {
        Err(Error::NormalizedCounts)
    }
}

fn validate_total_counts_are_different_sparse(counts: &UmiCounts) -> Result<(), Error> {
    let mut total_counts_per_cell = counts
        .outer_iterator()
        .map(|cell| cell.data().iter().sum::<i32>());
    let Some(first) = total_counts_per_cell.next() else {
        return Err(Error::EmptyCounts);
    };

    if total_counts_per_cell.any(|s| s != first) {
        Ok(())
    } else {
        Err(Error::NormalizedCounts)
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

    use crate::reference_dataset::umi_counts::{RawUmiCounts, read_umi_counts_from_h5ad};

    impl RawUmiCounts {
        fn data(&self) -> &[i32] {
            self.0.data()
        }
    }

    #[test]
    fn read_sparse_h5ad_files() {
        // The last one apparently has some compression applied to it because of scanpy.write's default behavior
        let files = ["csr_adata", "csc_adata", "WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix"]
            .map(|fname| format!("test-data/{fname}.h5ad"))
            .map(|path| File::open(path).unwrap());

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
        }
    }

    #[test]
    fn read_dense_h5ad_file() {
        let filename = "test-data/dense_adata.h5ad";
        let file = File::open(filename).unwrap();
        let counts = read_umi_counts_from_h5ad(&file).unwrap();

        assert_eq!(
            counts.data()[0],
            10,
            "first entry in UMI counts of {filename} != 10"
        );
    }
}
