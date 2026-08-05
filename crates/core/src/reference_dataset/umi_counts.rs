use std::str::FromStr;

use hdf5_metno::{File, types::VarLenUnicode};
pub use matrix::RawCscUmiCounts;
use serde::Serialize;

use crate::reference_dataset::{
    h5::{FieldType, ReadFieldError, read_attribute, read_dataset_raw},
    umi_counts::encoding_type::{DenseEncodingType, EncodingType, SparseEncodingType},
};

mod encoding_type;
mod matrix;

pub fn read_umi_counts_from_h5ad(file: &File) -> Result<RawCscUmiCounts, Error> {
    let encoding_type: VarLenUnicode = read_attribute(file, "X/encoding-type")?;

    let encoding_type =
        EncodingType::from_str(&encoding_type).map_err(|()| Error::UnknownEncodingType {
            found: encoding_type.to_string(),
        })?;

    match encoding_type {
        EncodingType::Sparse(enc) => read_x_sparse(file, enc),
        EncodingType::Dense(enc) => read_x_dense(file, enc),
    }
}

fn read_x_sparse(file: &File, encoding_type: SparseEncodingType) -> Result<RawCscUmiCounts, Error> {
    let data = read_dataset_raw(file, "X/data")?;
    let indptr = read_dataset_raw(file, "X/indptr")?;
    let indices = read_dataset_raw(file, "X/indices")?;

    // It's very nice that scanpy decides to store the shape as an attribute rather
    // than the following 10x Genomics and storing it as a dataset. It's great when
    // a library built to analyze data ends up changing the format of the data :)
    let shape = file
        .attr("X/shape")
        .and_then(|sh| sh.read_1d())
        .map_err(|_| Error::MalformedCounts {
            error: ReadFieldError::DataTypeOrMissing {
                field_type: FieldType::Attribute,
                path: "X/shape".to_owned(),
            },
        })?;
    let shape = (shape[0], shape[1]);

    RawCscUmiCounts::from_sparse_matrix(shape, indptr, indices, data, encoding_type)
}

fn read_x_dense(file: &File, encoding_type: DenseEncodingType) -> Result<RawCscUmiCounts, Error> {
    let counts =
        file.dataset("X")
            .and_then(|ds| ds.read_2d())
            .map_err(|_| Error::MalformedMatrix {
                reason: "failed to read counts in dataset 'X' as 2D array".to_owned(),
            })?;

    RawCscUmiCounts::from_dense_matrix(&counts, encoding_type)
}

#[derive(Debug, thiserror::Error, Clone, Serialize)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum Error {
    #[error(transparent)]
    MalformedCounts {
        #[from]
        error: ReadFieldError,
    },
    #[error(
        "unknown encoding type: '{found}'. Expected one of: {:?}",
        EncodingType::VARIANTS
    )]
    UnknownEncodingType { found: String },
    #[error("empty counts")]
    EmptyCounts,
    #[error("transformed counts")]
    TransformedCounts,
    #[error("normalized counts")]
    NormalizedCounts,
    #[error("malformed matrix: {reason}")]
    MalformedMatrix { reason: String },
}

impl From<sprs::errors::StructureError> for Error {
    fn from(err: sprs::errors::StructureError) -> Self {
        Self::MalformedMatrix {
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
        // The last one apparently has some compression applied to it because of
        // scanpy.write's default behavior. It's really nice that scanpy's default
        // behavior differs from anndata's default behavior :)
        let files = ["csr_adata", "csc_adata", "dense_adata", "10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx_count_sample_filtered_feature_bc_matrix"]
            .map(|fname| format!("test-data/{fname}.h5ad"))
            .map(|path| File::open(path).unwrap());

        let mut all_counts = Vec::with_capacity(files.len());

        for f in files {
            let filename = f.filename();
            let counts = read_umi_counts_from_h5ad(&f)
                .unwrap_or_else(|e| panic!("failed to read UMI counts from {filename}: {e}"));

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
