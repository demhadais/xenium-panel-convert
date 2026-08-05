use std::str::FromStr;

use hdf5_metno::{
    File, H5Type,
    types::{FixedAscii, VarLenUnicode},
};
use ndarray::Array1;
use serde::Serialize;
use strum::VariantNames;

pub fn read_attribute<T: H5Type>(file: &File, path: &str) -> Result<T, ReadFieldError> {
    file.attr(path)
        .and_then(|a| a.read_scalar())
        .map_err(|_| ReadFieldError::DataTypeOrMissing {
            field_type: FieldType::Attribute,
            path: path.to_owned(),
        })
}

pub fn read_1d_dataset<T: H5Type>(file: &File, path: &str) -> Result<Array1<T>, ReadFieldError> {
    file.dataset(path)
        .and_then(|ds| ds.read_1d())
        .map_err(|_| ReadFieldError::DataTypeOrMissing {
            field_type: FieldType::Dataset,
            path: path.to_owned(),
        })
}

pub fn read_dataset_raw<T: H5Type>(file: &File, path: &str) -> Result<Vec<T>, ReadFieldError> {
    file.dataset(path)
        .and_then(|ds| ds.read_raw())
        .map_err(|_| ReadFieldError::DataTypeOrMissing {
            field_type: FieldType::Dataset,
            path: path.to_owned(),
        })
}

pub fn read_1d_string_dataset(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadFieldError> {
    let encoding_type: VarLenUnicode = read_attribute(file, &format!("{path}/encoding-type"))?;

    let encoding_type = StringEncodingType::from_str(&encoding_type).map_err(|_| {
        ReadFieldError::UnknownEncodingType {
            path: path.to_owned(),
            found: encoding_type.to_string(),
            expected: StringEncodingType::VARIANTS,
        }
    })?;

    match encoding_type {
        StringEncodingType::Categorical => read_categorical_array(file, path),
        StringEncodingType::StringArray => read_string_array(file, path),
        StringEncodingType::NullableStringArray => read_nullable_string_array(file, path),
    }
}

pub fn to_ascii<const N: usize>(s: &VarLenUnicode) -> FixedAscii<N> {
    FixedAscii::from_ascii(&s).expect("all strings are ASCII in this context")
}

fn read_categorical_array(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadFieldError> {
    let codes = read_1d_dataset::<i32>(file, &format!("{path}/codes"))?;
    let categories = read_1d_dataset::<VarLenUnicode>(file, &format!("{path}/categories"))?;

    codes
        .iter()
        .enumerate()
        .map(|(i, code)| {
            if *code == -1 {
                return Err(ReadFieldError::NullValue {
                    index: i,
                    dataset_path: path.to_owned(),
                });
            }

            Ok(categories[*code as usize].clone())
        })
        .collect()
}

fn read_string_array(file: &File, path: &str) -> Result<Array1<VarLenUnicode>, ReadFieldError> {
    read_1d_dataset(file, path)
}

fn read_nullable_string_array(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadFieldError> {
    let is_null_array = read_1d_dataset::<bool>(file, &format!("{path}/mask"))?;
    if let Some((index, _)) = is_null_array
        .iter()
        .enumerate()
        .find(|(_, is_null)| **is_null)
    {
        return Err(ReadFieldError::NullValue {
            index,
            dataset_path: path.to_owned(),
        });
    }

    read_string_array(file, &format!("{path}/values"))
}

#[derive(Clone, Copy, Debug, strum::EnumString, strum::VariantNames)]
#[strum(serialize_all = "kebab-case")]
pub enum StringEncodingType {
    Categorical,
    StringArray,
    NullableStringArray,
}

#[derive(Debug, thiserror::Error, Clone, Serialize)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum ReadFieldError {
    #[error(
        "failed to read {path} as {field_type} - it may be nonexistent or a different type than expected"
    )]
    DataTypeOrMissing { field_type: FieldType, path: String },
    #[error("missing value at index {index} in dataset {dataset_path}")]
    NullValue { index: usize, dataset_path: String },
    #[error("failed to read {path} due to unknown encoding-type {found}")]
    UnknownEncodingType {
        path: String,
        found: String,
        expected: &'static [&'static str],
    },
}

#[derive(Debug, Clone, Serialize, strum::Display)]
#[serde(tag = "type", rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum FieldType {
    Attribute,
    Dataset,
    Group,
}
