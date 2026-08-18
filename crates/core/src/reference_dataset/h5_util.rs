use std::str::FromStr;

use hdf5_metno::{
    Container, File, Group, H5Type,
    types::{FixedAscii, VarLenUnicode},
};
use ndarray::{Array1, ArrayView, Dimension};
use serde::Serialize;
use strum::VariantNames;

pub(crate) fn read_container(file: &File, path: &str) -> Result<Container, ReadH5FieldError> {
    let container = match file.group(path) {
        Ok(g) => g.as_container().unwrap(),
        Err(_) => file
            .dataset(path)
            .and_then(|ds| ds.as_container())
            .map_err(|_| ReadH5FieldError::DataTypeOrMissing {
                field_type: FieldType::Container,
                object_path: path.to_owned(),
            })?,
    };

    Ok(container)
}

pub(crate) fn read_attribute<T: H5Type>(
    container: &Container,
    path: &str,
) -> Result<T, ReadH5FieldError> {
    container
        .attr(path)
        .and_then(|a| a.read_scalar())
        .map_err(|_| ReadH5FieldError::DataTypeOrMissing {
            field_type: FieldType::Attribute,
            object_path: path.to_owned(),
        })
}

pub(crate) fn read_1d_dataset<T: H5Type>(
    file: &File,
    path: &str,
) -> Result<Array1<T>, ReadH5FieldError> {
    file.dataset(path).and_then(|ds| ds.read_1d()).map_err(|_| {
        ReadH5FieldError::DataTypeOrMissing {
            field_type: FieldType::Dataset,
            object_path: path.to_owned(),
        }
    })
}

pub(crate) fn read_dataset_raw<T: H5Type>(
    file: &File,
    path: &str,
) -> Result<Vec<T>, ReadH5FieldError> {
    file.dataset(path)
        .and_then(|ds| ds.read_raw())
        .map_err(|_| ReadH5FieldError::DataTypeOrMissing {
            field_type: FieldType::Dataset,
            object_path: path.to_owned(),
        })
}

pub(crate) fn read_1d_string_dataset(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    let encoding_type: VarLenUnicode =
        read_attribute(&read_container(file, path)?, "encoding-type")?;

    let encoding_type = StringEncodingType::from_str(&encoding_type).map_err(|_| {
        ReadH5FieldError::UnknownEncodingType {
            object_path: path.to_owned(),
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

fn read_categorical_array(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    let codes = read_1d_dataset::<i32>(file, &format!("{path}/codes"))?;
    let categories = read_1d_dataset::<VarLenUnicode>(file, &format!("{path}/categories"))?;

    codes
        .iter()
        .enumerate()
        .map(|(i, code)| {
            if *code == -1 {
                return Err(ReadH5FieldError::NullValue {
                    index: i,
                    object_path: path.to_owned(),
                });
            }

            #[allow(clippy::cast_sign_loss)]
            Ok(categories[*code as usize].clone())
        })
        .collect()
}

fn read_string_array(file: &File, path: &str) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    read_1d_dataset(file, path)
}

fn read_nullable_string_array(
    file: &File,
    path: &str,
) -> Result<Array1<VarLenUnicode>, ReadH5FieldError> {
    let is_null_array = read_1d_dataset::<bool>(file, &format!("{path}/mask"))?;
    if let Some((index, _)) = is_null_array
        .iter()
        .enumerate()
        .find(|(_, is_null)| **is_null)
    {
        return Err(ReadH5FieldError::NullValue {
            index,
            object_path: path.to_owned(),
        });
    }

    read_string_array(file, &format!("{path}/values"))
}

pub(crate) fn to_ascii<const N: usize>(s: &VarLenUnicode) -> FixedAscii<N> {
    FixedAscii::from_ascii(&s).expect("all strings are ASCII in this context")
}

pub(crate) fn create_h5_group(file: &File, path: &str) -> Result<Group, CreateH5GroupError> {
    file.create_group(path).map_err(|e| CreateH5GroupError {
        object_path: path.to_owned(),
        reason: e.to_string(),
    })
}

pub(crate) fn write_dataset_to_h5_group<'d, A, T, D>(
    group: &Group,
    path: &str,
    data: A,
) -> Result<(), WriteH5DatasetError>
where
    A: Into<ArrayView<'d, T, D>>,
    T: H5Type,
    D: Dimension,
{
    group
        .new_dataset_builder()
        .with_data(data)
        .create(path)
        .map_err(|e| WriteH5DatasetError {
            object_path: path.to_owned(),
            reason: e.to_string(),
        })?;

    Ok(())
}

#[derive(Clone, Copy, Debug, strum::EnumString, strum::VariantNames)]
#[strum(serialize_all = "kebab-case")]
enum StringEncodingType {
    Categorical,
    StringArray,
    NullableStringArray,
}

#[derive(Debug, Clone, Serialize, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum ReadH5FieldError {
    #[error("wrong data type for or missing {field_type} - {object_path}")]
    DataTypeOrMissing {
        field_type: FieldType,
        object_path: String,
    },
    #[error("null value at index {index} - {object_path}")]
    NullValue { index: usize, object_path: String },
    #[error("unknown encoding type {found} at {object_path} - expected one of {expected:?}")]
    UnknownEncodingType {
        object_path: String,
        found: String,
        expected: &'static [&'static str],
    },
}

#[derive(Debug, Clone, Copy, Serialize, strum::Display)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum FieldType {
    Attribute,
    Container,
    Dataset,
}

#[derive(Debug, Clone, Serialize, thiserror::Error)]
#[error("failed to create H5 group {object_path} - {reason}")]
pub struct CreateH5GroupError {
    pub object_path: String,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, thiserror::Error)]
#[error("failed to write H5 dataset {object_path} - {reason}")]
pub struct WriteH5DatasetError {
    pub object_path: String,
    pub reason: String,
}
