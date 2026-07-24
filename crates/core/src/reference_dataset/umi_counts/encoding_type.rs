use std::{ops::Deref, str::FromStr};

use hdf5_metno::{Dataset, Group, Location, types::VarLenUnicode};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, strum::EnumString)]
#[strum(serialize_all = "snake_case")]
pub enum DenseEncodingType {
    Array,
}

impl DenseEncodingType {
    pub fn from_dataset(ds: &Dataset) -> Result<Self, super::Error> {
        encoding_type(&**ds)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, strum::EnumString)]
#[strum(serialize_all = "snake_case")]
pub enum SparseEncodingType {
    CsrMatrix,
    CscMatrix,
}

impl SparseEncodingType {
    pub fn from_group(group: &Group) -> Result<Self, super::Error> {
        encoding_type(group)
    }
}

fn encoding_type<T>(x: &impl Deref<Target = Location>) -> Result<T, super::Error>
where
    T: FromStr,
{
    let encoding_type: VarLenUnicode = x.attr("encoding-type").and_then(|a| a.read_scalar())?;

    T::from_str(encoding_type.as_str()).map_err(|_| super::Error::UnknownEncodingType)
}
