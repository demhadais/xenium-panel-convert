use hdf5_metno::{File, types::VarLenUnicode};
use ndarray::Array1;

pub fn read_string_array_from_file(
    file: &File,
    path: &str,
) -> hdf5_metno::Result<Array1<VarLenUnicode>> {
    file.dataset(path).and_then(|ds| ds.read_1d())
}
