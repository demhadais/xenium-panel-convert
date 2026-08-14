use std::{
    fs,
    io::{self, Write},
    path::Path,
};

use anyhow::Context;

pub(crate) fn write_map_to_file(
    path: &Path,
    map_name: &str,
    map: &phf_codegen::Map<'_, &str>,
) -> anyhow::Result<()> {
    let file = fs::File::create(path)
        .with_context(|| format!("failed to write file {}", path.to_str().unwrap()))?;
    let mut file_writer = io::BufWriter::new(file);

    writeln!(
        file_writer,
        "pub(super) static {map_name}: phf::Map<&'static str, &'static str> = {};",
        map.build()
    )?;

    Ok(())
}
