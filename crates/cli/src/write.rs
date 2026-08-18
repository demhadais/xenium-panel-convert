use std::{fs, io};

use anyhow::{Context, ensure};
use camino::Utf8Path;
use serde::Serialize;

pub(crate) fn write_to_file(data: &impl Serialize, path: &Utf8Path) -> anyhow::Result<()> {
    ensure!(!path.exists(), "cannot overwrite file at {path}");

    let error_message = || format!("failed to write error report to {path}");
    let file = fs::File::create(path).with_context(error_message)?;

    serde_json::to_writer(io::BufWriter::new(file), data).with_context(error_message)?;

    Ok(())
}
