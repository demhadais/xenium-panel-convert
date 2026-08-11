use std::{fs, io};

use anyhow::Context;
use camino::Utf8Path;
use serde::Serialize;

pub fn write_error_report(data: &impl Serialize, path: &Utf8Path) -> anyhow::Result<()> {
    let error_message = || format!("failed to write error report to {path}");
    let file = fs::File::create(path).with_context(error_message)?;

    serde_json::to_writer(io::BufWriter::new(file), data).with_context(error_message)?;

    Ok(())
}
