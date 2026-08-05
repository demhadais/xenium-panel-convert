use std::{
    fs,
    io::{self, Write},
};

use anyhow::Context;
use camino::Utf8Path;
use clap::Parser;
use serde::Serialize;
use xenium_panel_validate::{
    reference_datasets::{self, validate_reference_datasets},
    targets::{self, Chemistry, Species, parse_target_list_from_file},
};
use xenium_panel_validate_core::reference_dataset::write_reference_dataset;

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli {
        Cli::Targets {
            args,
            common:
                CommonOptions {
                    species,
                    chemistry,
                    output_format,
                    output,
                },
        } => {
            let parsed_targets = parse_target_list_from_file(&args, species, chemistry)?;

            match output_format {
                Format::Json => write_json_report(&parsed_targets, output.as_deref())?,
            }
        }
        Cli::References {
            args,
            common:
                CommonOptions {
                    species,
                    chemistry,
                    output_format,
                    output,
                },
        } => {
            let results = validate_reference_datasets(&args, species, chemistry)?;

            for (result, path) in results
                .iter()
                .zip(args.paths().iter().map(|p| p.with_extension("")))
            {
                match result {
                    Ok(ds) => write_reference_dataset(&path, ds)
                        .with_context(|| format!("failed to write dataset to {path}"))?,
                    Err(e) => match output_format {
                        Format::Json => write_json_report(e, output.as_deref())?,
                    },
                }
            }
        }
    }

    Ok(())
}

fn write_json_report(data: &impl Serialize, output_path: Option<&Utf8Path>) -> anyhow::Result<()> {
    if let Some(path) = output_path {
        let error_message = || format!("failed to write report to {path}");
        let file = fs::File::create(path).with_context(error_message)?;

        serde_json::to_writer(io::BufWriter::new(file), data).with_context(error_message)?;
    } else {
        let mut stdout = io::BufWriter::new(io::stdout().lock());

        serde_json::to_writer(&mut stdout, data)?;
        stdout
            .write_all(b"\n")
            .context("failed to write report to stdout")?;
    }

    Ok(())
}

#[derive(clap::Parser)]
enum Cli {
    Targets {
        #[clap(flatten)]
        args: targets::CommandlineArgs,
        #[clap(flatten)]
        common: CommonOptions,
    },
    References {
        #[clap(flatten)]
        args: reference_datasets::CommandlineArgs,
        #[clap(flatten)]
        common: CommonOptions,
    },
}

#[derive(Debug, Clone, clap::Args)]
struct CommonOptions {
    #[clap(long, short)]
    species: Species,
    #[clap(long, short)]
    chemistry: Chemistry,
    #[clap(long, short = 'f', default_value_t = Format::Json)]
    output_format: Format,
    #[clap(long, short)]
    output: Option<camino::Utf8PathBuf>,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
enum Format {
    Json,
}

impl std::fmt::Display for Format {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Json => "json".fmt(f),
        }
    }
}
