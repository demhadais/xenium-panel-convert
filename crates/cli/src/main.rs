use std::fs;

use anyhow::Context;
use camino::Utf8Path;
use clap::Parser;
use xenium_panel_validate::{
    reference_datasets::{self, validate_reference_datasets},
    targets::{self, parse_target_list_from_file},
};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli {
        Cli::Targets {
            args,
            common:
                CommonOptions {
                    output_format,
                    output,
                },
        } => {
            let parsed_targets = parse_target_list_from_file(&args)?;

            let parsed_targets = match output_format {
                Format::Json => serde_json::to_string(&parsed_targets)?,
            };

            write_report(&parsed_targets, output.as_deref())?;
        }
        Cli::ReferenceDatasets { args, common } => {
            let results = validate_reference_datasets(&args)?;
            // REMOVE!
            for r in results {
                r.unwrap();
            }
        }
    }

    Ok(())
}

fn write_report(data: &str, output_path: Option<&Utf8Path>) -> anyhow::Result<()> {
    if let Some(path) = output_path {
        fs::write(path, data).context(format!("failed to write report to {path}"))?;
    } else {
        println!("{data}");
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
    ReferenceDatasets {
        #[clap(flatten)]
        args: reference_datasets::CommandlineArgs,
        #[clap(flatten)]
        common: CommonOptions,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, clap::Args)]
struct CommonOptions {
    #[clap(long, short = 'f', default_value_t = Format::Json)]
    output_format: Format,
    #[clap(long, short)]
    output: Option<camino::Utf8PathBuf>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
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
