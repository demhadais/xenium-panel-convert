use clap::Parser;
use xenium_panel_validate::{
    reference_datasets::{self, convert_reference_datasets},
    targets::{self, convert_target_list},
};
use xenium_panel_validate_core::{Chemistry, Species};

fn main() -> anyhow::Result<()> {
    let Cli {
        species,
        chemistry,
        output_dir,
        options:
            CliOptions {
                target_options,
                reference_dataset_options,
            },
    } = Cli::parse();

    match (target_options, reference_dataset_options) {
        (Some(opts), None) => {
            convert_target_list(&opts, species, chemistry, &output_dir)?;
        }
        (None, Some(opts)) => {
            convert_reference_datasets(&opts, species, &output_dir)?;
        }
        (None, None) => unreachable!(),
        (Some(_), Some(_)) => todo!(),
    }

    Ok(())
}

#[derive(clap::Parser)]
struct Cli {
    #[clap(long, short)]
    species: Species,
    #[clap(long, short)]
    chemistry: Chemistry,
    #[clap(long, short)]
    output_dir: camino::Utf8PathBuf,
    #[clap(flatten)]
    options: CliOptions,
}

#[derive(clap::Args)]
#[group(required = true)]
struct CliOptions {
    #[clap(flatten)]
    target_options: Option<targets::CliOptions>,
    #[clap(flatten)]
    reference_dataset_options: Option<reference_datasets::CliOptions>,
}
