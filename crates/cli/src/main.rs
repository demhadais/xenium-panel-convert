use clap::Parser;
use xenium_panel_validate::{
    reference_datasets::{ReferenceDatasetCliOptions, convert_reference_datasets},
    targets::{TargetListCliOptions, convert_target_list},
};
use xenium_panel_validate_core::Species;

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli {
        Cli::Targets {
            common:
                CommonOptions {
                    species,
                    output_dir,
                },
            options,
        } => convert_target_list(&options, species, &output_dir),
        Cli::References {
            common:
                CommonOptions {
                    species,
                    output_dir,
                },
            options,
        } => convert_reference_datasets(&options, species, &output_dir),
        Cli::Submission {
            common,
            targets_options,
            references_options,
        } => todo!(),
    }
}

#[derive(clap::Parser)]
enum Cli {
    /// Convert a target list to a format suitable for the Xenium Panel Designer.
    ///
    /// The target-list must be a CSV-file with the header: "ensembl_id,gene_name,group,priority". The column "priority" must be one of "must_have", "desired", or "backup". The file will be converted to a sorted CSV-file suitable for copy-pasting into the Xenium Panel Designer. If errors are encountered, they are collected and written to <OUTPUT_DIR>/target-list-errors.json.
    Targets {
        #[clap(flatten)]
        common: CommonOptions,
        #[clap(flatten)]
        options: TargetListCliOptions,
    },
    /// Convert scanpy-annotated single-cell RNA sequencing datasets to a format suitable for the Xenium Panel Designer.
    ///
    /// Each annotated dataset (generated with scanpy, typically with the file-extention .h5ad) is converted to a directory containing a matrix.h5 and annotations.csv. This directory can be archived, zipped, and uploaded directly to the panel designer as a reference dataset. If errors are encountered, they are written to <OUTPUT_DIR>/<DATASET_PATH>-errors.json
    References {
        #[clap(flatten)]
        common: CommonOptions,
        #[clap(flatten)]
        options: ReferenceDatasetCliOptions,
    },
    /// Validate both the target-list and the reference datasets and convert both to the appropriate formats for the Xenium Panel Designer
    Submission {
        #[clap(flatten)]
        common: CommonOptions,
        #[clap(flatten)]
        targets_options: TargetListCliOptions,
        #[clap(flatten)]
        references_options: ReferenceDatasetCliOptions,
    },
}

#[derive(Clone, clap::Args)]
struct CommonOptions {
    #[clap(long, short)]
    species: Species,
    #[clap(long, short)]
    output_dir: camino::Utf8PathBuf,
}
