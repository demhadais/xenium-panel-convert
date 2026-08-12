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
    /// Convert a target-list to a format suitable for the Xenium Panel
    /// Designer.
    ///
    /// The target-list must be a CSV-file with the header:
    /// "ensembl_id,gene_name,group,priority". The column "priority" must be one
    /// of "must_have", "desired", or "backup". A "cleaned" version of the file
    /// will be saved at <OUTPUT_DIR>/target-list.csv, and the version for the
    /// panel designer will be saved at
    /// <OUTPUT_DIR>/xenium-panel-designer-target-list.csv. If errors occur,
    /// they are collected and written to <OUTPUT_DIR>/target-list-errors.json.
    Targets {
        #[clap(flatten)]
        common: CommonOptions,
        #[clap(flatten)]
        options: TargetListCliOptions,
    },
    /// Convert scanpy-annotated single-cell RNA sequencing datasets to a format
    /// suitable for the Xenium Panel Designer.
    ///
    /// Each annotated dataset (h5ad files generated with scanpy) is converted
    /// to a directory containing a matrix.h5 and annotations.csv. This
    /// directory can be fed to tar to create a compressed archive and uploaded
    /// directly to the panel designer as a reference dataset. If errors occur,
    /// they are collected and written to
    /// <OUTPUT_DIR>/<DATASET_PATH>-errors.json
    References {
        #[clap(flatten)]
        common: CommonOptions,
        #[clap(flatten)]
        options: ReferenceDatasetCliOptions,
    },
    /// Convert both a target-list and reference datasets to formats suitable
    /// for the Xenium Panel Designer.
    ///
    /// This is the equivalent of running both commands at once with additional
    /// checks that are only possible with both a target-list and a reference
    /// dataset. Namely
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
