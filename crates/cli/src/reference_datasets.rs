use camino::{Utf8Path, Utf8PathBuf};
use xenium_panel_validate_core::reference_dataset::{
    self, ReferenceDataset, validate_reference_dataset,
};

use crate::targets::{Chemistry, Species};

#[allow(clippy::missing_errors_doc)]
pub fn validate_reference_datasets(
    CommandlineArgs {
        paths,
        cell_barcode_cols,
        cell_annotation_cols,
        ensembl_id_cols,
        gene_name_cols,
    }: &CommandlineArgs,
    species: Species,
    chemistry: Chemistry,
) -> anyhow::Result<Vec<Result<ReferenceDataset, DatasetErrors>>> {
    anyhow::ensure!(
        paths.len() == cell_barcode_cols.len()
            && cell_barcode_cols.len() == cell_annotation_cols.len()
            && cell_annotation_cols.len() == ensembl_id_cols.len()
            && ensembl_id_cols.len() == gene_name_cols.len(),
        "the number of dataset paths, cell-barcode columns, cell-annotation columns, Ensembl ID columns, and gene-name columns must all be equal"
    );

    let mut results = Vec::with_capacity(paths.len());

    for ((((path, cell_barcode_col), cell_annotation_col), ensembl_id_col), gene_name_col) in paths
        .iter()
        .zip(cell_barcode_cols)
        .zip(cell_annotation_cols)
        .zip(ensembl_id_cols)
        .zip(gene_name_cols)
    {
        match validate_reference_dataset(
            path,
            cell_barcode_col,
            cell_annotation_col,
            ensembl_id_col,
            gene_name_col,
        ) {
            Ok(ds) => results.push(Ok(ds)),
            Err(errors) => results.push(Err(DatasetErrors {
                path: path.to_owned(),
                errors,
            })),
        };
    }

    Ok(results)
}

pub fn write_reference_datasets() {}

#[derive(Clone, Debug, PartialEq, Eq, clap::Args)]
pub struct CommandlineArgs {
    #[clap(short, long)]
    paths: Vec<Utf8PathBuf>,
    #[clap(short = 'b', long)]
    cell_barcode_cols: Vec<String>,
    #[clap(short = 'a', long)]
    cell_annotation_cols: Vec<String>,
    #[clap(short, long)]
    ensembl_id_cols: Vec<String>,
    #[clap(short, long)]
    gene_name_cols: Vec<String>,
}

#[derive(Debug)]
pub struct DatasetErrors {
    path: Utf8PathBuf,
    errors: Vec<reference_dataset::Error>,
}
