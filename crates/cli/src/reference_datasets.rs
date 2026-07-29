use xenium_panel_validate_core::reference_dataset::{
    self, ReferenceDataset, validate_reference_dataset,
};

#[allow(clippy::missing_errors_doc)]
pub fn validate_reference_datasets(
    CommandlineArgs {
        paths,
        cell_annotations_cols,
        ensembl_id_cols,
        gene_name_cols,
    }: &CommandlineArgs,
) -> anyhow::Result<Vec<Result<ReferenceDataset, DatasetErrors>>> {
    anyhow::ensure!(
        paths.len() == cell_annotations_cols.len()
            && cell_annotations_cols.len() == ensembl_id_cols.len()
            && ensembl_id_cols.len() == gene_name_cols.len(),
        "must provide the same number of dataset paths, cell-annotations columns, Ensembl ID columns, and gene-name columns"
    );

    let mut results = Vec::with_capacity(8);

    for (((path, cell_annotations_col), ensembl_id_col), gene_name_col) in paths
        .iter()
        .zip(cell_annotations_cols)
        .zip(ensembl_id_cols)
        .zip(gene_name_cols)
    {
        match validate_reference_dataset(path, cell_annotations_col, ensembl_id_col, gene_name_col)
        {
            Ok(ds) => results.push(Ok(ds)),
            Err(errors) => results.push(Err(DatasetErrors {
                path: path.to_owned(),
                errors,
            })),
        };
    }

    Ok(results)
}

#[derive(Clone, Debug, PartialEq, Eq, clap::Args)]
pub struct CommandlineArgs {
    #[clap(short, long)]
    paths: Vec<String>,
    #[clap(short, long)]
    cell_annotations_cols: Vec<String>,
    #[clap(short, long)]
    ensembl_id_cols: Vec<String>,
    #[clap(short, long)]
    gene_name_cols: Vec<String>,
}

#[derive(Debug)]
pub struct DatasetErrors {
    path: String,
    errors: Vec<reference_dataset::Error>,
}
