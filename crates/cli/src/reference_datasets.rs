use std::{
    collections::{HashMap, HashSet},
    str::FromStr,
};

use anyhow::{Context, anyhow, bail, ensure};
use camino::{Utf8Path, Utf8PathBuf};
use serde::Serialize;
use xenium_panel_convert_core::reference_dataset::{
    self,
    columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneNameCol},
    feature_set::FeatureSet,
    read_reference_dataset, write_reference_dataset,
};

use crate::write::write_to_file;

pub fn convert_reference_datasets(
    ReferenceDatasetCliOptions { reference_datasets }: &ReferenceDatasetCliOptions,
    output_dir: &Utf8Path,
) -> anyhow::Result<()> {
    let output_path = |filename: &str| output_dir.join(filename);

    for ReferenceDatasetSpec {
        path,
        cell_barcode_col,
        cell_annotation_col,
        ensembl_id_col,
        gene_name_col,
        transcriptome,
        rename,
    } in reference_datasets
    {
        let dataset_name = rename
            .as_deref()
            .and_then(Utf8Path::file_stem)
            .or(path.file_stem())
            .ok_or_else(|| anyhow!("invalid filename: {}", rename.as_ref().unwrap_or(path)))?;

        match read_reference_dataset(
            path,
            cell_barcode_col,
            cell_annotation_col,
            ensembl_id_col,
            gene_name_col,
            *transcriptome,
        ) {
            Ok(ds) => {
                write_reference_dataset(&output_path(dataset_name), &ds)?;
            }
            Err(errors) => {
                let errors = DatasetErrors {
                    path: path.to_owned(),
                    errors,
                };

                let error_path = format!("{dataset_name}-errors.json");
                write_to_file(&errors, &output_path(&error_path))?;
            }
        }
    }

    Ok(())
}

#[derive(Debug, Serialize)]
struct DatasetErrors {
    path: Utf8PathBuf,
    errors: Vec<reference_dataset::Error>,
}

#[derive(Clone, Debug)]
struct ReferenceDatasetSpec {
    path: Utf8PathBuf,
    cell_barcode_col: CellBarcodeCol,
    cell_annotation_col: CellAnnotationCol,
    ensembl_id_col: EnsemblIdCol,
    gene_name_col: GeneNameCol,
    transcriptome: FeatureSet,
    rename: Option<Utf8PathBuf>,
}

impl ReferenceDatasetSpec {
    fn parse_commandline(s: &str) -> anyhow::Result<Self> {
        const EXAMPLE: &str = "path=matrix.h5ad,barcode-col=barcodes,annotation-col=annotations,\
                               ensembl-id-col=gene_ids\nmatrix.h5ad,b=barcodes,a=annotations,\
                               e=gene_ids";

        fn get_spec_value_default<'a, T: Default + FromStr>(
            spec: &'a HashMap<&str, &str>,
            key: &str,
        ) -> anyhow::Result<T>
        where
            Result<T, T::Err>: Context<T, T::Err>,
        {
            let Some(val) = spec.get(key) else {
                return Ok(T::default());
            };

            T::from_str(*val)
                .with_context(|| format!("failed to parse value '{val}' for key '{key}'"))
        }

        fn get_spec_value<'a, T: FromStr>(
            spec: &'a HashMap<&str, &str>,
            key: &str,
        ) -> anyhow::Result<T>
        where
            Result<T, T::Err>: Context<T, T::Err>,
        {
            let val = spec
                .get(key)
                .ok_or_else(|| anyhow!("key '{key}' is required"))?;

            T::from_str(*val)
                .with_context(|| format!("failed to parse value '{val}' for key '{key}'"))
        }

        let key_aliases: HashMap<_, _> = [
            ("p", "path"),
            ("b", "barcode-col"),
            ("a", "annotation-col"),
            ("e", "ensembl-id-col"),
            ("g", "gene-name-col"),
            ("t", "transcriptome"),
            ("f", "flex"),
            ("r", "rename"),
        ]
        .clone()
        .into_iter()
        .collect();
        let allowed_keys: HashSet<_> = key_aliases
            .iter()
            .map(|(s1, s2)| [*s1, *s2])
            .flatten()
            .collect();

        let mut spec = HashMap::with_capacity(6);

        for (i, kv_pair) in s.split(',').enumerate() {
            let (key, value) = match (i, kv_pair.split_once('=')) {
                (_, Some((key, value))) => (key, value),
                (0, None) => ("path", kv_pair),
                (_, None) => {
                    bail!(
                        "reference dataset specification must be provided like one of the \
                         following (mixing aliases and full-names is allowed):\n{EXAMPLE}"
                    );
                }
            };

            ensure!(
                allowed_keys.contains(key),
                "key '{key}' not recognized in reference dataset specification"
            );

            let key = key_aliases.get(key).unwrap_or(&key);

            ensure!(
                spec.insert(*key, value).is_none(),
                "key '{key}' may not be specified more than once"
            );
        }

        Ok(Self {
            path: get_spec_value(&spec, "path")?,
            cell_barcode_col: get_spec_value_default(&spec, "barcode-col")?,
            cell_annotation_col: get_spec_value(&spec, "annotation-col")?,
            ensembl_id_col: get_spec_value_default(&spec, "ensembl-id-col")?,
            gene_name_col: get_spec_value_default(&spec, "gene-name-col")?,
            transcriptome: FeatureSet::new(
                get_spec_value(&spec, "transcriptome")?,
                get_spec_value_default(&spec, "flex")?,
            ),
            rename: get_spec_value(&spec, "rename").ok(),
        })
    }
}

#[derive(Clone, Debug, clap::Args)]
pub struct ReferenceDatasetCliOptions {
    #[clap(value_parser = ReferenceDatasetSpec::parse_commandline)]
    reference_datasets: Vec<ReferenceDatasetSpec>,
}
