use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, bail, ensure};
use camino::{Utf8Path, Utf8PathBuf};
use serde::Serialize;
use xenium_panel_validate_core::{
    Species,
    reference_dataset::{
        self,
        columns::{CellAnnotationCol, CellBarcodeCol, EnsemblIdCol, GeneNameCol},
        read_reference_dataset,
    },
};

use crate::error::write_error_report;

pub fn convert_reference_datasets(
    CliOptions {
        reference: reference_datasets,
    }: &CliOptions,
    species: Species,
    output_dir: &Utf8Path,
) -> anyhow::Result<()> {
    let mut parsed_datasets = Vec::with_capacity(reference_datasets.len());
    let mut errors = Vec::with_capacity(reference_datasets.len());

    for ReferenceDatasetSpecification {
        path,
        cell_barcode_col,
        cell_annotation_col,
        ensembl_id_col,
        gene_name_col,
    } in reference_datasets
    {
        match read_reference_dataset(
            path,
            cell_barcode_col,
            cell_annotation_col,
            ensembl_id_col,
            gene_name_col,
            species,
        ) {
            Ok(ds) => {
                parsed_datasets.push(ds);
            }
            Err(errs) => {
                errors.push(DatasetErrors {
                    path: path.to_owned(),
                    errors: errs,
                });
            }
        }
    }

    for _ds in &parsed_datasets {
        todo!()
    }

    for error_set in &errors {
        write_error_report(
            error_set,
            &output_dir.join(
                error_set
                    .path
                    .file_name()
                    .expect("a filename shouldn't terminate in '..'"),
            ),
        )?;
    }

    Ok(())
}

#[derive(Debug, Serialize)]
pub struct DatasetErrors {
    path: Utf8PathBuf,
    errors: Vec<reference_dataset::Error>,
}

#[derive(Clone, Debug)]
struct ReferenceDatasetSpecification {
    path: Utf8PathBuf,
    cell_barcode_col: CellBarcodeCol,
    cell_annotation_col: CellAnnotationCol,
    ensembl_id_col: EnsemblIdCol,
    gene_name_col: GeneNameCol,
}

impl ReferenceDatasetSpecification {
    fn parse_commandline(s: &str) -> anyhow::Result<Self> {
        const EXAMPLE: &str = "xp-prep --reference \
                               'path=matrix.h5ad,barcode_col=barcodes,annotation_col=annotations,\
                               ensembl_id_col=gene_ids\nxp-prep --r \
                               'matrix.h5ad,b=barcodes,a=annotations,e=gene_ids";

        fn get_spec_value_default<'a, T: Default + From<&'a str>>(
            spec: &'a HashMap<&str, &str>,
            key: &str,
        ) -> T {
            spec.get(key).map(|v| T::from(*v)).unwrap_or_default()
        }

        fn get_spec_value<'a, T: From<&'a str>>(
            spec: &'a HashMap<&str, &str>,
            key: &str,
        ) -> anyhow::Result<T> {
            spec.get(key).map(|v| T::from(*v)).ok_or_else(|| {
                anyhow!(
                    "{key} not found - specify key-value pairs as one of the following (mixing \
                     allowed):\n{EXAMPLE}"
                )
            })
        }

        let key_aliases: HashMap<_, _> = [
            ("p", "path"),
            ("b", "barcode_col"),
            ("a", "annotation_col"),
            ("e", "ensembl_id_col"),
            ("g", "gene_name_col"),
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
                        "reference dataset specification must be provided as one of the following \
                         (mixing allowed):\n{EXAMPLE}"
                    );
                }
            };

            ensure!(
                allowed_keys.contains(key),
                "key {key} not recognized in reference dataset specification"
            );

            let key = key_aliases.get(key).unwrap_or(&key);

            ensure!(
                spec.insert(*key, value).is_none(),
                "{key} may not be specified more than once"
            );
        }

        Ok(Self {
            path: get_spec_value(&spec, "path")?,
            cell_barcode_col: get_spec_value_default(&spec, "barcode_col"),
            cell_annotation_col: get_spec_value(&spec, "annotation_col")?,
            ensembl_id_col: get_spec_value_default(&spec, "ensembl_id_col"),
            gene_name_col: get_spec_value_default(&spec, "gene_name_col"),
        })
    }
}

#[derive(Clone, Debug, clap::Args)]
pub struct CliOptions {
    #[clap(long, value_parser = ReferenceDatasetSpecification::parse_commandline, help = "Reference dataset specification.")]
    reference: Vec<ReferenceDatasetSpecification>,
}
