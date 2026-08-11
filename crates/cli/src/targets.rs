use std::{collections::HashMap, fs};

use anyhow::{Context, anyhow};
use camino::{Utf8Path, Utf8PathBuf};
use xenium_panel_validate_core::{
    Chemistry, Species,
    target_list::{
        ValidTarget, XeniumPanelDesignerGeneList,
        chemistry::{
            xenium_prime_human_ensembl_id_to_gene, xenium_prime_mouse_ensembl_id_to_gene,
            xenium_v1_human_ensembl_id_to_gene, xenium_v1_mouse_ensembl_id_to_gene,
        },
        parse_target_list,
    },
};

use crate::error::write_error_report;

pub fn convert_target_list(
    TargetListCliOptions {
        targets_path,
        field_alias_path,
        field_aliases,
        chemistry,
    }: &TargetListCliOptions,
    species: Species,
    output_dir: &Utf8Path,
) -> anyhow::Result<()> {
    let target_list = fs::read_to_string(targets_path)
        .with_context(|| format!("failed to read target-list from {targets_path}"))?;

    let field_alias_file_contents = field_alias_path
        .as_ref()
        .map(fs::read)
        .transpose()
        .with_context(|| {
            if let Some(field_alias_path) = field_alias_path {
                format!("failed to read field aliases from {field_alias_path}")
            } else {
                "failed to read field aliases".to_owned()
            }
        })?;

    let field_aliases = combine_field_aliases(field_alias_file_contents.as_deref(), field_aliases)
        .with_context(|| {
            if let Some(path) = field_alias_path {
                format!("failed to read field aliases from {path}")
            } else {
                "failed to construct field aliases".to_owned()
            }
        })?;

    let result = match (species, chemistry) {
        (Species::HomoSapiens, Chemistry::V1) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_v1_human_ensembl_id_to_gene,
        ),
        (Species::HomoSapiens, Chemistry::Prime) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_prime_human_ensembl_id_to_gene,
        ),
        (Species::MusMusculus, Chemistry::V1) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_v1_mouse_ensembl_id_to_gene,
        ),
        (Species::MusMusculus, Chemistry::Prime) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_prime_mouse_ensembl_id_to_gene,
        ),
    };

    match result {
        Ok(targets) => {
            write_valid_targets(&targets, &output_dir.join("validated-targets.csv"))?;
            write_targets_for_xenium_panel_designer(
                &XeniumPanelDesignerGeneList::from_valid_targets(targets),
                &output_dir.join("xenium-panel-designer-targets.csv"),
            )
        }
        Err(e) => write_error_report(&e, &output_dir.join("target-list.errors.json")),
    }
}

fn write_valid_targets(targets: &[ValidTarget], path: &Utf8Path) -> anyhow::Result<()> {
    todo!()
}

fn write_targets_for_xenium_panel_designer(
    targets: &XeniumPanelDesignerGeneList,
    path: &Utf8Path,
) -> anyhow::Result<()> {
    todo!()
}

#[derive(Debug, Clone, clap::Args)]
pub struct TargetListCliOptions {
    #[clap(long, short)]
    targets_path: Utf8PathBuf,
    #[clap(long, short = 'p')]
    field_alias_path: Option<Utf8PathBuf>,
    #[clap(long, short = 'a', value_parser = parse_field_aliases)]
    field_aliases: Vec<(String, String)>,
    #[clap(long, short)]
    chemistry: Chemistry,
}

fn parse_field_aliases(s: &str) -> anyhow::Result<(String, String)> {
    s.split_once('=')
        .map(|(alias, field)| (alias.to_owned(), field.to_owned()))
        .ok_or_else(|| anyhow!("field aliases must be specified as '<ALIAS>=<FIELD>'"))
}

fn combine_field_aliases<'a>(
    field_alias_file_contents: Option<&'a [u8]>,
    field_aliases: &'a [(String, String)],
) -> anyhow::Result<HashMap<&'a str, &'a str>> {
    let mut field_aliases: HashMap<_, _> = field_aliases
        .iter()
        .map(|(s1, s2)| (s1.as_str(), s2.as_str()))
        .collect();

    let Some(aliases_from_file) = field_alias_file_contents else {
        return Ok(field_aliases);
    };

    let aliases_from_file: HashMap<_, _> = toml::from_slice(aliases_from_file)?;

    for (alias, field) in aliases_from_file {
        // We want field-aliases from the CLI to take precedence
        field_aliases.entry(alias).or_insert(field);
    }

    Ok(field_aliases)
}

#[cfg(test)]
mod tests {
    use crate::targets::combine_field_aliases;

    #[test]
    fn field_aliases_are_combined_correctly() {
        let field_aliases = ["alias1", "field1", "alias2", "field2"];

        let field_aliases: Vec<(String, String)> = field_aliases
            .chunks(2)
            .map(|alias_field| (alias_field[0].to_owned(), alias_field[1].to_owned()))
            .collect();

        let field_aliases =
            combine_field_aliases(Some(br#"alias1 = "field2""#), &field_aliases).unwrap();

        assert_eq!(
            field_aliases,
            [("alias1", "field1"), ("alias2", "field2")]
                .into_iter()
                .collect()
        );
    }
}
