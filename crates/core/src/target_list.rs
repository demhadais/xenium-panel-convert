use std::collections::{HashMap, HashSet};

use chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId};
pub use error::Error;
pub(crate) use error::ErrorInner;
pub(crate) use target::{Priority as TargetPriority, ValidGene, ValidTarget};
pub use xenium_panel_designer::XeniumPanelDesignerGeneList;

use crate::target_list::{
    csv_util::{extract_record, read_csv_trimmed, rename_fields},
    target::UnvalidatedTarget,
};

pub mod chemistry;
mod csv_util;
mod error;
mod target;
mod xenium_panel_designer;

#[allow(clippy::missing_errors_doc)]
pub fn parse_target_list(
    target_list: &str,
    field_aliases: &HashMap<&str, &str>,
    ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)> + Copy,
) -> Result<Vec<ValidTarget>, Vec<Error>> {
    const N_GENES: usize = 500;

    let mut reader = read_csv_trimmed(target_list);
    // If we can't get headers, just return early
    let headers = reader.headers().map_err(|e| {
        vec![Error {
            errors: vec![e.into()],
            line_number: None,
            submitted_target: None,
        }]
    })?;

    // We initialize the list of errors from the field-renaming, but it doesn't
    // prevent us from continuing the parsing
    let (fieldnames, error) = rename_fields(headers, field_aliases);
    let mut errors = error.map(|e| vec![e]).unwrap_or_default();

    let mut valid_targets = Vec::with_capacity(N_GENES);
    let mut seen_genes = HashSet::with_capacity(N_GENES);

    for record in reader.records() {
        let Some(record) = extract_record(record.as_ref(), &mut errors) else {
            continue;
        };

        let line_number = record.position().map(csv::Position::line);

        let submitted_target = UnvalidatedTarget::from_record(record, &fieldnames);
        let validation_result = submitted_target.validate(ensembl_id_to_gene);

        push_validation_result(
            line_number,
            submitted_target,
            validation_result,
            &mut seen_genes,
            &mut valid_targets,
            &mut errors,
        );
    }

    if !errors.is_empty() {
        return Err(errors);
    }

    Ok(valid_targets)
}

fn push_validation_result(
    line_number: Option<u64>,
    submitted_target: UnvalidatedTarget,
    validation_result: Result<ValidTarget, Vec<ErrorInner>>,
    seen_genes: &mut HashSet<ValidGene>,
    valid_targets: &mut Vec<ValidTarget>,
    errors: &mut Vec<Error>,
) {
    let submitted_target = Some(submitted_target);

    match validation_result {
        Ok(vt) => {
            let is_new = seen_genes.insert(vt.gene);

            if is_new {
                valid_targets.push(vt);
            } else {
                errors.push(Error {
                    line_number,
                    submitted_target,
                    errors: vec![ErrorInner::DuplicateGene],
                });
            }
        }
        Err(errs) => errors.push(Error {
            line_number,
            submitted_target,
            errors: errs,
        }),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::target_list::{
        ErrorInner,
        chemistry::{tests::tp53_ensembl_id, xenium_v1_human_ensembl_id_to_gene},
        parse_target_list,
    };

    #[test]
    fn duplicate_targets() {
        let ensembl_id = tp53_ensembl_id();
        let ensembl_id_str = ensembl_id.as_str();

        // Two rows with the same Ensembl ID/gene-name pair but differing other fields
        let gene_list = format!(
            "ensembl_id,gene_name,group,priority\n{ensembl_id_str},TP53,group0,must_have\\
             n{ensembl_id_str},TP53,group1,backup"
        );

        let errors = parse_target_list(
            &gene_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(errors.len(), 1, "did not find exactly 1 error: {errors:?}");
        assert_eq!(errors[0].errors, [ErrorInner::DuplicateGene]);
    }

    #[test]
    fn error_reports_correct_file_line_number() {
        let gene_list = "ensembl_id,gene_name,group,priority\nid,gene,0,must_have";

        let errors = parse_target_list(
            gene_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(errors.len(), 1);
        assert_eq!(errors[0].line_number, Some(2));
    }
}
