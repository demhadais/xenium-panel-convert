use serde::{Deserialize, Serialize};

use crate::gene_list::{
    ErrorInner,
    chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName},
};

pub fn validate_target(
    UnvalidatedTarget {
        gene,
        group,
        is_backup,
        must_have,
    }: &UnvalidatedTarget,
    ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<ValidTarget, Vec<ErrorInner>> {
    let mut errors = Vec::new();

    if group.is_none() {
        errors.push(ErrorInner::MissingField { fieldname: "group" });
    }

    let is_backup = match parse_bool_from_str(is_backup.as_deref(), "is_backup").map(IsBackup) {
        Ok(is_backup) => Some(is_backup),
        Err(err) => {
            errors.push(err);
            None
        }
    };

    let must_have = match parse_bool_from_str(must_have.as_deref(), "must_have").map(MustHave) {
        Ok(must_have) => Some(must_have),
        Err(err) => {
            errors.push(err);
            None
        }
    };

    if let (Some(is_backup), Some(must_have)) = (is_backup, must_have)
        && is_backup.0
        && must_have.0
    {
        errors.push(ErrorInner::BackupAndMustHave);
    }

    let valid_gene = match validate_ensembl_id_gene_name_pair(gene, ensembl_id_to_gene) {
        Ok(vg) => Some(vg),
        Err(err) => {
            errors.push(err);
            None
        }
    };

    match (valid_gene, group, is_backup, must_have) {
        (Some(valid_gene), Some(group), Some(is_backup), Some(must_have)) if errors.is_empty() => {
            Ok(ValidTarget {
                gene: valid_gene,
                group: group.to_ascii_lowercase(),
                is_backup,
                must_have,
            })
        }
        _ => Err(errors),
    }
}

pub fn validate_ensembl_id_gene_name_pair(
    unvalidated_gene: &UnvalidatedGene,
    ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<ValidGene, ErrorInner> {
    let UnvalidatedGene {
        ensembl_id,
        gene_name: submitted_gene_name,
    } = unvalidated_gene;

    let Some(ensembl_id) = ensembl_id else {
        return Err(ErrorInner::NoEnsemblId);
    };

    let map_valid_gene = |(ensembl_id, gene_name)| ValidGene {
        ensembl_id,
        gene_name,
    };

    if !ensembl_id.is_versionless_and_uppercase() {
        let correct_gene =
            ensembl_id_to_gene(&ensembl_id.to_versionless_uppercase()).map(map_valid_gene);

        return Err(ErrorInner::VersionedOrLowercaseEnsemblId { correct_gene });
    }

    let valid_gene = ensembl_id_to_gene(ensembl_id)
        .map(map_valid_gene)
        .ok_or(ErrorInner::GeneNotFound)?;

    let Some(submitted_gene_name) = submitted_gene_name else {
        return Err(ErrorInner::NoGeneName {
            probable_gene_name: valid_gene.gene_name,
        });
    };

    if *submitted_gene_name == valid_gene.gene_name {
        Ok(valid_gene)
    } else {
        Err(ErrorInner::EnsemblIdGeneNameMismatch {
            correct_gene_name: valid_gene.gene_name,
        })
    }
}

fn parse_bool_from_str(s: Option<&str>, fieldname: &'static str) -> Result<bool, ErrorInner> {
    let Some(s) = s else {
        return Err(ErrorInner::MissingField { fieldname });
    };

    if s.eq_ignore_ascii_case("true") {
        Ok(true)
    } else if s.eq_ignore_ascii_case("false") {
        Ok(false)
    } else {
        Err(ErrorInner::ParseBool {
            value: s.to_owned(),
        })
    }
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct ValidTarget {
    #[serde(flatten)]
    pub gene: ValidGene,
    pub group: String,
    pub is_backup: IsBackup,
    pub must_have: MustHave,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub struct ValidGene {
    pub ensembl_id: EnsemblId,
    pub gene_name: GeneName,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[serde(transparent)]
pub struct IsBackup(pub bool);

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[serde(transparent)]
pub struct MustHave(pub bool);

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedTarget {
    #[serde(flatten)]
    gene: UnvalidatedGene,
    group: Option<String>,
    is_backup: Option<String>,
    must_have: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedGene {
    ensembl_id: Option<UnvalidatedEnsemblId>,
    gene_name: Option<UnvalidatedGeneName>,
}

#[cfg(test)]
mod tests {
    use crate::gene_list::{
        ErrorInner,
        chemistry::{
            UnvalidatedEnsemblId, UnvalidatedGeneName, tests::tp53_ensembl_id,
            xenium_v1_human_ensembl_id_to_gene,
        },
        csv_util::read_csv_trimmed,
        parse_unvalidated_target_from_record,
        target::{
            UnvalidatedGene, UnvalidatedTarget, ValidGene, validate_ensembl_id_gene_name_pair,
            validate_target,
        },
    };

    #[test]
    fn valid_gene() {
        let _ = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap();
    }

    #[test]
    fn unvalidated_target_deserializes_with_invalid_fields() {
        let data = "field1,field2,ensembl_id\nvalue1,value2,id";
        let mut reader = read_csv_trimmed(data);

        let deserialized =
            parse_unvalidated_target_from_record(&reader.records().next().unwrap().unwrap(), None);

        assert_eq!(
            deserialized,
            UnvalidatedTarget {
                gene: UnvalidatedGene {
                    ensembl_id: Some(UnvalidatedEnsemblId::new("id".to_owned())),
                    gene_name: None
                },
                group: None,
                is_backup: None,
                must_have: None
            }
        )
    }

    #[test]
    fn ensembl_id_gene_name_mismatch() {
        let ensembl_id = tp53_ensembl_id();
        let gene_name = UnvalidatedGeneName::new(String::new());

        let err = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id.clone()),
                gene_name: Some(gene_name.clone()),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        let (_correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        assert_eq!(
            err,
            ErrorInner::EnsemblIdGeneNameMismatch { correct_gene_name },
            "failed to create Ensembl ID-gene name mismatch error"
        );
    }

    #[test]
    fn backup_and_must_have_is_rejected() {
        let ensembl_id = tp53_ensembl_id();
        let ensembl_id_str = ensembl_id.as_str();

        let target = UnvalidatedTarget {
            gene: UnvalidatedGene {
                ensembl_id: Some(UnvalidatedEnsemblId::new("TP53".to_owned())),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            group: Some("group0".to_owned()),
            is_backup: Some("true".to_owned()),
            must_have: Some("true".to_owned()),
        };

        let gene_list = format!(
            "ensembl_id,gene_name,group,is_backup,must_have\n{ensembl_id_str},TP53,group0,true,\
             true"
        );

        let errors = validate_target(&target, xenium_v1_human_ensembl_id_to_gene).unwrap_err();

        assert_eq!(errors, [ErrorInner::BackupAndMustHave]);
    }

    #[test]
    fn versioned_or_lowercase_ensembl_id_suggests_correct_gene() {
        let ensembl_id = tp53_ensembl_id();
        let (correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        let versioned =
            UnvalidatedEnsemblId::new(format!("{}.1", ensembl_id.as_str().to_lowercase()));

        let err = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(versioned),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(
            err,
            ErrorInner::VersionedOrLowercaseEnsemblId {
                correct_gene: Some(ValidGene {
                    ensembl_id: correct_ensembl_id,
                    gene_name: correct_gene_name,
                }),
            }
        );
    }

    #[test]
    fn unavailable_ensembl_id_is_gene_not_found() {
        let err = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(UnvalidatedEnsemblId::new("ENSG00000273816".to_owned())),
                gene_name: Some(UnvalidatedGeneName::new(String::new())),
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(err, ErrorInner::GeneNotFound);
    }

    #[test]
    fn missing_gene_name_suggests_probable_name() {
        let ensembl_id = tp53_ensembl_id();
        let (_, correct_gene_name) = xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        let err = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id),
                gene_name: None,
            },
            xenium_v1_human_ensembl_id_to_gene,
        )
        .unwrap_err();

        assert_eq!(
            err,
            ErrorInner::NoGeneName {
                probable_gene_name: correct_gene_name
            }
        );
    }
}
