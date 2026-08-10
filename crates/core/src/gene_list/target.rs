use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::gene_list::{
    ErrorInner,
    chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName},
};

impl UnvalidatedTarget {
    pub fn validate(
        &self,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
    ) -> Result<ValidTarget, Vec<ErrorInner>> {
        let Self {
            gene,
            group,
            priority,
        } = self;

        let mut errors = Vec::new();

        if group.is_none() {
            errors.push(ErrorInner::MissingField { fieldname: "group" });
        }

        let priority = match parse_priority_from_str(priority.as_deref(), "is_backup") {
            Ok(is_backup) => Some(is_backup),
            Err(err) => {
                errors.push(err);
                None
            }
        };

        let valid_gene = match validate_ensembl_id_gene_name_pair(gene, ensembl_id_to_gene) {
            Ok(vg) => Some(vg),
            Err(err) => {
                errors.push(err);
                None
            }
        };

        match (valid_gene, group, priority) {
            (Some(valid_gene), Some(group), Some(priority)) => Ok(ValidTarget {
                gene: valid_gene,
                group: group.to_ascii_lowercase(),
                priority,
            }),
            _ => Err(errors),
        }
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

fn parse_priority_from_str(
    s: Option<&str>,
    fieldname: &'static str,
) -> Result<Priority, ErrorInner> {
    let Some(s) = s else {
        return Err(ErrorInner::MissingField { fieldname });
    };

    Priority::from_str(s).map_err(|_| ErrorInner::ParsePriority {
        value: s.to_owned(),
    })
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct ValidTarget {
    #[serde(flatten)]
    pub gene: ValidGene,
    pub group: String,
    pub priority: Priority,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, strum::EnumString, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum Priority {
    MustHave,
    Desired,
    Backup,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub struct ValidGene {
    pub ensembl_id: EnsemblId,
    pub gene_name: GeneName,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedTarget {
    #[serde(flatten)]
    gene: UnvalidatedGene,
    group: Option<String>,
    priority: Option<String>,
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
                priority: None,
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
