use std::str::FromStr;

use csv::StringRecord;
use serde::{Deserialize, Serialize};

use crate::gene_list::{
    ErrorInner,
    chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName},
};

impl UnvalidatedTarget {
    pub(super) fn from_record(record: &StringRecord, fieldnames: &StringRecord) -> Self {
        // Unwrapping is fine because extra fields won't cause a failure, nor will
        // missing fields
        record.deserialize(Some(fieldnames)).unwrap()
    }

    pub(super) fn validate(
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

        let priority = match parse_priority_field(priority.as_deref()) {
            Ok(priority) => Some(priority),
            Err(err) => {
                errors.push(err);
                None
            }
        };

        let valid_gene = match ValidGene::from_unvalidated(gene, ensembl_id_to_gene) {
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

fn parse_priority_field(s: Option<&str>) -> Result<Priority, ErrorInner> {
    let Some(s) = s else {
        return Err(ErrorInner::MissingField {
            fieldname: "priority",
        });
    };

    Priority::from_str(s).map_err(|_| ErrorInner::ParsePriority {
        value: s.to_owned(),
    })
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub(crate) struct ValidTarget {
    #[serde(flatten)]
    pub(crate) gene: ValidGene,
    pub(crate) group: String,
    pub(crate) priority: Priority,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, strum::EnumString, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub(crate) enum Priority {
    MustHave,
    Desired,
    Backup,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub(crate) struct ValidGene {
    pub(crate) ensembl_id: EnsemblId,
    pub(crate) gene_name: GeneName,
}

impl ValidGene {
    fn from_unvalidated(
        UnvalidatedGene {
            ensembl_id,
            gene_name: submitted_gene_name,
        }: &UnvalidatedGene,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
    ) -> Result<Self, ErrorInner> {
        let Some(ensembl_id) = ensembl_id else {
            return Err(ErrorInner::NoEnsemblId);
        };

        let map_valid_gene = |(ensembl_id, gene_name)| Self {
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
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub(crate) struct UnvalidatedTarget {
    #[serde(flatten)]
    gene: UnvalidatedGene,
    group: Option<String>,
    priority: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
struct UnvalidatedGene {
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
        target::{UnvalidatedGene, UnvalidatedTarget, ValidGene},
    };

    #[test]
    fn valid_gene() {
        ValidGene::from_unvalidated(
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

        let fieldnames = reader.headers().unwrap().clone();
        let record = reader.records().next().unwrap().unwrap();
        let deserialized = UnvalidatedTarget::from_record(&record, &fieldnames);

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
        );
    }

    #[test]
    fn ensembl_id_gene_name_mismatch() {
        let ensembl_id = tp53_ensembl_id();
        let gene_name = UnvalidatedGeneName::new(String::new());

        let err = ValidGene::from_unvalidated(
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

        let err = ValidGene::from_unvalidated(
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
        let err = ValidGene::from_unvalidated(
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

        let err = ValidGene::from_unvalidated(
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
