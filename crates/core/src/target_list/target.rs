use std::str::FromStr;

use csv::StringRecord;
use serde::{Deserialize, Serialize};
use strum::VariantNames;

use crate::{
    common::ErrorVecExt,
    target_list::{
        TargetErrorInner,
        chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName},
    },
};

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedGene {
    pub ensembl_id: Option<UnvalidatedEnsemblId>,
    pub gene_name: Option<UnvalidatedGeneName>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct UnvalidatedTarget {
    #[serde(flatten)]
    pub gene: UnvalidatedGene,
    pub group: Option<String>,
    pub priority: Option<String>,
}

impl UnvalidatedTarget {
    pub(super) fn from_record(record: &StringRecord, fieldnames: &StringRecord) -> Self {
        // Unwrapping is fine because extra fields won't cause a failure, nor will
        // missing fields
        record.deserialize(Some(fieldnames)).unwrap()
    }

    pub(super) fn validate(
        &self,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
    ) -> Result<ValidTarget, Vec<TargetErrorInner>> {
        let Self {
            gene,
            group,
            priority,
        } = self;

        let mut errors = Vec::new();

        if group.is_none() {
            errors.push(TargetErrorInner::MissingField { fieldname: "group" });
        }

        let priority =
            parse_priority_field(priority.as_deref()).map_or_else(|err| errors.push_err(err), Some);

        let valid_gene = ValidGene::from_unvalidated(gene, ensembl_id_to_gene)
            .map_or_else(|err| errors.push_err(err), Some);

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

fn parse_priority_field(s: Option<&str>) -> Result<Priority, TargetErrorInner> {
    let Some(s) = s else {
        return Err(TargetErrorInner::MissingField {
            fieldname: "priority",
        });
    };

    Priority::from_str(s).map_err(|_| TargetErrorInner::InvalidPriority {
        value: s.to_owned(),
        allowed: Priority::VARIANTS,
    })
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub struct ValidGene {
    pub ensembl_id: EnsemblId,
    pub gene_name: GeneName,
}

impl ValidGene {
    fn from_unvalidated(
        UnvalidatedGene {
            ensembl_id,
            gene_name: submitted_gene_name,
        }: &UnvalidatedGene,
        ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
    ) -> Result<Self, TargetErrorInner> {
        let Some(ensembl_id) = ensembl_id else {
            return Err(TargetErrorInner::NoEnsemblId);
        };

        let map_valid_gene = |(ensembl_id, gene_name)| Self {
            ensembl_id,
            gene_name,
        };

        if !ensembl_id.is_versionless_and_uppercase() {
            let correct_gene =
                ensembl_id_to_gene(&ensembl_id.to_versionless_uppercase()).map(map_valid_gene);

            return Err(TargetErrorInner::VersionedOrLowercaseEnsemblId { correct_gene });
        }

        let valid_gene = ensembl_id_to_gene(ensembl_id)
            .map(map_valid_gene)
            .ok_or(TargetErrorInner::GeneNotFound)?;

        let Some(submitted_gene_name) = submitted_gene_name else {
            return Err(TargetErrorInner::NoGeneName {
                probable_gene_name: valid_gene.gene_name,
            });
        };

        if *submitted_gene_name == valid_gene.gene_name {
            Ok(valid_gene)
        } else {
            Err(TargetErrorInner::EnsemblIdGeneNameMismatch {
                ensembl_id: valid_gene.ensembl_id,
                correct_gene_name: valid_gene.gene_name,
            })
        }
    }
}

#[derive(
    Clone,
    Copy,
    Debug,
    Serialize,
    PartialEq,
    Eq,
    strum::EnumString,
    PartialOrd,
    Ord,
    strum::VariantNames,
)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub(crate) enum Priority {
    MustHave,
    Desired,
    Backup,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct ValidTarget {
    #[serde(flatten)]
    pub(crate) gene: ValidGene,
    pub(crate) group: String,
    pub(crate) priority: Priority,
}

#[cfg(test)]
mod tests {
    use crate::target_list::{
        TargetErrorInner,
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

        let (correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene(&ensembl_id).unwrap();

        assert_eq!(
            err,
            TargetErrorInner::EnsemblIdGeneNameMismatch {
                ensembl_id: correct_ensembl_id,
                correct_gene_name
            },
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
            TargetErrorInner::VersionedOrLowercaseEnsemblId {
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

        assert_eq!(err, TargetErrorInner::GeneNotFound);
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
            TargetErrorInner::NoGeneName {
                probable_gene_name: correct_gene_name
            }
        );
    }
}
