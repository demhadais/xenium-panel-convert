use serde::Serialize;

use crate::{
    common::ErrorVecExt,
    target_list::{
        chemistry::GeneName,
        target::{UnvalidatedTarget, ValidGene},
    },
};

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct TargetErrorSet {
    pub(super) line_number: Option<u64>,
    pub(super) submitted_target: Option<UnvalidatedTarget>,
    pub(super) errors: Vec<TargetErrorInner>,
}

#[derive(Clone, Debug, Serialize, PartialEq, thiserror::Error)]
#[serde(rename_all = "snake_case", tag = "type")]
pub(crate) enum TargetErrorInner {
    #[error("malformed CSV - {reason}")]
    MalformedCsv {
        reason: String,
    },
    #[error("missing field - {fieldname}")]
    MissingField {
        fieldname: &'static str,
    },
    #[error("invalid priority - {value}")]
    InvalidPriority {
        value: String,
    },

    #[error(
        "versioned or lowercase Ensembl ID{}",
        correct_gene.as_ref().map_or_else(
            String::new,
            |gene| format!(" - the correct gene is {} ({})", gene.gene_name, gene.ensembl_id),
        )
    )]
    VersionedOrLowercaseEnsemblId {
        correct_gene: Option<ValidGene>,
    },
    #[error("no Ensembl ID")]
    NoEnsemblId,
    #[error("no gene name - the gene name is probably {probable_gene_name}")]
    NoGeneName {
        probable_gene_name: GeneName,
    },
    #[error("renamed field - {original_fieldname} should be {correct_fieldname}")]
    RenamedField {
        original_fieldname: String,
        correct_fieldname: String,
    },
    #[error("Ensembl ID and gene name mismatch - the correct gene name is {correct_gene_name}")]
    EnsemblIdGeneNameMismatch {
        correct_gene_name: GeneName,
    },
    #[error("gene not found")]
    GeneNotFound,
    #[error("duplicate gene")]
    DuplicateGene,
}

impl From<csv::Error> for TargetErrorInner {
    fn from(err: csv::Error) -> Self {
        Self::MalformedCsv {
            reason: err.to_string(),
        }
    }
}

impl<'a> From<&'a csv::Error> for TargetErrorInner {
    fn from(err: &'a csv::Error) -> Self {
        Self::MalformedCsv {
            reason: err.to_string(),
        }
    }
}

impl ErrorVecExt<TargetErrorInner> for Vec<TargetErrorInner> {
    fn push_err<T>(&mut self, err: TargetErrorInner) -> Option<T> {
        self.push(err);

        None
    }
}
