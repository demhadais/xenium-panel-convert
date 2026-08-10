use serde::Serialize;

use crate::gene_list::{
    chemistry::GeneName,
    target::{UnvalidatedTarget, ValidGene},
};

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Error {
    pub(super) line_number: Option<u64>,
    pub(super) submitted_target: Option<UnvalidatedTarget>,
    pub(super) errors: Vec<ErrorInner>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum ErrorInner {
    MalformedCsv {
        reason: String,
    },
    MissingField {
        fieldname: &'static str,
    },
    ParseBool {
        value: String,
    },
    VersionedOrLowercaseEnsemblId {
        correct_gene: Option<ValidGene>,
    },
    NoEnsemblId,
    NoGeneName {
        probable_gene_name: GeneName,
    },
    RenamedField {
        original_fieldname: String,
        correct_fieldname: String,
    },
    EnsemblIdGeneNameMismatch {
        correct_gene_name: GeneName,
    },
    BackupAndMustHave,
    GeneNotFound,
    DuplicateGene,
}

impl From<csv::Error> for ErrorInner {
    fn from(err: csv::Error) -> Self {
        Self::MalformedCsv {
            reason: err.to_string(),
        }
    }
}

impl<'a> From<&'a csv::Error> for ErrorInner {
    fn from(err: &'a csv::Error) -> Self {
        Self::MalformedCsv {
            reason: err.to_string(),
        }
    }
}
