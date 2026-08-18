use std::{convert::Infallible, fmt::Display, str::FromStr};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellBarcodeCol(String);

impl From<&str> for CellBarcodeCol {
    fn from(value: &str) -> Self {
        Self(value.into())
    }
}

impl Default for CellBarcodeCol {
    fn default() -> Self {
        Self(String::from("_index"))
    }
}

impl Display for CellBarcodeCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for CellBarcodeCol {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        String::from_str(s).map(Self)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellAnnotationCol(String);

impl From<&str> for CellAnnotationCol {
    fn from(value: &str) -> Self {
        Self(value.into())
    }
}

impl Display for CellAnnotationCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for CellAnnotationCol {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        String::from_str(s).map(Self)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EnsemblIdCol(String);

impl From<&str> for EnsemblIdCol {
    fn from(value: &str) -> Self {
        Self(value.into())
    }
}

impl EnsemblIdCol {
    #[must_use]
    pub(super) fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for EnsemblIdCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Default for EnsemblIdCol {
    fn default() -> Self {
        Self(String::from("gene_ids"))
    }
}

impl FromStr for EnsemblIdCol {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        String::from_str(s).map(Self)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneNameCol(String);

impl From<&str> for GeneNameCol {
    fn from(value: &str) -> Self {
        Self(value.into())
    }
}

impl GeneNameCol {
    #[must_use]
    pub(super) fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for GeneNameCol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Default for GeneNameCol {
    fn default() -> Self {
        Self(String::from("_index"))
    }
}

impl FromStr for GeneNameCol {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        String::from_str(s).map(Self)
    }
}
