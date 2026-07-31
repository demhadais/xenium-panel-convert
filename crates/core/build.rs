use std::str::FromStr;

use crate::xenium_panel_allowed_genes::write_gene_maps;

#[path = "build/complete_feature_sets.rs"]
mod complete_feature_sets;
#[path = "build/xenium_panel_allowed_genes.rs"]
mod xenium_panel_allowed_genes;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    if !option_env!("BUILD_XENIUM_PANEL_VALIDATE")
        .map(bool::from_str)
        .transpose()?
        .is_some_and(|build| build)
    {
        return Ok(());
    }

    write_gene_maps().await?;

    Ok(())
}
