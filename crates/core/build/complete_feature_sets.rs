use std::path::Path;

use hdf5_metno::types::FixedAscii;

use crate::common::write_map_to_file;

pub(crate) fn write_complete_feature_sets() -> anyhow::Result<()> {
    let grch38 = (
        "datasets/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5",
        "GRCH38_2020_A",
        "grch38_2020_a.rs",
    );
    let grch38_flex = (
        "datasets/4plex_human_colorectal_kidney_scFFPE_multiplex_Kidney_Manual_BC3_count_sample_filtered_feature_bc_matrix.h5",
        "GRCH38_2020_A_FLEX",
        "grch38_2020_a_flex.rs",
    );
    let mm10 = (
        "datasets/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5",
        "MM10",
        "mm10_2020_a.rs",
    );
    let mm10_flex = (
        "datasets/10k_mouse_spleen_scFFPE_singleplex_10k_mouse_spleen_scFFPE_singleplex_count_sample_filtered_feature_bc_matrix.h5",
        "MM10_FLEX",
        "mm10_2020_a_flex.rs",
    );

    for (h5_path, map_name, map_path) in [grch38, grch38_flex, mm10, mm10_flex] {
        let (ensembl_ids, gene_names) = read_feature_set(h5_path)?;

        let map = construct_map(&ensembl_ids, &gene_names);

        let map_path = format!("src/reference_dataset/feature_set/{map_path}");
        write_map_to_file(Path::new(&map_path), map_name, &map)?;
    }

    Ok(())
}

// Human Ensembl IDs are 15 characters while mouse Ensembl IDs are 18
type EnsemblId = FixedAscii<18>;

// No gene name is likely to exceed 32 characters
type GeneName = FixedAscii<32>;

fn read_feature_set(file_path: &str) -> anyhow::Result<(Vec<EnsemblId>, Vec<GeneName>)> {
    let file = hdf5_metno::File::open(file_path)?;

    let ensembl_ids = file
        .dataset("matrix/features/id")
        .and_then(|ds| ds.read_raw())?;
    let gene_names = file
        .dataset("matrix/features/name")
        .and_then(|ds| ds.read_raw())?;

    Ok((ensembl_ids, gene_names))
}

fn construct_map<'a>(
    ensembl_ids: &'a [EnsemblId],
    gene_names: &'a [GeneName],
) -> phf_codegen::Map<'a, &'a str> {
    let mut map = phf_codegen::Map::new();

    for (id, name) in ensembl_ids.iter().zip(gene_names) {
        map.entry(id.as_str(), format!(r#""{name}""#));
    }

    map
}
