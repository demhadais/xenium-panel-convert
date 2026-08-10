use std::{
    fs,
    io::{self, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use bytes::Bytes;
use serde::Deserialize;
use url::Url;

pub(crate) async fn write_gene_maps() -> anyhow::Result<()> {
    let Config {
        xenium_v1_gene_list,
        xenium_prime_gene_list,
    } = toml::from_slice(include_bytes!("../genes.toml"))
        .context("failed to parse config from genes.toml")?;

    let http_client = reqwest::Client::new();
    let raw_gene_lists = [xenium_v1_gene_list, xenium_prime_gene_list]
        .map(|url| fetch_raw_gene_list(&http_client, url));

    let raw_gene_lists = futures::future::try_join_all(raw_gene_lists).await?;
    let mut gene_lists: Vec<_> = raw_gene_lists
        .iter()
        .map(|raw| csv::Reader::from_reader(raw.as_ref()))
        .collect();

    let gene_lists: Vec<Vec<_>> = gene_lists
        .iter_mut()
        .map(|list| list.deserialize().map(|res| res.unwrap()))
        .map(Iterator::collect)
        .collect();

    let gene_lists: Vec<_> = gene_lists.iter().map(|list| construct_maps(list)).collect();
    let [
        GeneMaps {
            homo_sapiens: v1_human,
            mus_musculus: v1_mouse,
        },
        GeneMaps {
            homo_sapiens: prime_human,
            mus_musculus: prime_mouse,
        },
    ] = gene_lists.as_array().unwrap();

    for (filename, map_name, gene_map) in [
        ("xenium_v1_human.rs", "XENIUM_V1_HUMAN_GENES", v1_human),
        ("xenium_v1_mouse.rs", "XENIUM_V1_MOUSE_GENES", v1_mouse),
        (
            "xenium_prime_human.rs",
            "XENIUM_PRIME_HUMAN_GENES",
            prime_human,
        ),
        (
            "xenium_prime_mouse.rs",
            "XENIUM_PRIME_MOUSE_GENES",
            prime_mouse,
        ),
    ] {
        write_map_to_file(
            &PathBuf::from(format!("src/gene_list/chemistry/{filename}")),
            map_name,
            gene_map,
        )?;
    }

    Ok(())
}

#[derive(Clone, Debug, Deserialize)]
struct Config {
    xenium_v1_gene_list: Url,
    xenium_prime_gene_list: Url,
}

fn construct_maps(gene_list: &[Gene]) -> GeneMaps<'_> {
    fn insert_gene<'a>(
        ensembl_id: &'a str,
        gene_symbol: &'a str,
        map: &mut phf_codegen::Map<'a, &'a str>,
    ) {
        map.entry(ensembl_id, format!(r#""{gene_symbol}""#));
    }

    let mut homo_sapiens = phf_codegen::Map::new();
    let mut mus_musculus = phf_codegen::Map::new();

    for Gene {
        species,
        ensembl_id,
        gene_symbol,
    } in gene_list
    {
        match species.as_str() {
            "Homo sapiens" => insert_gene(ensembl_id, gene_symbol, &mut homo_sapiens),
            "Mus musculus" => insert_gene(ensembl_id, gene_symbol, &mut mus_musculus),
            s => panic!("species {s} not expected"),
        }
    }

    GeneMaps {
        homo_sapiens,
        mus_musculus,
    }
}

struct GeneMaps<'a> {
    homo_sapiens: phf_codegen::Map<'a, &'a str>,
    mus_musculus: phf_codegen::Map<'a, &'a str>,
}

fn write_map_to_file(
    path: &Path,
    map_name: &str,
    map: &phf_codegen::Map<'_, &str>,
) -> anyhow::Result<()> {
    let file = fs::File::create(path)
        .with_context(|| format!("failed to write file {}", path.to_str().unwrap()))?;
    let mut file_writer = io::BufWriter::new(file);

    writeln!(
        file_writer,
        "pub(super) static {map_name}: phf::Map<&'static str, &'static str> = {};",
        map.build()
    )?;

    Ok(())
}

async fn fetch_raw_gene_list(client: &reqwest::Client, url: Url) -> anyhow::Result<Bytes> {
    let response = client.get(url).send().await?;

    let raw = response.bytes().await?;

    Ok(raw)
}

#[derive(Deserialize)]
struct Gene {
    #[serde(rename = "Species")]
    species: String,
    #[serde(rename = "Ensembl ID")]
    ensembl_id: String,
    #[serde(rename = "Gene symbol")]
    gene_symbol: String,
}
