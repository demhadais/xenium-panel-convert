# Single-cell RNA sequencing Fixed RNA Capture (Flex) Probe-sets

[The build script](../build.rs) reads the probe-set lists in this directory to construct the expected set of genes, so the package can validate whether a user has filtered out genes from their anndata (H5AD). The probe-sets are versioned, so a user's data might have been generated with any of the following probe-sets:

- **Human Transcriptome v2.0.0 (GRCh38-2024-A)**: https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv

- **Mouse Transcriptome v2.0.0 (GRCm39-2024-A)**: https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv

- **Human Transcriptome v1.1.0 (GRCh38-2024-A)**: https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv

- **Mouse Transcriptome v1.1.1 (GRCm39-2024-A)**: https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.csv

- **Human Transcriptome v1.0.1 (GRCh38-2020-A)**: https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv

- **Mouse Transcriptome v1.0.1 (mm10-2020-A)**: https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv

A listing of the probe-sets can be found on [10x Genomics' downloads page](https://www.10xgenomics.com/support/software/cell-ranger/downloads#probe-set-downloads), and more information about them can be found in [their documentation](https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets?tag%5BstepSlugs%5D=probe-sets&tag%5BproductSlugs%5D=flex-gene-expression&tag%5BisAssayDocumentStr%5D=true).
