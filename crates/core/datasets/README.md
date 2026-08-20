# 3' Gene Expression (polyA-capture) Datasets

[The build script](../build.rs) reads the datasets in this directory to extract transcriptomes so the package can validate whether a user has filtered out genes from their anndata (H5AD). Each dataset corresponds to a given reference transcriptome:

- **GRCh38-2020-A**: https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5

- **mm10-2020-A**: https://cf.10xgenomics.com/samples/cell-exp/7.0.0/1k_mouse_kidney_CNIK_3pv3/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.h5

- **GRCh38-2024-A**: https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5

- **GRCm39-2024-A**: https://cf.10xgenomics.com/samples/cell-exp/10.0.0/SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1/SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5

More information about the reference transcriptomes is available in [10x Genomics' documentation](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads).
