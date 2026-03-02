
# Microbial single-cell RNA sequencing for subgingival microbiome

This repository contains metadata and analysis scripts for the single-cell RNA sequencing study of the subgingival microbiome in health control and periodontitis, as presented in [Ding et al.](https://www.cell.com/cell-host-microbe/abstract/S1931-3128(26)00042-9) (*Cell Host & Microbe*, 2026).

For questions regarding this dataset or analysis, please contact:

- Yongcheng Wang (yongcheng@zju.edu.cn)
- Pei-Hui Ding (phding@zju.edu.cn)

## Data Availability

### Repositary structure

All data files referenced in the analysis scripts (except RDS files) are located in the `data/` directory.

### Data Access

- **Processed data**: RDS files, expression matrices, and taxonomy files are available via [Figshare](https://figshare.com/s/867dbfdd7997e78d6e71).

- **Sequencing data** ([Genome Sequence Archive](https://ngdc.cncb.ac.cn/gsa/)):
  - scRNA-seq: CRA023798
  - *N.elongata* whole-genome sequencing: CRA026826
  - 16S rRNA sequencing: CRA030816

All sequencing data will be made publicly available upon final publication of the article.

### Analysis Scripts

**Included scripts**:

`figure2.R`: Workflow of QC; clustering; box plot of gene counts; correlation analysis; comparisons between HC and Pd samples.

`figure3.R`: Clustering; WGCNA; marker gene identification.

`figure4.R`: Transcriptional alterations in subgingival bacteria under periodontitis. 

`figure5.R`: Functional modules and heterogeneity of *N. elongata* in HC and Pd groups.

`figure6.R`: Metabolic signatures of *P. gingivalis*, *P. intermedia*, and *T. denticola* in periodontitis.



## Requirement

R packages:

- Seurat (v4.4.0)
- SeuratObject (v5.0.2)
- ggplot2
- data.table
- Hmisc
- monocle
- hdWGCNA
- SingleCellExperiment
- DESeq2
- harmony



