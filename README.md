
# Microbial single-cell RNA sequencing for subgingival microbiome

### 1. Data Preprocessing

`01_data_loading_QC.R`

- Performs data loading and quality control filtering
- Applies batch correction (Harmony)
- Conducts cell clustering using Seurat

### 2. Functional Annotation

`02_annotation.R`

- Performs co-expression network analysis using hdWGCNA
- Visualizes key marker gene expression patterns
- Assigns functional annotations to cell clusters

### 3. *N. elongata* Subpopulation Analysis

`03_N.elongata_standard_analysis.R`

- Isolates *N. elongata* subpopulation
- Performs sub-clustering and cluster-specific DEG analysis
- Annotates functions for each subcluster

### 4. *N. elongata* Trajectory Inference

`04_N.elongata_monocle.R`

- Consctructs single-cell trajectories using Monocle2
- Identifies pseudotime-dependent gene expression changes

### 5. Cross-Species Differential Analysis

`05_pseudobulk_DEG_species.R`

- Performs pseudobulk-based differential expression analysis
- Compares *T. denticola* vs *P. gingivalis* transcriptomes

### 5. Condition-Specific Differential Analysis

`06_pseudobulk_DEG_groups.R`

- Conducts pseudobulk DEG analysis between:
  - Healty control (HC) group
  - Periodontitis (Pd) group 




