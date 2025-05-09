# ====== 1. Load Date and Libraries ======
print(Sys.Date())

suppressPackageStartupMessages({
  library(harmony)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(Seurat)
  library(tidyverse)
  library(gplots)
  library(cowplot)
  library(Hmisc)
  library(broom)
  library(tidyr)
  library(dplyr)
})

# ====== 2. Parameters ======
sp_rate_threshold <- 0.01

name_list <- c(paste0("HC", 1:8), paste0("Pd", 1:8))
seurat_list <- list()

# ====== 3. Sample Processing Loop ======

for (sample_name in name_list) {
  
  # Read taxonomy info
  bacteria_info <- read.csv(paste0(species_path, sample_name, "_sc_taxonomy_S.report"), sep = "\t", header = TRUE)
  bacteria_info_G <- read.csv(paste0(species_path, sample_name, "_sc_taxonomy_G.report"), sep = "\t", header = TRUE)
  
  # Load and preprocess raw data
  sample <- Read10X(paste0(data_path, '/', sample_name))
  sample <- rename_genes(sample, map_df)
  sample <- CreateSeuratObject(counts = sample, min.cells = 3)
  
  # MAD-based QC filtering
  mad_nCount <- mad(sample$nCount_RNA)
  mad_nFeature <- mad(sample$nFeature_RNA)
  median_nCount <- median(sample$nCount_RNA)
  median_nFeature <- median(sample$nFeature_RNA)

  # Filter out outliers 
  sample <- subset(sample, subset = (
    nCount_RNA > (median_nCount - 3 * mad_nCount) &
      nCount_RNA < (median_nCount + 3 * mad_nCount) &
      nFeature_RNA > (median_nFeature - 3 * mad_nFeature) &
      nFeature_RNA < (median_nFeature + 3 * mad_nFeature)
  ))
  
  # Merge species annotations
  BC <- colnames(sample) %>% as.data.frame()
  colnames(bacteria_info) <- c("BC", "index", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads", "all_reads")
  BC_info <- merge(BC, bacteria_info, all.x = TRUE)
  sample$species_info <- BC_info[, 2]
  
  # Remove cells without species info
  valid_cells <- !is.na(sample$species_info)
  sample <- sample[, valid_cells]
  
  # Filter low-frequency species
  species_counts <- table(sample$species_info)
  num_threshold <- ncol(sample) * sp_rate_threshold
  abundant_species <- names(species_counts[species_counts >= num_threshold])
  sample <- subset(sample, subset = species_info %in% abundant_species)
  
  # Remove artificial SEQ genes
  sample <- subset(sample, features = rownames(sample)[!grepl("^SEQ", rownames(sample))])
  
  # Drop empty cells
  sample <- sample[, colSums(GetAssayData(sample)) != 0]
  
  # Rename cells with sample prefix
  sample <- RenameCells(sample, new.names = paste(sample_name, colnames(sample), sep = "_"))
  
  # Add sample info
  sample$orig.ident <- sample_name
  set.seed(123)
  
  # Save into list
  seurat_list[[sample_name]] <- sample
  
  # Print progress
  print(sample_name)
  print(dim(sample))
}

# ====== 4. Combine All Samples ======
combined_seurat_object <- if (length(seurat_list) > 1) {
  Reduce(function(x, y) merge(x, y), seurat_list)
} else {
  seurat_list[[1]]
}

# ====== 5. Add Metadata and Normalize ======
sample_to_group <- setNames(rep(c("HC", "Pd"), each = 8), name_list)
combined_seurat_object$group <- sapply(combined_seurat_object$orig.ident, function(x) sample_to_group[[x]])

# SCTransform & Harmony Integration
combined_seurat_object <- SCTransform(
  combined_seurat_object, 
  vars.to.regress = NULL, 
  verbose = FALSE, 
  return.only.var.genes = FALSE, 
  variable.features.n = 2000
)
combined_seurat_object <- RunPCA(combined_seurat_object)
combined_seurat_object <- RunHarmony(
  combined_seurat_object, 
  reduction = "pca", 
  group.by.vars = "orig.ident", 
  reduction.save = "harmony"
)

# ====== 6. UMAP and Clustering ======
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.4)

combined_seurat_object$BC <- rownames(combined_seurat_object@meta.data)
combined_seurat_object$orig.ident <- factor(combined_seurat_object$orig.ident, levels = name_list)

# ====== 7. Save Seurat Object ======
saveRDS(combined_seurat_object, "Combined_alter.rds")
