# ====== Load Required Packages ======
library(harmony)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyr)
library(muscat)
library(ggsignif)

# ====== Set Working Directory and Load Data ======
setwd('D:/AA_R_analysis/a_subgingival/')
Combined_seurat_object <- readRDS('Combined_alter.rds')

# ====== Prepare COG Functional Annotation Table ======
cog.20.def <- read.delim("cog-20.def.tab", header = FALSE)
colnames(cog.20.def) = c('gene','COG_class','Func_anno','name','Func','time','Anno')
# Ensure 'name' field is populated
cog.20.def <- cog.20.def %>% mutate(name = ifelse(name == '', gene, name))

# ====== Convert to SingleCellExperiment and Format Metadata ======
DefaultAssay(Combined_seurat_object) <- "RNA"
Combined_seurat_object_sce <- as.SingleCellExperiment(Combined_seurat_object)
Combined_seurat_object_sce <- prepSCE(
  Combined_seurat_object_sce,
  kid = "celltype",      # Cluster info
  gid = "group",         # Group: HC vs Pd
  sid = "orig.ident",    # Sample ID
  drop = TRUE
)
Combined_seurat_object_sce$group_id <- factor(Combined_seurat_object_sce$group_id, levels = c("HC", "Pd"))

# ====== Pseudo-bulk Aggregation by Cell Type and Sample ======
pb <- aggregateData(
  Combined_seurat_object_sce,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)

# ====== Run Differential Expression Analysis using muscat ======
res <- muscat::pbDS(pb)

# ====== Extract Results for Pd vs HC Comparison ======
test <- res[["table"]][["Pd"]]

# ====== Combine and Annotate DEGs ======
combined_DEG <- do.call(rbind, res[["table"]][["Pd"]]) %>%
  filter(p_val < 0.05) %>%
  rename(
    avg_log2FC = logFC,
    cluster = cluster_id,
    p_val_adj = p_adj.loc,
    name = gene
  ) %>%
  select(cluster, name, avg_log2FC, p_val, p_val_adj) %>%
  arrange(cluster, avg_log2FC) %>%
  left_join(cog.20.def, by = 'name', relationship = "many-to-many")

# ====== Save DEG Table ======
write.csv(combined_DEG, 'DEG_bulk.csv', row.names = FALSE)