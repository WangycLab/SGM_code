# Load required libraries
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(paletteer)
library(plot1cell)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
library(ggcor)
library(ggridges)
library(gmp)
library(ClusterGVis)
library(viridis)
library(gridExtra)
library(monocle)
library(ggpubr)
library(scales)
library(ggnewscale)

# Load Seurat object and set cell identity
Species_sub_mono <- readRDS('./Ne.rds')
Idents(Species_sub_mono) <- 'celltype_sub'

# Add group information to metadata
Species_sub_mono@meta.data <- Species_sub_mono@meta.data %>%
  mutate(group = case_when(
    grepl('HC', orig.ident) ~ 'HC',
    grepl('Pd', orig.ident) ~ 'Pd'
  ))

# Prepare Monocle input
expression_matrix <- as(as.matrix(Species_sub_mono@assays$SCT@counts), 'sparseMatrix')
cell_metadata <- new('AnnotatedDataFrame', data = Species_sub_mono@meta.data)
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
feature_data <- new('AnnotatedDataFrame', data = gene_annotation)

# Create Monocle CellDataSet
cds <- newCellDataSet(expression_matrix,
                      phenoData = cell_metadata,
                      featureData = feature_data,
                      expressionFamily = negbinomial.size())

# Preprocessing
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# Filter expressed genes
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# Differential expression test across sub-celltypes
diff <- differentialGeneTest(cds[expressed_genes,],
                             fullModelFormulaStr = "~celltype_sub",
                             cores = 1)

# Select ordering genes
ordering_genes <- rownames(subset(diff, qval < 0.01))

# Dimensionality reduction and trajectory inference
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 1)

# Extract pseudotime and component data
data_df <- t(reducedDimS(cds)) %>%
  as.data.frame() %>%
  dplyr::select(`Component_1` = 1, `Component_2` = 2) %>%
  rownames_to_column("cells") %>%
  left_join(pData(cds) %>%
              dplyr::select(BC, State, Pseudotime, orig.ident, celltype_sub),
            by = c('cells' = 'BC')) %>%
  mutate(group = case_when(
    grepl('HC', orig.ident) ~ 'HC',
    grepl('Pd', orig.ident) ~ 'Pd'
  ))

# Plot trajectory
plot_cell_trajectory(cds)

# BEAM test at branch point 1
BEAM_res <- BEAM(cds[ordering_genes,],
                 branch_point = 1,
                 cores = 2,
                 progenitor_method = 'duplicate') %>%
  arrange(qval) %>%
  dplyr::select(gene_short_name, pval, qval)

# Branched heatmap visualization
heatmap_df <- plot_genes_branched_heatmap2(
  cds[row.names(subset(BEAM_res, pval < 0.05)),],
  branch_point = 1,
  num_clusters = 3,
  cores = 3,
  use_gene_short_name = FALSE,
  show_rownames = FALSE
)

# Save branched heatmap visualization
pdf('Ne_Monocle_BeamHeatmap.pdf', height = 10, width = 8)
visCluster(
  object = heatmap_df,
  pseudotime_col = c("#9E0142", "#377EB8", "#D53E4F"),
  plot.type = "both",
  ctAnno.col = c("#9E0142", "#377EB8", "#D53E4F"),
  cluster.order = c(1, 2, 3)
)
dev.off()

# Extract cluster info from heatmap
heatmap <- plot_genes_branched_heatmap(
  cds[row.names(subset(BEAM_res, pval < 0.05)),],
  branch_point = 1,
  num_clusters = 3,
  cores = 3,
  use_gene_short_name = FALSE,
  return_heatmap = TRUE,
  show_rownames = FALSE
)
clusters <- cutree(heatmap$ph[["tree_row"]], k = 3)
clustering <- data.frame(clusters = as.character(clusters),
                         gene = rownames(heatmap$ph[["tree_row"]]))

# Merge results with annotation (assumes cog.20.def is loaded)
Time_genes <- BEAM_res %>%
  filter(pval < 0.05) %>%
  arrange(desc(qval)) %>%
  pull(gene_short_name)

branched_gene_df <- clustering %>%
  left_join(cog.20.def, by = 'gene') %>%
  left_join(BEAM_res, by = c('gene' = 'gene_short_name')) %>%
  arrange(clusters, qval)

# Save final gene table
write.csv(branched_gene_df, 'Ne_diffgene_state.csv', row.names = FALSE)