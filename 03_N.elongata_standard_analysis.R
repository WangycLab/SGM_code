#============================#
#     Required Libraries     #
#============================#
library(Seurat)
library(harmony)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
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

#============================#
#        Load Data           #
#============================#
combined_seurat <- readRDS('./Combined_alter.rds')
species_id <- 'Neisseria elongata'
species_label <- 'N.elongata'

species_sub <- subset(combined_seurat, subset = species_info %in% species_id)

#===============================#
#   Preprocessing & Clustering  #
#===============================#
species_sub <- SCTransform(species_sub, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = 2000)
species_sub <- RunPCA(species_sub)
species_sub <- RunUMAP(species_sub, dims = 1:10)
species_sub <- FindNeighbors(species_sub, dims = 1:10)
species_sub <- FindClusters(species_sub, resolution = 0.35)
species_sub$celltype_sub <- species_sub$seurat_clusters

#============================#
#        UMAP Plot           #
#============================#
DimPlot(species_sub, reduction = "umap", label = TRUE, group.by = 'seurat_clusters', split.by = 'group',
        pt.size = 0.5, cols = cluster_color_sub) +
  labs(x = "u-MAP-Dim1", y = "u-MAP-Dim2") + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", size = 0.4)
  )

#============================#
#     Find Cluster Markers   #
#============================#
cluster_markers <- FindAllMarkers(species_sub, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.table(cluster_markers, file = "Ne.cluster.markers.xls", sep = "\t", row.names = FALSE, quote = TRUE)

#============================#
#   Define Functional Types  #
#============================#
functional_order <- c(
  "Twitching motility", "Translation regulation", "Nicotinamide metabolism",
  "Cell adhesion", "Signal transduction", "Binding-effector",
  "Polysaccharides degradation", "Antioxidant defense"
)
names(cluster_color_sub) <- functional_order

species_sub <- species_sub %>%
  mutate(
    celltype_sub = case_when(
      seurat_clusters == 0 ~ "Twitching motility",
      seurat_clusters == 1 ~ "Translation regulation",
      seurat_clusters == 2 ~ "Nicotinamide metabolism",
      seurat_clusters == 3 ~ "Cell adhesion",
      seurat_clusters == 4 ~ "Signal transduction",
      seurat_clusters == 5 ~ "Binding-effector",
      seurat_clusters == 6 ~ "Polysaccharides degradation",
      seurat_clusters == 7 ~ "Antioxidant defense"
    ),
    celltype_sub = factor(celltype_sub, levels = functional_order)
  )

#============================#
#         Plot UMAP          #
#============================#
umap_df <- Embeddings(species_sub, "umap") %>%
  as.data.frame() %>%
  cbind(
    cellType = species_sub$celltype_sub,
    sample = species_sub$orig.ident,
    group = species_sub$group
  ) %>%
  mutate(cellType = factor(cellType, levels = functional_order))

pdf(paste0(format(Sys.Date(), "%y%m%d"), "_", species_id, "-umap.pdf"), width = 8, height = 4)
ggplot(umap_df, aes(x = umap_1, y = umap_2, fill = cellType, color = cellType)) +
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2, 
                   arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  stat_ellipse(geom = "polygon", linetype = 2, alpha = 0.05, show.legend = FALSE, level = 0.96) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cluster_color_sub) +
  scale_fill_manual(values = cluster_color_sub) +
  ggtitle(paste0(nrow(umap_df), " cells")) +
  guides(
    fill = guide_legend(title = paste0(species_label, " subcluster")),
    color = guide_legend(title = paste0(species_label, " subcluster"))
  ) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = "bold", hjust = 0.03),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )
dev.off()

#============================#
#   Save Marker & Seurat Obj #
#============================#
Idents(species_sub) <- 'celltype_sub'
final_markers <- FindAllMarkers(species_sub, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.table(final_markers, file = "Ne.celltype.markers.xls", sep = "\t", row.names = FALSE, quote = TRUE)

saveRDS(species_sub, file = "Ne.rds")
