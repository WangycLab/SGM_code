# ====== 1. Load Required Libraries ======
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
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
  library(clustree)
})

# Load preprocessed Seurat object
Combined_seurat_object <- readRDS('Combined_alter.rds')

# ====== Step 1: WGCNA Preparation ====== (Recommend running on the server/cluster)
Idents(Combined_seurat_object) <- 'seurat_clusters'
set.seed(12345)

# Prepare Seurat object for WGCNA, with gene selection based on fraction of cells expressed
seurat_obj  <- SetupForWGCNA(
  Combined_seurat_object,
  gene_select = "fraction", 
  fraction = 0.005, 
  wgcna_name = "tutorial"
)

# Construct Metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("seurat_clusters"), 
  k = 25, 
  max_shared = 10, 
  ident.group = 'seurat_clusters'
)

# Normalize metacell expression matrix
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars = 'orig.ident')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction = 'harmony', dims = 1:15)

# Gene Expression Setup for WGCNA
# Set up the expression matrix for WGCNA analysis
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c(0:14),
  group.by='seurat_clusters',
  assay = 'SCT', 
  slot = 'data'
)

# Select Soft Threshold for Network Construction
seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'signed')

# Extract the power table and select the optimal soft threshold
power_table <- GetPowerTable(seurat_obj)

# Function to find the optimal soft threshold (based on a minimum R^2 value)
find_soft_threshold <- function(power_table, threshold = 0.8) {
  valid_rows <- power_table[power_table$SFT.R.sq >= threshold, ] 
  if (nrow(valid_rows) == 0) {
    stop("No results found that meet the condition SFT.R.sq >= ", threshold)
  }
  optimal_power <- valid_rows[1, "Power"] 
  return(optimal_power)
}

# Find the soft threshold with R^2 value >= 0.8
soft_threshold <- find_soft_threshold(power_table, threshold = 0.8)

# Construct Co-expression Network
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = soft_threshold, 
  setDatExpr=FALSE, 
  overwrite_tom = TRUE, 
  tom_name = 'PD'
)

# Module Eigengenes and Connectivity
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "orig.ident")
hMEs <- GetMEs(seurat_obj)

# Calculate Connectivity and Hub Genes
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'seurat_clusters', 
  group_name = c(0:16)
)
modules <- GetModules(seurat_obj)

# Identify the hub genes in each module
hub_df <- GetHubGenes(seurat_obj, n_hubs = 20)

# Apply UMAP to Co-expression Network
# Run UMAP dimensionality reduction for the module eigengene network
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, 
  n_neighbors = 15, 
  min_dist = 0.1
)
umap_df <- GetModuleUMAP(seurat_obj)

# Extract the UMAP coordinates of nodes (genes) in the network
g <- ModuleUMAPPlot(seurat_obj, return_graph = TRUE)

# Create a dataframe with gene names and corresponding UMAP coordinates
umap_df <- data.frame(
  gene = igraph::V(g)$name,  
  UMAP1 = igraph::V(g)$UMAP1, 
  UMAP2 = igraph::V(g)$UMAP2, 
  module = igraph::V(g)$module,  
  kME = igraph::V(g)$kME,    
  hub = igraph::V(g)$hub
)

# Extract edge information, including the start and end points of each edge
edges_df <- data.frame(
  from = igraph::ends(g, igraph::E(g))[, 1], 
  to = igraph::ends(g, igraph::E(g))[, 2]
)

# Merge UMAP coordinates with edge data to display connectivity
edges_df <- merge(edges_df, umap_df, by.x = "from", by.y = "gene", all.x = TRUE)
edges_df <- merge(edges_df, umap_df, by.x = "to", by.y = "gene", suffixes = c(".from", ".to"))

# ====== Step 2: DotPlot of Module Features ======
color_map <- c("green" = "ME1","black" = "ME2",'yellow' = "ME3", "turquoise" = "ME4", 
               "blue" = "ME5", "brown" = "ME6","red" = "ME7","grey" = "grey")
colnames(Feature) <- color_map[colnames(Feature)]
colnames(Feature)[1] <- 'BC'
Feature <- Feature[, !colnames(Feature) %in% "grey"]
tmp <- left_join(Combined_seurat_object@meta.data, Feature, by = 'BC')
rownames(tmp) <- tmp$BC
Combined_seurat_object <- AddMetaData(Combined_seurat_object, metadata = tmp)

mods <- colnames(Feature)
mods <- mods[mods %in% c('grey','BC') == F]

pdf('WGCNA_dotlot.pdf', width = 6, height = 5)
DotPlot(Combined_seurat_object, features= paste0('ME',1:7), group.by = 'seurat_clusters') +
  scale_color_gradient2(low = "#7a0177", mid = "white", high = "darkred", midpoint = 0) +
  labs(x = "Cluster", y = "Module")+
  theme(
    axis.text.x = element_text(color = 'black'),
    text = element_text(size = 12),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")
  )
dev.off()

umap_nodes$module <- color_map[umap_nodes$module]
umap_edges$module.from <- color_map[umap_edges$module.from]
umap_edges$module.to <- color_map[umap_edges$module.to]
umap_nodes <- umap_nodes %>%
  group_by(module) %>%
  arrange(desc(kME)) %>%
  mutate(hub = case_when(
    module == 'ME1' & gene %in% c('Pgu1', 'BglS', 'LytB', 'COG3942','ManB2') ~ 'hub',
    module == 'ME2' & gene %in% c('OprB', 'DR0291', 'GumC', 'LpxD','SleL') ~ 'hub',
    module == 'ME3' & gene %in% c('PorA-2', 'BisC', 'CytC552', 'YjgC', 'DcuA' ) ~ 'hub',
    module == 'ME4' & gene %in% c( 'OmpC', 'PutA', 'PutA1','Hia','PutA2') ~ 'hub',
    module == 'ME5' & gene %in% c('NuoI', 'NapF', 'CdhA', 'PreA', 'PorD') ~ 'hub',
    module == 'ME6' & gene %in% c('BtuB', 'FepA', 'CirA', 'SpeA', 'COG3291') ~ 'hub',
    module == 'ME7' & gene %in% c('TufA','FusA', 'RpsC','SecY','RplB') ~ 'hub',
    TRUE ~ 'other'
  )) %>% ungroup()

pdf('dotumap.pdf', width = 10, height = 8)
ggplot() +
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  geom_segment(data = umap_edges %>% dplyr::filter(module.from == module.to), 
               aes(x = UMAP1.from, y = UMAP2.from, xend = UMAP1.to, yend = UMAP2.to, color = module.from), 
               alpha = 0.3, size = 0.05) +
  geom_point(data = umap_nodes %>% dplyr::filter(hub == 'other'),
             aes(x = UMAP1, y = UMAP2, color = module, size = kME)) +
  geom_point(data = umap_nodes %>% dplyr::filter(hub == 'hub'), 
             aes(x = UMAP1, y = UMAP2, fill = module, size = 3),
             shape = 21, stroke = 0.5, color = "black") + 
  scale_color_viridis_d(guide = "none") + scale_fill_viridis_d(name = "Module") + 
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    aspect.ratio = 1, 
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = "bold", hjust = 0.03),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  scale_size_continuous(name = "kME", range = c(0.005, 5)) +
  guides(fill = guide_legend(title = "Module"), size = 'none')
dev.off()

# ====== Step 3: Define Cluster Annotations ======
order = c('Polysaccharide degradation cells',
          'C4-dicarboxylate anaerobe cells',
          'Anaerobic ferredoxin-utilizing cells',
          'Agmatine biosynthesis cells',
          'Transposase+ cells',
          'Signal communication cells',
          'NorM efflux cells',
          'Outer membrane biogenesis cells') 
names(cluster_color) <- order
Combined_seurat_object@meta.data <- Combined_seurat_object@meta.data %>% 
  mutate(
    celltype = case_when(
      seurat_clusters %in% c(6) ~  "Polysaccharide degradation cells",
      seurat_clusters %in% c(7,11) ~  'C4-dicarboxylate anaerobe cells',
      seurat_clusters %in% c(2,12) ~  'Anaerobic ferredoxin-utilizing cells',
      seurat_clusters %in% c(0,8,9,13,15) ~  "Agmatine biosynthesis cells",
      seurat_clusters %in% c(1) ~  'Transposase+ cells',
      seurat_clusters %in% c(3,4,14) ~  "Signal communication cells",                                
      seurat_clusters %in% c(5) ~  "NorM efflux cells",
      seurat_clusters %in% c(10) ~ "Outer membrane biogenesis cells"
    ), 
    celltype = factor(celltype,levels = order)
  )
Idents(Combined_seurat_object) <- 'celltype'
all.cluster.markers <- FindAllMarkers(Combined_seurat_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.cluster.markers %>% group_by(cluster), file ="celltype_marker.xls", sep = "\t")




