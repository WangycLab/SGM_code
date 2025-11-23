library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(igraph)
library(ggpubr)
library(gplots)
library(tidydr)
library(paletteer)

# set the working directory and load the Seurat object
seu <- readRDS("data/seurat_object_celltype.rds")
Idents(seu) <- "seurat_clusters"

# set color palettes
cluster_color <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
  "#dda0dd", "#8c564b", "#e377c2", "#7f7f7f", 
  "#636363", "#17becf", "#ffd700", "#bcbd22", 
  "#dcdcdc", "#ff9896", "#c49c94", "#ffbb78"
)

celltype_color <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
  "#dda0dd", "#8c564b", "#e377c2", "#7f7f7f"
)
order <- c("Polysaccharide degradation cells",
           "C4-dicarboxylate anaerobe cells",
           "Anaerobic ferredoxin-utilizing cells",
           "Agmatine biosynthesis cells",
           "Transposase+ cells",
           "Adhesive cells",
           "NorM efflux cells",
           "Outer membrane biogenesis cells")
names(celltype_color) <- order

color_map <- c(
  "green" = "ME1", "black" = "ME2", "yellow" = "ME3", "turquoise" = "ME4", 
  "blue" = "ME5",  "brown" = "ME6", "red" = "ME7", "grey" = "grey"
)

## Figure 3A -------------------------------------------------------------------
## Generate the UMAP plot for cluster

p <- DimPlot(seu, reduction = "umap", label = FALSE, group.by = "seurat_clusters",
             pt.size = 0.1, cols = cluster_color, raster = FALSE) +
  labs(x = "UMAP1", y = "UMAP2") + 
  theme(axis.text.y = element_blank(),
        text = element_text(size = 15, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
ggsave("plots/Figure3A.pdf", plot = p, width = 9, height = 8)

## Figures 3B & 3C -------------------------------------------------------------

# Step 1: data preparation for Figures 3B and 3C 
# (Recommend running on the server)

Idents(seu) <- "seurat_clusters"
set.seed(12345)

# prepare Seurat object for WGCNA, with gene selection based on fraction of cells expressed
seurat_obj <- SetupForWGCNA(seu, gene_select = "fraction", fraction = 0.005, wgcna_name = "tutorial")

# construct metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("seurat_clusters"), 
  k = 25, 
  max_shared = 10, 
  ident.group = "seurat_clusters"
)

# perform normalization on metacells
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars = "orig.ident")
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction = "harmony", dims = 1:15)

# set up the expression matrix for WGCNA analysis
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c(0:15),
  group.by="seurat_clusters",
  assay = "SCT", 
  slot = "data"
)

# select soft threshold for network construction
seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")

# extract the power table and select the optimal soft threshold
power_table <- GetPowerTable(seurat_obj)

# find the optimal soft threshold (based on a minimum R^2 value)
find_soft_threshold <- function(power_table, threshold = 0.8) {
  valid_rows <- power_table[power_table$SFT.R.sq >= threshold, ] 
  if (nrow(valid_rows) == 0) {
    stop("No results found that meet the condition SFT.R.sq >= ", threshold)
  }
  optimal_power <- valid_rows[1, "Power"] 
  return(optimal_power)
}
soft_threshold <- find_soft_threshold(power_table, threshold = 0.8)

# construct co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = soft_threshold, 
  setDatExpr = FALSE, 
  overwrite_tom = TRUE, 
  tom_name = "PD")

# identify module eigengenes and connectivity
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "orig.ident")
hMEs <- GetMEs(seurat_obj)
fwrite(hMEs, "data/Figure3C_WCGNA_feature.csv.gz", sep = ",")

# extract the UMAP coordinates of nodes (genes) in the network
g <- ModuleUMAPPlot(seurat_obj, return_graph = TRUE)

# create a dataframe with gene names and corresponding UMAP coordinates
umap_df <- data.frame(
  gene = igraph::V(g)$name,  
  UMAP1 = igraph::V(g)$UMAP1, 
  UMAP2 = igraph::V(g)$UMAP2, 
  module = igraph::V(g)$module,  
  kME = igraph::V(g)$kME,    
  hub = igraph::V(g)$hub
)

# extract edge information, including the start and end points of each edge
edges_df <- data.frame(
  from = igraph::ends(g, igraph::E(g))[, 1], 
  to = igraph::ends(g, igraph::E(g))[, 2]
)

# merge UMAP coordinates with edge data to display connectivity
edges_df <- merge(edges_df, umap_df, by.x = "from", by.y = "gene", all.x = TRUE)
edges_df <- merge(edges_df, umap_df, by.x = "to", by.y = "gene", suffixes = c(".from", ".to"))
write.csv(umap_df, "data/Figure3C_umap_nodes.csv", row.names = FALSE)
write.csv(edges_df, "data/Figure3C_umap_edges.csv", row.names = FALSE)

# Step 2: plotting

# load module eigengenes
Feature <- fread("data/Figure3C_WCGNA_feature.csv.gz", sep = ",", header = TRUE)
colnames(Feature) <- color_map[colnames(Feature)]
colnames(Feature)[1] <- "BC"
Feature <- Feature %>% dplyr::select(-grey)

# add labels to the Seurat object
seu@meta.data <- seu@meta.data %>%
  left_join(., Feature, by = "BC") %>%
  mutate(BC2 = BC) %>%
  tibble::column_to_rownames("BC2")

mods <- colnames(Feature)
mods <- mods[!mods %in% c("grey","BC")]

# generate the plot Figure 3B
p <- DotPlot(seu, features = paste0("ME", 1:7), group.by = "seurat_clusters") +
  scale_color_gradient2(low = "#7a0177", mid = "white", high = "darkred", midpoint = 0) +
  labs(x = "Cluster", y = "Module")+
  theme(
    axis.text.x = element_text(color = "black"),
    text = element_text(size = 12),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")
  ) 
ggsave("plots/Figure3B.pdf", plot = p, width = 6, height = 5)

# load the network information
umap_nodes <- read.csv("data/Figure3C_umap_nodes.csv")
umap_edges <- read.csv("data/Figure3C_umap_edges.csv")

# format the data for plotting
umap_nodes$module <- color_map[umap_nodes$module]
umap_edges$module.from <- color_map[umap_edges$module.from]
umap_edges$module.to <- color_map[umap_edges$module.to]
umap_nodes <- umap_nodes %>%
  dplyr::group_by(module) %>%
  dplyr::arrange(desc(kME)) %>%
  mutate(hub = case_when(
    module == "ME1" & gene %in% c("Pgu1", "BglS", "ManB2", "COG3942", "LytB") ~ "hub",
    module == "ME2" & gene %in% c("OprB", "GumC", "LpxD", "DR0291","SleL") ~ "hub",
    module == "ME3" & gene %in% c("BisC", "DcuA", "NadB", "AspA","SdhA" ) ~ "hub",
    module == "ME4" & gene %in% c( "OmpC", "Hia", "PutA1", "PutA2", "PutA") ~ "hub",
    module == "ME5" & gene %in% c("NuoI", "NapF", "PorD","CdhA","PreA") ~ "hub",
    module == "ME6" & gene %in% c("BtuB", "FepA", "CirA", "SpeA", "COG3291") ~ "hub",
    module == "ME7" & gene %in% c("TufA","FusA", "RpsC","SecY","RplB") ~ "hub",
    TRUE ~ "other"
  )) %>% dplyr::ungroup()

# generate the plot Figure 3C
p <- ggplot() +
  tidydr::theme_dr(
    xlength = 0.2, ylength = 0.2, 
    arrow = arrow(length = unit(0.2, "inches"), type = "closed")
  ) +
  geom_segment(
    data = umap_edges %>% dplyr::filter(module.from == module.to), 
    aes(x = UMAP1.from, y = UMAP2.from, xend = UMAP1.to, yend = UMAP2.to, color = module.from), 
    alpha = 0.3, size = 0.05
  ) +
  geom_point(
    data = umap_nodes %>% dplyr::filter(hub == "other"),
    aes(x = UMAP1, y = UMAP2, color = module, size = kME)
  ) +
  geom_point(
    data = umap_nodes %>% dplyr::filter(hub == "hub"), 
    aes(x = UMAP1, y = UMAP2, fill = module, size = 3),
    shape = 21, stroke = 0.5, color = "black"
  ) + 
  scale_color_viridis_d(guide = "none") + 
  scale_fill_viridis_d(name = "Module") + 
  geom_text_repel(
    data = umap_nodes %>% dplyr::filter(hub == "hub"), 
    aes(x = UMAP1, y = UMAP2, label = gene),
    size = 3, color = "black", max.overlaps = 10, box.padding = 0.3,  
   point.padding = 0.2,  
   segment.size = 0.2, 
   bg.color = "white",
   bg.r = 0.15) + 
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = "bold", hjust = 0.03),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  )+
  scale_size_continuous(name = "kME", range = c(0.005, 5)) +
  guides(fill = guide_legend(title = "Module"), size = "none")
ggsave("plots/Figure3C.pdf", plot = p, width = 10, height = 8)

## Figures 3D ------------------------------------------------------------------
## Plot marker genes

# list cell types
celltypes <- c("Polysaccharide degradation cells",
               "C4-dicarboxylate anaerobe cells",
               "Anaerobic ferredoxin-utilizing cells",
               "Agmatine biosynthesis cells",
               "Transposase+ cells",
               "Adhesive cells",
               "NorM efflux cells",
               "Outer membrane biogenesis cells") 
names(celltype_color) <- celltypes

# add annotations to the Seurat object
seu@meta.data <- seu@meta.data %>% 
  mutate(
    celltype = case_when(
      seurat_clusters %in% c(6) ~  "Polysaccharide degradation cells",
      seurat_clusters %in% c(7,11) ~  "C4-dicarboxylate anaerobe cells",
      seurat_clusters %in% c(2,12) ~  "Anaerobic ferredoxin-utilizing cells",
      seurat_clusters %in% c(0,8,13,9,15) ~  "Agmatine biosynthesis cells",
      seurat_clusters %in% c(1) ~  "Transposase+ cells",
      seurat_clusters %in% c(3,4,14) ~  "Adhesive cells",                                
      seurat_clusters %in% c(5) ~  "NorM efflux cells",
      seurat_clusters %in% c(10) ~ "Outer membrane biogenesis cells"
    ), 
    celltype = factor(celltype, levels = celltypes)
  )
Idents(seu) <- "celltype"

# generate the dot plot of marker genes
marker_list <- c(
  "BglS", "Pgu1", "ManB2", "DcuA", "AspA", "SdhA", "NosZ", "NapF", "NuoI", "PorB", 
  "SpeA", "BtuB",  "Tra8", "OmpC",  "Hia", "PilA", "NorM", "LpxD", "GumC"
)
dp <- DotPlot(seu, features = marker_list, assay = "SCT") + RotatedAxis()
df <- dp$data[,c(3,4,1,2,5)] %>% 
  dplyr::rename(gene = features.plot) %>%
  mutate(pct.exp = case_when(pct.exp < 25 ~ 0, TRUE ~ pct.exp))
p <- ggplot(df, aes(x = gene, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() + 
  scale_color_gradientn(
    colours = c("#2E5A87FF", "#6191BBFF", "white", "#E8544FFF", "#A90C38FF"),
    values = scales::rescale(c(min(seu$RNA@data), 0, max(seu$RNA@data)))
  ) +
  scale_size_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c("<25", "25", "50", "75", "100")
  )+
  scale_y_discrete(limits = rev(levels(seu@meta.data$celltype_sub))) + 
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black", angle = 90, hjust = 1),
    text = element_text(size = 10),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))+
  labs(size = "Percentage Expression", color = "Scaled Expression")+
  scale_color_gradient(low = "white", high = "darkred", limits = c(-1.0, 3))
ggsave("plots/Figure3D.pdf", plot = p, width = 9.8, height = 4)

## Figures 3E -------------------------------------------------------------------
## Generate the UMAP visualization of functional clusters

p <- DimPlot(seu, reduction = "umap", label = FALSE, group.by = "celltype",
             pt.size = 0.1, cols = celltype_color, raster = FALSE) +
  labs(x = "UMAP1", y = "UMAP2") + 
  theme(axis.text.y = element_blank(),
        text = element_text(size = 15 , face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
ggsave("plots/Figure3E.pdf", plot = p, width = 9, height = 5)
