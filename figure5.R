# =====================================================================
# Figure 5 — Subcluster organization, transcriptional suppression, and 
# functional reprogramming of *Neisseria elongata* in health and 
# periodontitis
# =====================================================================

library(ggplot2)
library(Seurat)
library(dplyr)
library(paletteer)
library(tidyr)
library(ggrepel)
library(tibble)
library(ggpubr)
library(ggalluvial)

# ---------------------------------------------------------------------
# Color settings
# ---------------------------------------------------------------------
cluster_color_sub <- c("#ffbb78", "#1f77b4", "darkred", "#636363", "#ffd700",
                       "#e377c2", "#ff7f0e", "#7f7f7f", "#bcbd22", "#dcdcdc")
group_color <- c("HC" = "#3D98D3FF", "Pd" = "#96281BFF")
cog_labels <- c(
  "J" = "J: Translation, ribosomal structure and biogenesis",
  "A" = "A: RNA processing and modification",
  "K" = "K: Transcription",
  "L" = "L: Replication, recombination and repair",
  "B" = "B: Chromatin structure and dynamics",
  "D" = "D: Cell cycle control, cell division, chromosome partitioning",
  "Y" = "Y: Nuclear structure",
  "V" = "V: Defense mechanisms",
  "T" = "T: Signal transduction mechanisms",
  "M" = "M: Cell wall/membrane/envelope biogenesis",
  "N" = "N: Cell motility",
  "Z" = "Z: Cytoskeleton",
  "W" = "W: Extracellular structures",
  "U" = "U: Intracellular trafficking, secretion, vesicular transport",
  "O" = "O: Posttranslational modification, protein turnover, chaperones",
  "C" = "C: Energy production and conversion",
  "G" = "G: Carbohydrate transport and metabolism",
  "E" = "E: Amino acid transport and metabolism",
  "F" = "F: Nucleotide transport and metabolism",
  "H" = "H: Coenzyme transport and metabolism",
  "I" = "I: Lipid transport and metabolism",
  "P" = "P: Inorganic ion transport and metabolism",
  "Q" = "Q: Secondary metabolites biosynthesis, transport, catabolism",
  "R" = "R: General function prediction only",
  "S" = "S: Function unknown",
  "X" = "X: Mobilome: prophages, transposons",
  "Others" = "Others"
)
cog_label_color <- c(paletteer_c(`"grDevices::Berlin"`,7),
                     paletteer_c(`"grDevices::Blues 2"`,5)[1:4],
                     paletteer_c(`"grDevices::Blues 3"`,8)[1:6],
                     paletteer_c(`"grDevices::Greens 2"`,5)[1:4],
                     paletteer_c(`"grDevices::Greens 3"`,8)[1:5],"grey")
names(cog_label_color) <- cog_labels

# ---------------------------------------------------------------------
# Define Wilcoxon test function for group comparison
# ---------------------------------------------------------------------
perform_wilcox_test <- function(data) {
  group1 <- data %>% filter(group == "Pd") %>% filter(!is.na(n)) %>% pull(n)
  group2 <- data %>% filter(group == "HC") %>% filter(!is.na(n)) %>% pull(n)
  if (length(group1) > 0 & length(group2) > 0) {
    test_result <- wilcox.test(group1, group2)
    return(test_result$p.value)
  } else {
    return(NA)
  }
}

order_sub <- c(
  "Twitching motility",
  "Translation regulation",
  "Nicotinamide metabolism",
  "Cell adhesion",
  "Signal transduction",
  "Binding-effector",
  "Polysaccharides degradation",
  "Antioxidant defense"
) 
names(cluster_color_sub) <- order_sub

# ---------------------------------------------------------------------
# UMAP visualization
# ---------------------------------------------------------------------
Roi <- "N.elongata"
Species_sub <- readRDS("data/N.elongata.rds")

# Extract UMAP embeddings and metadata
umap <- Species_sub@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(
    cellType = Species_sub@meta.data$celltype_sub,
    sample = Species_sub@meta.data$orig.ident,
    group = Species_sub@meta.data$group
  ) %>%
  mutate(cellType = factor(cellType, levels = order_sub))

# Generate and save UMAP plot
p <- ggplot(umap, aes(x = umap_1, y = umap_2, fill = cellType, color = cellType)) +
  
  # Add directional axes with arrowheads
  tidydr::theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = arrow(length = unit(0.2, "inches"), type = "closed")
  ) +
  
  # Add cluster confidence ellipses
  stat_ellipse(
    aes(x = umap_1, y = umap_2, fill = cellType),
    geom = "polygon",
    linetype = 2,
    alpha = 0.05,
    show.legend = FALSE,
    level = 0.96
  ) +
  
  # Add single-cell points
  geom_point(size = 0.1, show.legend = TRUE) +
  
  # Apply publication-style formatting
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = 2, hjust = 0.03),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  ) +
  
  # Set cluster color scheme
  scale_color_manual(values = cluster_color_sub) +
  scale_fill_manual(values = cluster_color_sub) +
  
  # Add title and legend
  ggtitle(paste0(nrow(umap), " cells")) +
  guides(
    fill = guide_legend(title = paste0(Roi, " subcluster")),
    color = guide_legend(title = paste0(Roi, " subcluster"))
  )
ggsave("plots/Figure5A.pdf", plot = p, width = 8, height = 4)

gene_roi <- c(
  "PilA","PilE","TufA","FusA","RplF","RpsC","PncA","TerY","COG5373",
  "PAS","WalK","Hia","LomR","PAAR","YukC","HybA1","XynB2","MsrB","Bcp"
)

p <- DotPlot(Species_sub, features = gene_roi, group.by = "celltype_sub") +
  scale_color_gradient2(low = "#3D98D3FF", mid = "white", high = "darkred", midpoint = 0) +
  labs(x = "Cluster", y = "Module")+
  theme(
    axis.text.x = element_text(color = "black", angle = 90, hjust = 1),
    text = element_text(size = 12),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")
  ) 
ggsave("plots/Figure5B.pdf", plot = p, width = 9.1, height = 4)

consist <- FetchData(Species_sub, vars = c("orig.ident", "celltype_sub")) %>%
  group_by(orig.ident, celltype_sub) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  complete(orig.ident, celltype_sub, fill = list(count = 0)) %>%
  group_by(orig.ident) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(n = count / total_count,
         group = case_when(grepl("HC",orig.ident) ~ "HC",
                           grepl("Pd",orig.ident) ~ "Pd"))
results <- consist %>%
  group_by(celltype_sub) %>%
  dplyr::summarise(p_value = perform_wilcox_test(cur_data())) %>%
  ungroup() %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "bonferroni"))
significant_clusters <- results %>%
  filter(p_value < 0.05) %>%
  pull(celltype_sub)
consist_plot_data <- Species_sub@meta.data %>%
  dplyr::select(celltype_sub, group) %>%
  mutate(celltype_sub = factor(celltype_sub, levels = order_sub))

P <- ggplot(consist_plot_data, aes(x = as.factor(group), fill = celltype_sub)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cluster_color_sub) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", face = "bold"),
        axis.title.y = element_text(color = "black", face = "bold"),
        text = element_text(size = 12)) +
  labs(x = "", y = paste0("Proportion of ",Roi," cells(%)"))
proportions <- consist_plot_data %>%  
  mutate(celltype_sub = factor(celltype_sub, levels = order_sub[length(order_sub):1])) %>%
  group_by(group, celltype_sub) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  dplyr::summarise(celltype_sub = celltype_sub,
                   total = sum(count),
                   proportion = count / total,
                   cumsum_lower = lag(cumsum(proportion), default = 0),
                   cumsum_upper = cumsum(proportion),
                   mid_point = cumsum_lower + (proportion / 2))
for (cluster in significant_clusters) {
  cluster_data <- proportions %>% filter(celltype_sub == cluster)
  max_group <- cluster_data$group[which.max(cluster_data$proportion)]
  mid_position <- cluster_data %>% filter(group == max_group) %>% pull(mid_point)
  print(paste("Cluster:", cluster, "Group:", max_group, "Mid position:", mid_position))
  
  P <- P + annotate("text",
                    x = max_group,
                    y = mid_position,
                    label = "*",
                    size = 5,
                    color = "black",
                    vjust = 0.75)
}
ggsave("plots/Figure5C.pdf", plot = P, width = 3.5, height = 5)

## Step 1: Data preparation for Figures 5D–5F
## (Recommended to run on an HPC server due to memory requirements)

# Load Seurat object for N. elongata
seu_sub <- readRDS("data/N.elongata.rds")
Idents(seu_sub) <- "celltype_sub"

# Initialize WGCNA framework with gene selection by expression fraction
set.seed(12345)
seurat_obj  <- SetupForWGCNA(
  seu_sub,
  gene_select = "fraction",     
  fraction = 0.005,              
  wgcna_name = "tutorial"
)

# Construct metacells (pseudo-bulk units) for each subcluster
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype_sub"),  
  k = 25,                        
  max_shared = 15,               
  ident.group = "celltype_sub"   
)

# Normalize, scale, and perform PCA on metacells
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features = VariableFeatures(seurat_obj))

# Correct sample-specific effects using Harmony
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars = "orig.ident")

# UMAP embedding of metacells
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction = "harmony", dims = 1:10)

# Prepare expression matrix for WGCNA
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Twitching motility","Translation regulation","Nicotinamide metabolism",
                 "Cell adhesion","Signal transduction","Binding-effector",
                 "Polysaccharides degradation","Antioxidant defense"), 
  group.by = "celltype_sub",     
  assay = "SCT",
  slot = "data"
)

# Identify soft-thresholding power for constructing scale-free network
seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
power_table <- GetPowerTable(seurat_obj)

# Select optimal soft threshold based on R^2 >= 0.8
find_soft_threshold <- function(power_table, threshold = 0.8) {
  valid_rows <- power_table[power_table$SFT.R.sq >= threshold, ]
  if (nrow(valid_rows) == 0) stop("No soft-power satisfies SFT.R.sq ≥ ", threshold)
  return(valid_rows[1, "Power"])
}
soft_threshold <- find_soft_threshold(power_table, 0.8)

# Build WGCNA network using TOM similarity
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = soft_threshold,
  setDatExpr = FALSE,
  overwrite_tom = TRUE,
  tom_name = "PD"
)

# Module eigengenes (MEs) and module connectivity (kME)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "orig.ident"   
)

# Export harmonized MEs
hMEs <- GetMEs(seurat_obj)
write.csv(hMEs, "Figure5_WCGNA_feature_Ne.csv")

# Compute eigengene-based connectivity (hubness score)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "celltype_sub",
  group_name = c("Twitching motility","Translation regulation","Nicotinamide metabolism",
                 "Cell adhesion","Signal transduction","Binding-effector",
                 "Polysaccharides degradation","Antioxidant defense")
)

# Extract module genes
modules <- GetModules(seurat_obj)

# UMAP layout for module network visualization
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors = 15,
  min_dist = 0.1
)

# Prepare gene–gene edge list for network plotting
edges_df <- data.frame(
  from = igraph::ends(g, igraph::E(g))[, 1],
  to = igraph::ends(g, igraph::E(g))[, 2]
)
edges_df <- merge(edges_df, umap_df, by.x = "from", by.y = "gene", all.x = TRUE)
edges_df <- merge(edges_df, umap_df, by.x = "to", by.y = "gene", suffixes = c(".from", ".to"))

# Save node and edge tables for Cytoscape or ggplot visualization
write.csv(umap_df, "data/Figure5_umap_nodes_Ne.csv", row.names = FALSE)
write.csv(edges_df, "data/Figure5_umap_edges_Ne.csv", row.names = FALSE)

## Figure 5D — module eigengene expression across subclusters

color_map <- c("brown" = "ME1", "blue" = "ME2", "turquoise" = "ME3", "grey" = "grey")

Feature <- read.csv("data/Figure5_WCGNA_feature_Ne.csv")
colnames(Feature) <- color_map[colnames(Feature)]
colnames(Feature)[1] <- "BC"

# Remove grey module (unassigned)
Feature <- Feature[, !colnames(Feature) %in% "grey"]

# Merge module eigengenes into original metadata
tmp <- left_join(Species_sub@meta.data, Feature, by = "BC")
rownames(tmp) <- tmp$BC
Species_sub@meta.data <- tmp

mods <- colnames(Feature)
mods <- mods[mods != "grey" & mods != "BC"]

p <- DotPlot(Species_sub, features = paste0("ME", 1:3), group.by = "celltype_sub", dot.scale = 10) +
  scale_color_gradientn(colours = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 5)) +
  labs(x = "Cluster", y = "Module") +
  theme(
    axis.text.x = element_text(color = "black"),
    text = element_text(size = 12),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", size = 1)
  )
ggsave("plots/Figure5D.pdf", plot = p, width = 8.3, height = 6)

## Figure 5E — Differential module eigengene analysis between HC and PD

umap_nodes <- read.csv("data/Figure5_umap_nodes_Ne.csv")
umap_edges <- read.csv("data/Figure5_umap_edges_Ne.csv")

# Map original module colors to ME labels
umap_nodes$module <- color_map[umap_nodes$module]

# Prepare assay with module eigengenes
MEs = Feature %>% column_to_rownames("BC")
MEs <- MEs[, colnames(MEs) != "grey"]
MEs[MEs < 0] <- 0
MEs <- t(MEs)
ME_assay <- Seurat::CreateAssayObject(MEs)

# Extract sample group
BCs <- colnames(MEs)
orig.ident <- sub("_.*", "", BCs)
Group <- dplyr::case_when(
  orig.ident %in% paste0("Pd", 1:8) ~ "PD",
  orig.ident %in% paste0("HC", 1:8) ~ "HC",
  TRUE ~ "Unknown"
)

tmp <- data.frame(BCs, orig.ident, Group)
group1 <- tmp %>% dplyr::filter(Group == "PD") %>% pull(BCs)
group2 <- tmp %>% dplyr::filter(Group == "HC") %>% pull(BCs)

# Compare ME expression between PD and HC
DMEs <- FindMarkers(
  ME_assay, cells.1 = group1, cells.2 = group2,
  slot = "counts", test.use = "wilcox",
  only.pos = FALSE, logfc.threshold = 0, min.pct = 0.15
)
DMEs$module <- rownames(DMEs)

# Volcano-style visualization of ME differences
p <- ggplot(DMEs, aes(y = module)) +
  geom_segment(aes(x = avg_log2FC, xend = 0, yend = module), color = "black") +
  geom_point(
    aes(x = avg_log2FC, size = abs(avg_log2FC), color = -log10(p_val))
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_viridis_c(option = "C") +
  labs(
    x = "Log2FC",
    y = "Module",
    size = "Log2FC (abs)",
    color = "-log10(P-value)"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.title = element_text(face = "bold")
  )
ggsave("plots/Figure5E.pdf", plot = p, width = 6, height = 4)

## Figure 5F — COG functional enrichment of hub genes in ME3

# Load COG annotation table
cog.20.def <- read.delim("data/cog-20.def.tab", header = FALSE)
colnames(cog.20.def) <- c("gene","COG_class","Func_anno","name","Func","time","Anno")
cog.20.def <- cog.20.def %>% mutate(name = ifelse(name == "", gene, name))

# Select hub genes in module ME3
top_ME3_genes <- umap_nodes[umap_nodes$module == "ME3" & umap_nodes$kME > 0.5, ]
top_ME3_genes$name <- top_ME3_genes$gene

# Join with COG annotation
top_ME3_genes_anno <- top_ME3_genes %>%
  left_join(cog.20.def, by = "name") %>%
  dplyr::filter(!is.na(COG_class)) %>%
  mutate(
    COG_class_label = cog_labels[COG_class],
    COG_class_label = if_else(nchar(COG_class) == 1, COG_class_label, "Others")
  )

# Pie chart data: collapse into top 10 COG categories + Others
colnames(top_ME3_genes) <- "name"
test1 <- top_ME3_genes %>%
  select(name) %>%
  left_join(cog.20.def, by = "name") %>%
  dplyr::filter(!is.na(COG_class)) %>%
  group_by(COG_class) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    percentage = round(count / sum(count) * 100, 3),
    COG_class_label = cog_labels[COG_class]
  ) %>%
  dplyr::filter(nchar(COG_class) == 1) %>%
  arrange(desc(percentage)) %>%
  mutate(COG_class_label = if_else(row_number() <= 10, COG_class_label, "Others")) %>%
  group_by(COG_class_label) %>%
  summarise(count = sum(count), percentage = sum(percentage), .groups = "drop") %>%
  arrange(desc(percentage))

# Set factor ordering for plotting
test1$COG_class_label <- factor(
  test1$COG_class_label,
  levels = c(test1$COG_class_label[1:length(test1$COG_class_label)], "Others")
)

# Draw pie chart
p <- ggplot(test1, aes(x = 2, y = percentage, fill = COG_class_label)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = paste0(round(percentage, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(fill = "COG Class") +
  scale_fill_manual(values = cog_label_color)
ggsave("plots/Figure5F.pdf", plot = p, width = 8, height = 5)


## Figure 5G — Representative genes from ME3 showing differential expression between HC and Pd

# Define group labels for comparison
group1 <- "HC"
group2 <- "Pd"

# Assign colors used in the plotting scheme
colors <- c("HC" = "#3D98D3FF", "Pd" = "#96281BFF")

# Retrieve cell identities for the two groups
group1_cells <- WhichCells(Species_sub, expression = group == group1)
group2_cells <- WhichCells(Species_sub, expression = group == group2)

# Select representative ME3 marker genes for visualization
gene_list <- c("RplB", "RpsC", "SecY", "InfA", "NusA", "FusA")

# Loop through each selected gene and generate violin/box/jitter plots
for (gene in gene_list) {
  
  # Extract normalized expression (RNA slot "data") from both groups
  expr_matrix <- GetAssayData(Species_sub, slot = "data")[gene, c(group1_cells, group2_cells)]
  
  # Organize expression values into a plotting-friendly data frame
  expression_data <- data.frame(
    Expression = as.numeric(expr_matrix),
    group = factor(
      c(rep(group1, length(group1_cells)),
        rep(group2, length(group2_cells))),
      levels = c(group1, group2)
    )
  )
  
  # Compute Wilcoxon rank-sum test between HC and Pd
  wilcox_test_result <- wilcox.test(Expression ~ group, data = expression_data)
  
  # Generate violin plot with boxplot overlay and jittered non-zero points
  p <- ggplot(expression_data, aes(x = group, y = Expression, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.15, color = "black", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(
      data = expression_data,
      fill = "darkgray",     
      color = "black",    
      width = 0.33,
      size = 1.6,
      alpha = 0.4,
      stroke = 0.1,
      shape = 21          
    ) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(gene, " expression"),
      y = "Expression Level"
    ) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c(group1, group2)),
      label = "p.signif",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      )
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(expression_data$Expression) * 1.1,
      label = paste("p =", format(wilcox_test_result$p.value,
                                  scientific = TRUE, digits = 3)),
      size = 5
    )
  
  # Output plot and Wilcoxon statistics to console
  print(p)
  print(wilcox_test_result)
  
  # Save panel as a high-resolution PDF
  ggsave(
    filename = paste0("plots/Figure5G_", gene, ".pdf"),
    plot = p,
    width = 4.0, height = 5, dpi = 600,
    device = cairo_pdf
  )
}
