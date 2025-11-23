library(ggplot2)
library(Seurat)
library(dplyr)
library(tibble)
library(cowplot)
library(tidyr)
library(ggpubr)
library(muscat)
library(UpSetR)
library(paletteer)
library(pheatmap)
library(scRNAtoolVis)

# --- Set working directory and load Seurat object -----------------------------
seu <- readRDS("data/seurat_object_celltype.rds")

# --- Metadata and color definitions -------------------------------------------
seu@meta.data$sample <- seu@meta.data$orig.ident
name_list <- c(paste0("HC", 1:8), paste0("Pd", 1:8))

# Define cell type colors (consistent with main figures)
celltype_color <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
  "#dda0dd", "#8c564b", "#e377c2", "#7f7f7f"
)
celltypes <- c(
  "Outer membrane biogenesis cells", 
  "NorM efflux cells", 
  "Adhesive cells", 
  "Transposase+ cells",
  "Agmatine biosynthesis cells",
  "Anaerobic ferredoxin-utilizing cells",
  "C4-dicarboxylate anaerobe cells",
  "Polysaccharide degradation cells"
) 
names(celltype_color) <- rev(celltypes)

# Define group colors (Health vs Periodontitis)
group_color <- c("HC" = "#3D98D3FF", "Pd" = "#96281BFF")

# Define COG category colors 
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
  "Others" = "Others")
cog_label_color <- c(paletteer_c(`"grDevices::Berlin"`,7),
                     paletteer_c(`"grDevices::Blues 2"`,5)[1:4],
                     paletteer_c(`"grDevices::Blues 3"`,8)[1:6],
                     paletteer_c(`"grDevices::Greens 2"`,5)[1:4],
                     paletteer_c(`"grDevices::Greens 3"`,8)[1:5], "grey")
names(cog_label_color) <- cog_labels

# --- Define statistical test (Wilcoxon rank-sum test) --------------------------------
perform_wilcox_test <- function(data) {
  group1 <- data %>% dplyr::filter(group == "Pd") %>% dplyr::filter(!is.na(n)) %>% pull(n)
  group2 <- data %>% dplyr::filter(group == "HC") %>% dplyr::filter(!is.na(n)) %>% pull(n)
  if (length(group1) > 0 & length(group2) > 0) {
    test_result <- wilcox.test(group1, group2)
    return(test_result$p.value)
  } else {
    return(NA)
  }
}

########################################################################################
## Figure 4A–B. Visualization of functional subpopulations and compositional shifts
## -------------------------------------------------------------------------------------
## Panel A–B:  UMAP visualization of bacterial transcriptional subpopulations 
##             in health (HC) and periodontitis (Pd)
## Right panel: Proportion of each functional cluster across groups,
##              with statistical annotation (* p < 0.05)
## -------------------------------------------------------------------------------------
## Output: Figure4A_and_B.pdf (15 × 4.5 in)
########################################################################################

# --- 1. Prepare UMAP embedding and metadata ------------------------------------------

umap <- seu@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cellType = seu@meta.data$celltype,
        sample = seu@meta.data$orig.ident,
        group = seu@meta.data$group) %>%
  mutate(cellType = factor(cellType, levels = celltypes))
groups <- c("HC","Pd")
umap_list <- lapply(groups, function(grp) subset(umap, group == grp))
plots <- list()

# --- 2. Plot UMAP for each condition (HC and Pd) -------------------------------------

for (i in 1:length(groups)) {
  umap_group <- umap_list[[i]]
  
  p <- ggplot(umap_group, aes(x = UMAP_1, y = UMAP_2, fill = cellType, color = cellType)) +
    tidydr::theme_dr(xlength = 0.2,
                     ylength = 0.2,
                     arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
    stat_ellipse(aes(x = UMAP_1, y = UMAP_2, fill = cellType),
                 geom = "polygon",  
                 linetype = 2,      
                 alpha = 0.05,     
                 show.legend = FALSE, 
                 level = 0.96) +    
    geom_point(size = 0.1, show.legend = F) +
    theme(
      aspect.ratio = 1,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(arrow = arrow(type = "closed")),
      axis.title = element_text(face = 2, hjust = 0.03),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12) 
    ) +
    guides(fill = guide_legend(override.aes = list(colour = "white", linetype = 1, size = 1))) +
    scale_color_manual(values = celltype_color) +
    scale_fill_manual(values = celltype_color) +
    theme(legend.position = "none") +
    ggtitle(paste0(groups[i],":",nrow(umap_group)," cells"))
  
  plots[[i]] <- p
}

# --- 3. Generate unified legend -------------------------------------------------------

legend_plot <- ggplot(umap_list[[1]], aes(x = UMAP_1, y = UMAP_2, fill = cellType, color = cellType)) +
  geom_point(size = 1, show.legend = T) +
  scale_color_manual(values = celltype_color) +
  scale_fill_manual(values = celltype_color) +
  labs(color = "Functional cluster", fill = "Functional cluster") +
  theme(legend.position = "right",
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(face = "bold")) 
legend <- get_legend(legend_plot)

# --- 4. Calculate cluster proportions per sample -------------------------------------

consist <- FetchData(seu, vars = c("orig.ident", "celltype")) %>%
  group_by(orig.ident, celltype) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  complete(orig.ident, celltype, fill = list(count = 0)) %>%
  group_by(orig.ident) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(n = count / total_count)  %>%
  mutate(group = case_when(orig.ident %in% paste0("Pd", 1:8) ~ "Pd", TRUE ~ "HC"))

# --- 5. Perform statistical comparison (Wilcoxon test per cell type) -----------------

results <- consist %>% 
  group_by(celltype) %>%
  dplyr::summarise(p_value = perform_wilcox_test(cur_data())) %>%
  ungroup() %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "bonferroni"))

significant_clusters <- results %>%
  dplyr::filter(adjusted_p_value < 0.05) %>%
  pull(celltype)

# --- 6. Generate composition bar plot ------------------------------------------------

consist_plot_data <- seu@meta.data %>%
  dplyr::select(celltype, group) %>%
  mutate(celltype = factor(celltype, levels = celltypes)) 

P3 <- ggplot(consist_plot_data, aes(x = as.factor(group), fill = celltype)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = celltype_color) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", face = "bold"),
        text = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "Proportion(%)", fill = "Functional cluster")

# --- 7. Annotate significant clusters (*) --------------------------------------------

proportions <- consist_plot_data %>%  
  mutate(celltype = factor(celltype, levels = celltypes[length(celltypes):1])) %>%
  group_by(group, celltype) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  dplyr::summarise(celltype = celltype,
                   total = sum(count),
                   proportion = count / total,
                   cumsum_lower = lag(cumsum(proportion), default = 0),
                   cumsum_upper = cumsum(proportion),
                   mid_point = cumsum_lower + (proportion / 2))

for (cluster in significant_clusters) {
  cluster_data <- proportions %>% filter(celltype == cluster)
  max_group <- cluster_data$group[which.max(cluster_data$proportion)]
  mid_position <- cluster_data %>% filter(group == max_group) %>% pull(mid_point)
  # Debug: Check positions
  print(paste("Cluster:", cluster, "Group:", max_group, "Mid position:", mid_position))
  
  P3 <- P3 + annotate("text",
                      x = max_group,
                      y = mid_position,
                      label = "*",
                      size = 5,
                      color = "black",
                      vjust = 0.75)
}

# --- 8. Combine plots into final composite figure ------------------------------------

pdf("plots/Figure4A_B.pdf", width = 15, height = 4.5)
plot_grid(plot_grid(plotlist = plots, ncol = length(groups)), legend, P3, ncol = 3, rel_widths = c(1, 0.25,0.3))
dev.off()

#-------------------------------------------------------------
# Figure 4C: Comparison of expression levels of selected genes
# (AceE, AceF, AdhP, EutG, LdhA) between healthy (HC) and periodontitis (Pd) samples
#-------------------------------------------------------------

# Define genes of interest
gene_roi <- c("AceE", "AceF", "AdhP", "EutG", "LdhA")

# Set the default assay to RNA to extract normalized expression values
DefaultAssay(seu) <- "SCT"

# Fetch expression data for the selected genes and sample identifiers
data_expr <- FetchData(seu, vars = c("orig.ident", gene_roi))

# Assign sample groups (HC or Pd) based on sample IDs
data_expr$group <- ifelse(grepl("^HC", data_expr$orig.ident), "HC", "Pd")

# Reshape data into long format (one row per sample–gene pair)
data_expr_long <- data_expr %>%
  pivot_longer(cols = all_of(gene_roi), names_to = "Gene", values_to = "Expression") %>%
  # Calculate average expression per sample for each gene
  group_by(orig.ident, group, Gene) %>%
  summarise(Average_Expression = mean(Expression), .groups = "drop")

# Perform Kruskal–Wallis test (non-parametric) to compare HC vs Pd for each gene
signif_results <- data_expr_long %>%
  group_by(Gene) %>%
  summarise(p_value = kruskal.test(Average_Expression ~ group)$p.value) %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# Merge statistical significance information back into expression data
data_expr_long <- left_join(data_expr_long, signif_results, by = "Gene")

# Set gene factor order for plotting
data_expr_long$Gene <- factor(data_expr_long$Gene, levels = gene_roi)

# Compute vertical position for significance labels
label_pos <- data_expr_long %>%
  group_by(Gene) %>%
  summarise(y_pos = max(Average_Expression) * 1.05) %>%
  left_join(signif_results, by = "Gene") %>%
  dplyr::filter(Significance != "")

#-------------------------------------------------------------
# Plot: Boxplots of average gene expression per group
#-------------------------------------------------------------

p <- ggplot(data_expr_long, aes(x = group, y = Average_Expression, fill = group)) +
  # Boxplot (no outliers shown)
  geom_boxplot(outlier.shape = NA, width = 0.54, color = "black", size = 0.4) +
  # Overlay jittered points for individual samples
  geom_jitter(aes(color = group), width = 0.2, size = 3, shape = 1, stroke = 0.4) +
  # Facet by gene, allowing free y-axis scales
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  # Add white-filled points for group means
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") +
  # Clean minimalist theme
  # theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # Axis labels and title
  labs(
    title = "Average Gene Expression by Group (HC vs Pd)",
    x = "Group",
    y = "Average Expression"
  ) +
  # Group color scheme
  scale_fill_manual(values = c("HC" = "#3D98D3FF", "Pd" = "#96281BFF")) +
  scale_color_manual(values = c("HC" = "black", "Pd" = "black")) +
  # Add significance symbols
  geom_text(
    data = label_pos,
    aes(x = 1.5, y = y_pos, label = Significance),
    inherit.aes = FALSE,
    size = 6
  )
ggsave("plots/Figure4C.pdf", plot = p, width = 10.5, height = 6)

########################################################################################
## Figure 4D: Differential expression (muscat pbDS) and volcano-style visualization
## - Aggregates single-cell counts to pseudo-bulk per cluster×sample
## - Runs muscat pbDS differential testing (Pd vs HC)
## - Filters DEGs (p_val < 0.05), annotates with COG definitions, and plots with jjVolcano
## Output: Figure4D_DEG.csv and Figure4D.pdf
########################################################################################

# --- 1. Read COG annotation table -----------------------------------------------------
# The file contains mapping from gene -> functional COG annotation.
cog.20.def <- read.delim("data/cog-20.def.tab", header = FALSE)
colnames(cog.20.def) <- c("gene", "COG_class", "Func_anno", "name", "Func", "time", "Anno")

# Replace empty "name" entries with the gene identifier for robust joins later.
cog.20.def <- cog.20.def %>% mutate(name = ifelse(name == "", gene, name))

# --- 2. Prepare Seurat -> SingleCellExperiment for muscat ----------------------------
# Use RNA assay for DE testing 

DefaultAssay(seu) <- "RNA"

# Convert Seurat object to SingleCellExperiment
seu_sce <- as.SingleCellExperiment(seu)

# Prepare SCE for muscat: set cluster id (kid), group id (gid) and sample id (sid).
# drop = TRUE removes unused factor levels
seu_sce <- prepSCE(seu_sce,
                   kid = "celltype",    # cluster identifier column in colData
                   gid = "group",       # group identifier column (HC vs Pd)
                   sid = "orig.ident",  # sample identifier column
                   drop = TRUE)

# Ensure group factor levels are ordered consistently (HC then Pd)
seu_sce$group_id <- factor(seu_sce$group_id, levels = c("HC", "Pd"))

# --- 3. Aggregate counts to pseudo-bulk per cluster x sample --------------------------
# aggregateData returns a list-like object "pb" with aggregated counts for each cluster/sample
pb <- aggregateData(seu_sce,
                    assay = "counts",
                    fun = "sum",
                    by = c("cluster_id", "sample_id"))

# --- 4. Differential state testing using muscat pbDS ---------------------------------
# pbDS runs differential testing on pseudo-bulk data; result structure in `res`
res <- muscat::pbDS(pb)

# Extract differential testing table for the contrast "Pd" (muscat stores results by group)
test <- res[["table"]][["Pd"]]

# --- 5. Collect significant DEGs and annotate ---------------------------------------
# Combine returned tables (if multiple), filter by p-val, rename fields for clarity,
# and join with COG annotations.
combined_DEG <- do.call(rbind, res[["table"]][["Pd"]]) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::rename(
    avg_log2FC = logFC,    # rename muscat log fold-change column
    cluster    = cluster_id,
    p_val_adj  = p_adj.loc,
    name       = gene      # gene name column to match cog table
  ) %>%
  dplyr::select(cluster, name, avg_log2FC, p_val, p_val_adj) %>%
  dplyr::arrange(cluster, avg_log2FC) %>%
  # Join to COG annotation; many-to-many relationship allowed
  left_join(cog.20.def, by = "name", relationship = "many-to-many")

# Save DEG table for downstream inspection and reporting
write.csv(combined_DEG, "data/Figure4D_DEG.csv", row.names = FALSE)

# Ensure cluster factor levels follow the pre-defined "order" used for figures
combined_DEG$cluster <- factor(combined_DEG$cluster, levels = celltypes)

# Re-assign names for celltype_color to guarantee consistent mapping
names(celltype_color) <- celltypes

# --- 6. Volcano-style plotting using jjVolcano --------------------------------------
# Plot sized for publication (12.8 × 8 in)
pdf("plots/Figure4D.pdf", width = 12.8, height = 8)

# jjVolcano visualizes DEGs per cluster; parameters:
jjVolcano(diffData = combined_DEG,
          log2FC.cutoff = 0.25,
          size = 2.6,
          aesCol = c("#0072B2", "#D55E00"),
          tile.col = celltype_color[celltypes],
          topGeneN = 3) +
  theme(legend.position = c(0.8, 0.9))

dev.off()


# Define a list to store DEGs for each functional cluster
lt <- list()
for (i in celltypes) {
  tmp <- combined_DEG %>% dplyr::filter(cluster == i) %>% pull(gene)
  lt[[i]] <- unique(tmp)
}

# Remove NAs
lt_clean <- lapply(lt, function(x) x[!is.na(x)])

# Assign cluster names to the color vector
names(celltype_color) <- celltypes

# Generate an UpSet plot showing intersection of DEGs among clusters
pdf("plots/Figure4E_upset.pdf", width = 18, height = 10)
upset(fromList(lt_clean),
      nsets = length(lt),
      nintersects = 30,
      keep.order = TRUE,
      number.angles = 0,
      point.size = 4,
      line.size = 0.6,
      mainbar.y.label = "Intersection size",
      main.bar.color = "black",
      matrix.color = "black",
      sets.x.label = "Set size",
      sets.bar.color = celltype_color[c(
        "Polysaccharide degradation cells",
        "C4-dicarboxylate anaerobe cells",
        "Anaerobic ferredoxin-utilizing cells",
        "Agmatine biosynthesis cells",
        "Transposase+ cells",
        "Adhesive cells",
        "NorM efflux cells",
        "Outer membrane biogenesis cells"
      )] %>% unname(),
      mb.ratio = c(0.6, 0.4),
      order.by = "freq")
dev.off()

#-------------------------------------------
# Pie chart of COG functional distribution
#-------------------------------------------

# Summarize DEGs by COG class and calculate percentages
test1 <- combined_DEG %>%
  dplyr::select(cluster, name) %>%
  group_by(name) %>%
  dplyr::summarise(count = n()) %>%
  left_join(cog.20.def, by = "name") %>%
  dplyr::filter(!is.na(COG_class)) %>%
  group_by(COG_class) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round(count / sum(count) * 100, 3),
         COG_class_label = cog_labels[COG_class]) %>%
  dplyr::filter(nchar(COG_class) == 1) %>%
  arrange(desc(percentage)) %>%
  mutate(COG_class_label = if_else(row_number() <= 10, COG_class_label, "Others")) %>%
  group_by(COG_class_label) %>%
  dplyr::summarise(count = sum(count),
                   percentage = sum(percentage),
                   .groups = "drop") %>%
  arrange(desc(percentage))

# Set factor levels to ensure “Others” appears last
test1$COG_class_label <- factor(
  test1$COG_class_label,
  levels = c(test1$COG_class_label[2:length(test1$COG_class_label)], "Others")
)

# Generate a pie plot of COG functional categories among DEGs
p <- ggplot(test1, aes(x = 2, y = percentage, fill = COG_class_label)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3) +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "", fill = "COG Class") +
  scale_fill_manual(values = cog_label_color)
ggsave("plots/Figure4E_pie.pdf", plot = p, width = 8, height = 5)

tmp <- combined_DEG %>% 
  dplyr::filter(p_val< 0.05) %>%
  dplyr::select(name, cluster, avg_log2FC, COG_class, Func_anno, Func) %>%
  arrange(cluster, avg_log2FC) 

## Figure 4F --------------------------------------------------------------------------
# Expression levels of TCA cycle genes across functional clusters

# Select TCA cycle genes from DEGs
gene_roi <- unique(c(tmp %>% dplyr::filter(grepl("TCA",Func)) %>% pull(name)))

# change order
gene_roi
gene_roi <- c("FumC","Mdh","AcnB","AcnA","LeuB","SucA","GltA","SucC","SdhA","SdhB/FrdB")

# Prepare a matrix of avg_log2FC values for selected genes across clusters
matrix <- tmp %>%
  dplyr::filter(name %in% gene_roi) %>%
  arrange(Func) %>%
  complete(name, cluster = celltypes, fill = list(avg_log2FC = 0)) %>%
  dplyr::select(-COG_class, -Func_anno, -Func) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = list(avg_log2FC = 0)) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Reorder matrix rows and columns for consistent cluster order
matrix <- matrix[gene_roi, celltypes]

# Prepare column annotation (cluster identity)
annotation_col <- data.frame(CellType = colnames(matrix)) 
rownames(annotation_col) <- colnames(matrix)
CellType_colors <- celltype_color[1:length(celltypes)]
names(CellType_colors) <- celltypes

# Define annotation color list for pheatmap
annotation_colors <- list(CellType = CellType_colors)

# Generate color gradient for heatmap
annotation_colors <- list(CellType = CellType_colors)
min_value <- min(matrix, na.rm = TRUE)
max_value <- max(matrix, na.rm = TRUE)
gradient_colors <- colorRampPalette(c("darkblue", "white", "salmon", "darkred"))(100)

# Define breaks for color scaling and ensure zero is included
breaks <- seq(min_value, max_value, length.out = 54.6)
zero_index <- which.min(abs(breaks - 0)) 
breaks[zero_index] <- 0  

# Plot heatmap without clustering, row/column names, or legend
pdf("plots/Figure4F.pdf", width = 4, height = 7)
pheatmap(
  matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_colors = annotation_colors,
  annotation_col = annotation_col,
  color = gradient_colors,
  breaks = breaks,
  annotation_legend = FALSE
)
dev.off()

## Figure 4G --------------------------------------------------------------------------
# TCA cycle module score across samples

# Define TCA cycle gene set for module score calculation
TCA_feature <- cog.20.def %>% dplyr::filter(grepl("TCA", Func)) %>% pull(name)

# Add TCA cycle module score to Seurat object
seu <- AddModuleScore(seu,
                      features = list(TCA_feature),
                      ctrl = 100,
                      name = "TCA_features")

# Aggregate module scores by sample
test <- seu@meta.data %>% 
  dplyr::select(sample, TCA_features1) %>%
  group_by(sample) %>%
  dplyr::summarise(TCA = mean(TCA_features1))

# Annotate group (HC or Pd) based on sample ID
test$group <- ifelse(grepl("^Pd", test$sample), "Pd", "HC")

# Genearate the box plot of TCA cycle module scores by group
p <- ggplot(test, aes(x = group, y = TCA, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, shape = 16, size = 2) +
  geom_signif(comparisons = list(c("HC", "Pd")),
              map_signif_level = TRUE,
              y_position = max(test$TCA) + 0.05) +
  labs(title = "", x = "", y = "TCA cycle module score") +
  theme_classic() +
  scale_fill_manual(values = group_color) +
  guides(fill = guide_legend(title = "Group"))
ggsave("plots/Figure4G.pdf", plot = p, width = 4, height = 4)

## Figure 4H --------------------------------------------------------------------------
# Expression levels of GHs across functional clusters

# Define gene set for GHs
gene_roi <- c("BglS","ManB2","Pgu1","YliI","AmyA","AraH","MalK","RbsB","UgpA","UgpB", "GDB1","GlgB")

# Prepare matrix of avg_log2FC values for selected genes across clusters
matrix <- tmp %>%
  dplyr::filter(name %in% gene_roi) %>%
  arrange(Func) %>%
  complete(name, cluster = celltypes, fill = list(avg_log2FC = 0)) %>%
  dplyr::select(-COG_class, -Func_anno, -Func) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = list(avg_log2FC = 0)) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Reorder matrix rows and columns according to the defined cluster order
matrix <- matrix[gene_roi, celltypes]

# Prepare column annotation (cluster identity)
annotation_col <- data.frame(CellType = colnames(matrix)) 
rownames(annotation_col) <- colnames(matrix)
CellType_colors <- celltype_color[1:length(celltypes)]
names(CellType_colors) <- celltypes

# Define annotation color list for pheatmap
annotation_colors <- list(CellType = CellType_colors)

# Generate color gradient for heatmap
min_value <- min(matrix, na.rm = TRUE)
max_value <- max(matrix, na.rm = TRUE)
gradient_colors <- colorRampPalette(c("darkblue", "white", "salmon", "red"))(100)

# Define breaks for color scaling and ensure zero is included
breaks <- seq(min_value, max_value, length.out = 86)
zero_index <- which.min(abs(breaks - 0)) 
breaks[zero_index] <- 0  

# Plot heatmap without clustering, row/column names, or legend
pdf("plots/Figure4H.pdf", width = 4, height = 7.2)
pheatmap(
  matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_colors = annotation_colors,
  annotation_col = annotation_col,
  color = gradient_colors,
  breaks = breaks,
  annotation_legend = FALSE
)
dev.off()

## Figure 4I --------------------------------------------------------------------------
# Expression levels of methionine and peptide transport-related genes across functional clusters

# Define gene set for methionine biosynthesis and peptide transport
gene_roi <- c(
  "MetH1","MetH2","MetE","GltB2","MetL1","DapA","Asd","COG4870","COG1572",
  "AprE","TldD","COG4412","PepF","OppA","AppF","GltS","DppB","DppD","GdhA",
  "ArcA","ArgF","AguA","TnaA"
)

# Prepare matrix of avg_log2FC values for selected genes across clusters
matrix <- tmp %>%
  dplyr::filter(name %in% gene_roi) %>%  
  arrange(Func) %>%  
  complete(name, cluster = celltypes, fill = list(avg_log2FC = 0)) %>%
  dplyr::select(-COG_class, -Func_anno, -Func) %>%
  distinct() %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = list(avg_log2FC = 0)) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Reorder rows to match the predefined gene order
matrix <- matrix[match(gene_roi, rownames(matrix)), ] 

# Prepare column annotation (cluster identity)
annotation_col <- data.frame(CellType = colnames(matrix)) 
rownames(annotation_col) <- colnames(matrix)
CellType_colors <- celltype_color[1:length(celltypes)]
names(CellType_colors) <- celltypes
matrix <- matrix[gene_roi, celltypes]

# Define annotation color list for pheatmap
annotation_col <- data.frame(CellType = colnames(matrix)) 
rownames(annotation_col) <- colnames(matrix)
CellType_colors <- celltype_color[1:length(celltypes)]
names(CellType_colors) <- celltypes
annotation_colors <- list(CellType = CellType_colors)

# Generate color gradient for heatmap
min_value <- min(matrix, na.rm = TRUE)
max_value <- max(matrix, na.rm = TRUE)
gradient_colors <- colorRampPalette(c("darkblue", "white", "salmon", "darkred"))(100)

# Define breaks for color scaling and ensure zero is included
breaks <- seq(min_value, max_value, length.out = 80.08)
zero_index <- which.min(abs(breaks - 0)) 
breaks[zero_index] <- 0  

# Plot heatmap without clustering, row/column names, or legend
pdf("plots/Figure4I.pdf",width = 5.8,height = 12)
pheatmap(
  matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = F,
  show_colnames = F,
  annotation_colors = annotation_colors,
  annotation_col = annotation_col,
  color = gradient_colors, 
  breaks = breaks,         
  annotation_legend = FALSE
)
dev.off()
