# =====================================================================
# Figure 6 — Single-cell transcriptomic analysis of pathogens
# =====================================================================

library(ggplot2)
library(Seurat)
library(dplyr)
library(paletteer)
library(dunn.test)
library(tidyr)
library(tidyverse)
library(hdWGCNA)
library(harmony)
library(ggrepel)
library(tibble)
library(patchwork)
library(SingleCellExperiment)
library(DESeq2)
library(Matrix)
library(scales)
library(ggvenn)
library(ggpubr)
library(stringr)

seu <- readRDS("data/seurat_object_celltype.rds")
DefaultAssay(seu) <- "RNA"
seu <- Seurat::NormalizeData(seu)
seu <- Seurat::FindVariableFeatures(seu)
all.genes <- rownames(seu)
seu <- Seurat::ScaleData(seu, features = all.genes)

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
  "Others" = "Others")
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
  "General metabolic",
  "Proteolysis",
  "Signal transduction",
  "Oxidative stress defense",
  "Anaerobic metabolism",
  "Polysaccharide degradation",
  "Transcriptional regulation"
)  
names(cluster_color_sub) <- order_sub

# ---------------------------------------------------------------------
# UMAP visualization
# ---------------------------------------------------------------------
Roi <- "P.gingivalis"
Species_sub <- readRDS("data/P.gingivalis.rds")
order_sub <- c(
  "General metabolic",
  "Proteolysis",
  "Signal transduction",
  "Oxidative stress defense",
  "Anaerobic metabolism",
  "Polysaccharide degradation",
  "Transcriptional regulation"
) 
names(cluster_color_sub) <- order_sub

# Extract UMAP embeddings and metadata
Umap <- Species_sub@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(
    cellType = Species_sub@meta.data$celltype_sub,
    sample = Species_sub@meta.data$orig.ident,
    group = Species_sub@meta.data$group
  ) %>%
  mutate(cellType = factor(cellType, levels = order_sub))

# ---------------------------------------------------------------------
# Figure6A
# ---------------------------------------------------------------------

p <- ggplot(Umap, aes(x = umap_1, y = umap_2, fill = cellType, color = cellType)) +
  tidydr::theme_dr(xlength = 0.2,
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  stat_ellipse(aes(x = umap_1, y = umap_2, fill = cellType),
               geom = "polygon", 
               linetype = 2,    
               alpha = 0.05,      
               show.legend = FALSE, 
               level = 0.96) +    
  geom_point(size = 0.1, show.legend = T) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = 2, hjust = 0.03),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
  ) +
  scale_color_manual(values = cluster_color_sub) +
  scale_fill_manual(values = cluster_color_sub) +
  ggtitle(paste0(nrow(Umap)," cells")) +
  guides(
    fill = guide_legend(title = NULL),
    color = guide_legend(title = NULL, override.aes = list(size = 3))
  )
ggsave("plots/Figure6A.pdf", plot = p, width = 8, height = 4)

## ---------------------------------------------------------------------
## Figure 6B – WGCNA analysis and module UMAP for P. gingivalis
## Step 1: data preparation (Recommended to run on server)
## ---------------------------------------------------------------------

## 1 Load single-cell Seurat object of P. gingivalis
seu_sub <- readRDS("data/P.gingivalis.rds")  
Idents(seu_sub) <- "celltype_sub"                 
set.seed(12345)                                 

## 2: Prepare object for WGCNA analysis
# Select genes expressed in at least 0.5% of cells for downstream WGCNA
seurat_obj <- SetupForWGCNA(
  seu_sub,
  gene_select = "fraction", 
  fraction = 0.005, 
  wgcna_name = "tutorial"
)

## 3: Construct metacells by cluster
# Aggregate similar cells into metacells to reduce noise and improve co-expression analysis
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype_sub"), 
  k = 25,
  max_shared = 15,
  ident.group = "celltype_sub"
)

## 4: Normalize and scale metacell expression
seurat_obj <- NormalizeMetacells(seurat_obj)           
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))  

## 5: Dimensionality reduction
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj)) 
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars = "orig.ident")         
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction = "harmony", dims = 1:10)        

## 6: Set up expression matrix for WGCNA modules
# Group clusters into functional categories for module analysis
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("General metabolic","Proteolysis","Oxidative stress defense",
                 "Signal transduction","Anaerobic metabolism",
                 "Polysaccharide degradation","Transcriptional regulation"),
  group.by="celltype_sub", 
  assay = "SCT", 
  slot = "data" 
)

## 7: Determine soft-thresholding power for network construction
seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
power_table <- GetPowerTable(seurat_obj)

# Function to select optimal soft power where scale-free topology fit R^2 >= 0.8
find_soft_threshold <- function(power_table, threshold = 0.8) {
  valid_rows <- power_table[power_table$SFT.R.sq >= threshold, ]
  if (nrow(valid_rows) == 0) {
    stop("No results found that meet the condition SFT.R.sq >= ", threshold)
  }
  optimal_power <- valid_rows[1, "Power"]
  return(optimal_power)
}
soft_threshold <- find_soft_threshold(power_table, threshold = 0.8)

## 8: Construct signed co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power = soft_threshold, 
  setDatExpr=FALSE, overwrite_tom = TRUE,
  minModuleSize = 25, deepSplit = 1,
  detectCutHeight = 0.995, pamRespectsDendro = FALSE,
  tom_name = "PD"
)

## 9: Compute module eigengenes and connectivity
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "orig.ident"
) 
hMEs <- GetMEs(seurat_obj)
write.csv(hMEs, "data/Figure6_WCGNA_feature_pg.csv")

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "celltype_sub",
  group_name = c("General metabolic","Proteolysis","Oxidative stress defense",
                 "Signal transduction","Anaerobic metabolism",
                 "Polysaccharide degradation","Transcriptional regulation")
)

##10: UMAP embedding of module genes and hub genes
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=15,
  min_dist=0.1
)
umap_df <- GetModuleUMAP(seurat_obj)

##11: Prepare node and edge tables for network visualization
g <- ModuleUMAPPlot(seurat_obj, return_graph = TRUE)
umap_df <- data.frame(
  gene = igraph::V(g)$name,  
  UMAP1 = igraph::V(g)$UMAP1, 
  UMAP2 = igraph::V(g)$UMAP2, 
  module = igraph::V(g)$module,  
  kME = igraph::V(g)$kME,   
  hub = igraph::V(g)$hub    
)
edges_df <- data.frame(
  from = igraph::ends(g, igraph::E(g))[, 1],  
  to = igraph::ends(g, igraph::E(g))[, 2]    
)
edges_df <- merge(edges_df, umap_df, by.x = "from", by.y = "gene", all.x = TRUE)
edges_df <- merge(edges_df, umap_df, by.x = "to", by.y = "gene", suffixes = c(".from", ".to"))

## 12: Save processed node and edge tables for plotting
write.csv(umap_df, "data/Figure6_umap_nodes_pg.csv", row.names = FALSE)
write.csv(edges_df, "data/Figure6_umap_edges_pg.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## Step2: Draw Plot
## ---------------------------------------------------------------------
## Figure6B
color_map <- c("brown" = "ME1",  "blue" = "ME2", "turquoise" = "ME3",
               "grey" = "grey", "yellow"= "ME4", "green"="ME5")
Feature <- read.csv("data/Figure6_WCGNA_feature_pg.csv")
colnames(Feature) <- color_map[colnames(Feature)]
colnames(Feature)[1] <- "BC"
Feature <- Feature[, !colnames(Feature) %in% "grey"]
Species_sub@meta.data <- Species_sub@meta.data %>%
  left_join(.,  Feature, by = "BC") %>%
  mutate(BC2 = BC) %>%
  tibble::column_to_rownames("BC2")
mods <- colnames(Feature)
mods <- mods[!mods %in% c("grey", "BC")]

p <- DotPlot(Species_sub, features= paste0("ME", 1:5), group.by = "celltype_sub", dot.scale = 10) + 
  scale_color_gradientn(colours = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 5)) +
  labs(x = "Cluster", y = "Module", title = "Module activity by clusters") +
  theme(
    axis.text.x = element_text(color = "black"),
    text = element_text(size = 12),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")
  )
ggsave("plots/Figure6B.pdf", plot = p, width = 8.6,height = 6)

## Figure6C
MEs <- Feature %>% column_to_rownames("BC")    
MEs <- MEs[,colnames(MEs) != "grey"]         
MEs[MEs < 0] <- 0                             
MEs <- t(MEs)                               
ME_assay <- Seurat::CreateAssayObject(MEs)  
BCs <- colnames(MEs)
orig.ident <- sub("_.*", "", BCs)
Group <- case_when(
  grepl("Pd", orig.ident) ~ "PD",
  grepl("HC", orig.ident) ~ "HC",
  TRUE ~ "Unknown"
)
tmp <- data.frame(BCs, orig.ident, Group)

group1 <- tmp %>% dplyr::filter(Group == "PD") %>% pull(BCs)  # Barcodes for PD group
group2 <- tmp %>% dplyr::filter(Group == "HC") %>% pull(BCs)  # Barcodes for HC group

# Compare module activity between PD and HC using Wilcoxon test
DMEs <- FindMarkers(
  ME_assay, 
  cells.1 = group1, 
  cells.2 = group2, 
  slot = "counts", 
  test.use = "wilcox", 
  only.pos = F, 
  logfc.threshold = 0,       # Keep all modules regardless of fold change
  min.pct = 0.15             # Modules expressed in at least 15% of cells
)
DMEs$module <- rownames(DMEs) # Save module names

p <- ggplot(DMEs, aes(y = module)) +
  geom_segment(                            
    aes(x = avg_log2FC, xend = 0, yend = module),
    color = "black"
  ) +
  geom_point(                            
    aes(
      x = avg_log2FC,                      
      size = abs(avg_log2FC),               
      color = -log10(p_val)                 
    )
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_viridis_c(option = "C") +   
  labs(
    x = "Avg. Log2(FC)", 
    y = "", 
    size = "Log2|FC|", 
    color = "-log10(p-value)"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold", hjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  )
ggsave("plots/Figure6C.pdf", plot = p, width = 6, height = 4)


# ---------------------------------------------------------------------
# Figure 6D – Ratio of expressed peptidase/protease genes across clusters
# ---------------------------------------------------------------------

# 1. Retrieve raw RNA counts from the Seurat object for all cells
raw_matrix <- GetAssayData(Species_sub, assay = "RNA", slot = "counts")

# 2. Load list of genes involved in amino acid transport and metabolism
Interested_gene <- as.matrix(read.csv("data/Figure6_genes_in_AA_transport_metabolism.csv", header=T))
unique(Interested_gene[, 5])  # Inspect functional categories in the gene list

# 3. Filter genes annotated as "peptidase" or "protease" (columns 4 = gene symbol)
filtered_genes <- subset(Interested_gene, (Interested_gene[, 5] == "peptidase") | (Interested_gene[, 5] == "protease"), select = 4)

# 4. Keep only genes actually expressed in our dataset
expression_matrix <- GetAssayData(Species_sub, assay = "RNA", slot = "counts")
filtered_genes <- filtered_genes[filtered_genes %in% rownames(expression_matrix)]

# 5. Calculate number of filtered genes expressed in each cell
filtered_genes_expressed <- apply(expression_matrix[filtered_genes, ], 2, function(x) sum(x > 0))

# 6. Calculate total number of genes expressed in each cell
all_genes_expressed <- apply(expression_matrix, 2, function(x) sum(x > 0))

# 7. Compute ratio: (# peptidase/protease genes expressed) / (total genes expressed)
filtered_genes_ratio <- filtered_genes_expressed / all_genes_expressed

# 8. Add this ratio as metadata to the Seurat object
Species_sub <- AddMetaData(Species_sub, metadata = filtered_genes_ratio, col.name = "DegradationGenesRatio")
head(Species_sub@meta.data$DegradationGenesRatio)  # Check result

# 9. Prepare data for plotting: select cluster and ratio, remove NA
degradation_data <- Species_sub@meta.data %>%
  select(celltype_sub, DegradationGenesRatio) %>%
  drop_na() 

# 10. Compute mean ratio per cluster and sort for plotting
celltype_mean <- degradation_data %>%
  group_by(celltype_sub) %>%
  summarise(Mean_DegradationGenesRatio = mean(DegradationGenesRatio, na.rm = TRUE)) %>%
  arrange(desc(Mean_DegradationGenesRatio))
print(celltype_mean)

# 12. Prepare for plotting
order_sub2 <- setNames(paste0("C", 0:(length(order_sub)-1)), order_sub)
order_sub3 <- setNames(c("b", "a", "d", "c", "c", "c", "c"), order_sub)
degradation_data <- degradation_data %>%
  mutate(DegradationGenesPCT = DegradationGenesRatio * 100,
         cluster = order_sub2[celltype_sub],
         label = order_sub3[celltype_sub],
         celltype_sub = factor(celltype_sub, levels = celltype_mean$celltype_sub))
text_data <- degradation_data %>%
  # mutate(celltype_sub = as.character(celltype_sub)) %>%
  group_by(celltype_sub) %>%
  summarise(
    mean_pct = mean(DegradationGenesPCT),
    label = unique(order_sub3[as.character(celltype_sub)])
  )

# 13. Plot mean Degradation Genes Ratio per cluster with error bars
p <- ggplot(degradation_data, aes(x = celltype_sub, y = DegradationGenesPCT, fill = celltype_sub)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7, color = "black") +  
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, size = 0.7) + 
  geom_text(
    data = text_data,
    aes(x = celltype_sub, y = mean_pct, label = label),
    vjust = -1.3, 
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(labels = order_sub2[levels(degradation_data$celltype_sub)]) + 
  scale_y_continuous(
    breaks = seq(0, 5, by = 1), 
    labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0")
  ) +
  labs(title = "Genes involved in\n protease/peptidase", 
       x = "",  
       y = "GFP of the reaction (%)",
       fill = "Cluster") +  
  scale_fill_manual(values = cluster_color_sub) +  
  theme_classic(base_size = 16) + 
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14), 
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,),  
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0.5), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = margin(1.5, 1, 1, 1, "cm")  
  )
ggsave("plots/Figure6D.pdf", plot = p, width = 10, height = 8)

# 13. Statistical comparison across clusters using Kruskal-Wallis test
kruskal_result <- kruskal.test(DegradationGenesRatio ~ celltype_sub, data = degradation_data)
print(kruskal_result)

# 14. Post-hoc pairwise comparison with Dunn"s test (Bonferroni correction)
dunn_result <- dunn.test(degradation_data$DegradationGenesRatio, 
                         g = degradation_data$celltype_sub, 
                         method = "bonferroni")
write.csv(dunn_result, file = "data/Figure6D_dunn_test_Pi.csv", row.names = FALSE)

# 15. Generate summary table of Dunn"s test
dunn_summary <- data.frame(
  Comparison = dunn_result$comparisons,
  Z.value = dunn_result$Z, 
  P.value = dunn_result$P,  
  stringsAsFactors = FALSE
)
write.csv(dunn_summary, file = "data/Figure6D_dunn_test_summary_Pi.csv", row.names = FALSE)

##Figure6E
cog.20.def <- read.delim("data/cog-20.def.tab", header=FALSE)
colnames(cog.20.def) <- c("gene","COG_class","Func_anno","name","Func","time","Anno")
cog.20.def <- cog.20.def %>% mutate(name = ifelse(name == "",gene,name))
umap_nodes <- read.csv("data/Figure6_umap_nodes_pg.csv")
umap_edges <- read.csv("data/Figure6_umap_edges_pg.csv")
umap_nodes$module <- color_map[umap_nodes$module]
umap_edges$module.from <- color_map[umap_edges$module.from]
umap_edges$module.to <- color_map[umap_edges$module.to]
top_ME3_genes <- umap_nodes[umap_nodes$module == "ME3" & umap_nodes$kME > 0.5, ]
top_ME3_genes$name <- top_ME3_genes$gene
top_ME3_genes_anno <- top_ME3_genes %>%
  left_join(cog.20.def, by = "name") %>%
  dplyr::filter(!is.na(COG_class))
top_ME3_genes_anno <- top_ME3_genes_anno %>%
  mutate(COG_class_label = cog_labels[COG_class]) %>%
  mutate(COG_class_label = if_else(nchar(COG_class) == 1, COG_class_label, "Others"))

p <- ggplot(top_ME3_genes_anno, aes(x = reorder(name, kME), y = kME, fill = COG_class_label)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  coord_flip() +
  labs(x = "", y = "kME Value", title = "ME3 module genes (kME > 0.5)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", face = "italic", size = 12),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 11)
  ) +
  scale_fill_manual(values = cog_label_color,
                    labels = function(x) stringr::str_wrap(x, width = 15)) 
ggsave("plots/Figure6E.pdf", plot = p, width = 5, height = 5)

## ---------------------------------------------------------------------
## Figure 6G – UMAP visualization of P.intermedia subclusters
## ---------------------------------------------------------------------

Roi <- "P.intermedia"
Species_sub <- readRDS("data/P.intermedia.rds")   # Load Seurat object for P. intermedia
order_sub <- c(
  "General metabolic",
  "Signal transduction",
  "Anaerobic metabolism",
  "Protein degradation",
  "Glycosyl transfer",
  "Sugar lactone hydrolysis",
  "Transcriptional regulator",
  "Formylglycine generation",
  "Chitin hydrolysis",
  "Fatty acid synthesis"
) 
names(cluster_color_sub) <- order_sub

# Extract UMAP embeddings and metadata
Umap <- Species_sub@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cellType = Species_sub@meta.data$celltype_sub, 
        sample = Species_sub@meta.data$orig.ident,     
        group = Species_sub@meta.data$group) %>%     
  mutate(cellType = factor(cellType, levels = order_sub)) 

# Plot UMAP with ellipses and colored by subcluster
p <- ggplot(Umap, aes(x = umap_1, y = umap_2, fill = cellType, color = cellType)) +
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  stat_ellipse(aes(x = umap_1, y = umap_2, fill = cellType),
               geom = "polygon", 
               linetype = 2, 
               alpha = 0.05, 
               show.legend = FALSE, 
               level = 0.96) +          
  geom_point(size = 0.1, show.legend = T) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(face = 2, hjust = 0.03),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
    legend.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = cluster_color_sub) + 
  scale_fill_manual(values = cluster_color_sub) +
  ggtitle(paste0(nrow(Umap)," cells")) +  
  guides(
    fill = guide_legend(title = NULL),
    color = guide_legend(title = NULL, override.aes = list(size = 3))
  )
ggsave("plots/Figure6G.pdf", plot = p, width = 8, height = 4)

## ---------------------------------------------------------------------
## Figure 6H – Ratio of degradation-related genes per subcluster
## ---------------------------------------------------------------------

raw_matrix <- GetAssayData(Species_sub, assay = "RNA", slot = "counts") # Raw counts

# Load list of genes involved in amino acid transport/metabolism
Interested_gene <- as.matrix(read.csv("data/Figure6_genes_in_AA_transport_metabolism.csv", header=T))
unique(Interested_gene[, 5])

# Select genes classified as "peptidase" or "protease"
filtered_genes <- subset(Interested_gene, (Interested_gene[, 5] == "peptidase") | (Interested_gene[, 5] == "protease"), select = 4)

# Filter genes present in dataset
expression_matrix <- GetAssayData(Species_sub, assay = "RNA", slot = "counts")
filtered_genes <- filtered_genes[filtered_genes %in% rownames(expression_matrix)]

# Compute number of expressed genes per cell
filtered_genes_expressed <- apply(expression_matrix[filtered_genes, ], 2, function(x) sum(x > 0))
all_genes_expressed <- apply(expression_matrix, 2, function(x) sum(x > 0))

# Compute ratio of expressed degradation genes per cell
filtered_genes_ratio <- filtered_genes_expressed / all_genes_expressed

# Add ratio to Seurat metadata
Species_sub <- AddMetaData(Species_sub, metadata = filtered_genes_ratio, col.name = "DegradationGenesRatio")
head(Species_sub@meta.data$DegradationGenesRatio)

# Prepare data for plotting
degradation_data <- Species_sub@meta.data %>%
  select(celltype_sub, DegradationGenesRatio) %>%
  drop_na() 

# Compute mean ratio per subcluster
celltype_mean <- degradation_data %>%
  group_by(celltype_sub) %>%
  summarise(Mean_DegradationGenesRatio = mean(DegradationGenesRatio, na.rm = TRUE)) %>%
  arrange(desc(Mean_DegradationGenesRatio))
print(celltype_mean)

# Set factor levels for plotting order
order_sub2 <- setNames(paste0("C", 0:(length(order_sub)-1)), order_sub)
order_sub3 <- setNames(c("b", "c", "b", "a", rep("b", 4), "d", "b"), order_sub)
degradation_data <- degradation_data %>%
  mutate(DegradationGenesPCT = DegradationGenesRatio * 100,
         cluster = order_sub2[celltype_sub],
         label = order_sub3[celltype_sub],
         celltype_sub = factor(celltype_sub, levels = celltype_mean$celltype_sub))
text_data <- degradation_data %>%
  group_by(celltype_sub) %>%
  summarise(
    mean_pct = mean(DegradationGenesPCT),
    label = unique(order_sub3[as.character(celltype_sub)])
  )

p <- ggplot(degradation_data, aes(x = celltype_sub, y = DegradationGenesPCT, fill = celltype_sub)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7, color = "black") +  
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, size = 0.7) + 
  geom_text(
    data = text_data,
    aes(x = celltype_sub, y = mean_pct, label = label),
    vjust = -1, 
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(labels = order_sub2[levels(degradation_data$celltype_sub)]) + 
  labs(title = "Genes involved in\n protease/peptidase", 
       x = "",  
       y = "GFP of the reaction (%)",
       fill = "Cluster") +  
  scale_fill_manual(values = cluster_color_sub) +  
  theme_classic(base_size = 16) + 
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14), 
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,),  
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0.5), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = margin(1.5, 1, 1, 1, "cm")  
  )
ggsave("plots/Figure6H.pdf", plot = p, width = 7, height = 7)

# Statistical comparison between subclusters
kruskal_result <- kruskal.test(DegradationGenesRatio ~ celltype_sub, data = degradation_data)
print(kruskal_result)

# Post-hoc pairwise comparison using Dunn"s test with Bonferroni correction
dunn_result <- dunn.test(degradation_data$DegradationGenesRatio, 
                         g = degradation_data$celltype_sub, 
                         method = "bonferroni")
write.csv(dunn_result, file = "data/Figure6G_dunn_test_Pi.csv", row.names = FALSE)

# Save summary table
dunn_summary <- data.frame(
  Comparison = dunn_result$comparisons,
  Z.value = dunn_result$Z, 
  P.value = dunn_result$P,  
  stringsAsFactors = FALSE
)
write.csv(dunn_summary, file = "data/Figure6G_dunn_test_summary_Pi.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## Figure 6I – Proportion of subclusters across samples
## ---------------------------------------------------------------------

# Calculate counts of each subcluster in each sample
consist <- FetchData(Species_sub, vars = c("orig.ident", "celltype_sub")) %>%
  group_by(orig.ident, celltype_sub) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  complete(orig.ident, celltype_sub, fill = list(count = 0)) %>%  # Fill missing combinations with 0
  group_by(orig.ident) %>%
  mutate(total_count = sum(count)) %>%  # Compute total cells per sample
  ungroup() %>%
  mutate(n = count / total_count,       # Compute proportion per subcluster
        group = case_when( grepl("HC", orig.ident) ~ "HC", TRUE ~ "Pd"))

# Perform Wilcoxon test across groups for each subcluster
results <- consist %>%
  group_by(celltype_sub) %>%
  dplyr::summarise(p_value = perform_wilcox_test(cur_data())) %>%  # Custom function for Wilcoxon
  ungroup() %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "bonferroni")) # Adjust p-values

# Extract significant clusters
significant_clusters <- results %>%
  dplyr::filter(p_value < 0.05) %>%
  pull(celltype_sub)

# Prepare metadata for plotting stacked barplot
consist_plot_data <- Species_sub@meta.data %>%
  dplyr::select(celltype_sub, group) %>%
  mutate(celltype_sub = factor(celltype_sub, levels = order_sub))

# Create stacked barplot
P <- ggplot(consist_plot_data, aes(x = as.factor(group), fill = celltype_sub)) +
  geom_bar(position = "fill") +                      # Show proportions
  scale_y_continuous(labels = scales::percent) +    # Percent scale
  scale_fill_manual(values = cluster_color_sub) +   # Subcluster colors
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", face = "bold"),
        text = element_text(size = 12)) +
  labs(x = "", y = paste0("Proportion of ", Roi, " cells(%)"))

# Compute midpoints for annotation of significant clusters
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

# Annotate significant clusters with asterisk
for (cluster in significant_clusters) {
  cluster_data <- proportions %>% filter(celltype_sub == cluster)
  max_group <- cluster_data$group[which.max(cluster_data$proportion)]
  mid_position <- cluster_data %>% filter(group == max_group) %>% pull(mid_point)
  
  P <- P + annotate("text",
                    x = max_group,
                    y = mid_position,
                    label = "*",
                    size = 5,
                    color = "black",
                    vjust = 0.75)
}
ggsave("plots/Figure6I.pdf", plot = P, width = 3.5, height = 5)


## ---------------------------------------------------------------------
## Figure 6J – Pseudobulk differential expression between Td and Pg
## ---------------------------------------------------------------------

# Subset the main Seurat object for the two species of interest
sub <- subset(seu, species_info %in% c("Treponema denticola","Porphyromonas gingivalis"))

# Extract raw counts and metadata
pesudo_counts <- sub@assays$RNA@counts
pesudo_metadata <- sub@meta.data
pesudo_metadata$cluster_id <- factor(rep("same", nrow(pesudo_metadata))) # Single cluster
pesudo_metadata$sample_id <- factor(sub@meta.data$orig.ident)

# Create SingleCellExperiment object for aggregation
sce_bulk <- SingleCellExperiment(assay = list(counts = pesudo_counts),
                                 colData = pesudo_metadata)

# Define aggregation groups by species and sample
group <- colData(sce_bulk)[,c("species_info","sample_id")]
group_factor <- interaction(group$species_info, group$sample_id, sep = "_", drop = TRUE)
design_matrix <- sparse.model.matrix(~ 0 + group_factor)
colnames(design_matrix) <- levels(group_factor)
pb_aggregated <- counts(sce_bulk) %*% design_matrix
pb <- t(pb_aggregated)
rownames(pb) <- levels(group_factor)
colnames(pb) <- rownames(sce_bulk)

# Prepare group metadata for DESeq2
group <- data.frame(Info = rownames(pb))
group$species_info <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), `[`, 1)
group$species_info <- factor(group$species_info, levels = c("Treponema denticola","Porphyromonas gingivalis"))
group$sample <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), `[`, 2)

pb_t <- t(pb)  # Transpose matrix for DESeq2 (genes x samples)

# Create DESeq2 dataset and run differential expression
dds <- DESeqDataSetFromMatrix(pb_t, colData = group, design = ~ species_info)
dds <- DESeq(dds)

# Extract DE results for Td vs Pg
res <- results(dds, 
               contrast = c("species_info","Treponema denticola","Porphyromonas gingivalis"),
               alpha = 0.05) %>% 
  as.data.frame() %>% arrange(desc(log2FoldChange))
res$gene <- rownames(res)

# Annotate genes with COG functional categories
cog.20.def <- read.delim("data/cog-20.def.tab", header=FALSE)
colnames(cog.20.def) <- c("COG_ID","COG_class","Func_anno","gene","Func","time","Anno")
cog.20.def <- cog.20.def %>% mutate(gene = ifelse(gene == "",COG_ID,gene))
res <- left_join(res, cog.20.def, by = "gene")

# Highlight selected genes in volcano plot
highlight_genes <- c(
  "grdB-2.1", "grdE-2","grdB-2","grdA1","grdA","DdpA","AppF",
  "DppD","DppC","DppB","PotB","fimA","mfa1","mfa3","AguA",
  "PotD","PotC","eMpr","COG1913","OPT","SpeA","PepO","PTR2","PepC","PepD2",
  "FadL","COG4874","COG3291","COG4412","ArgD","FadL","ArcA","ArgF","ArcC"
)
highlight_data <- subset(res, gene %in% highlight_genes)

# Volcano plot
p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(padj), color = log2FoldChange)) +
  geom_point(data = highlight_data, aes(size = -log10(padj), color = log2FoldChange), 
             stroke = 1.0, shape = 21, color = "black") +
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 4)) +
  geom_text_repel(
    data = highlight_data, aes(label = gene), size = 3.5, box.padding = 0.5, 
    point.padding = 0.3, segment.color = "grey50", segment.size = 0.5,  
    segment.alpha = 0.6, min.segment.length = 1, force = 1, force_pull = 1, 
    max.overlaps = 20, max.time = 1, max.iter = 10000, nudge_x = 0.1, 
    nudge_y = 0.1, direction = "both", seed = 123
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",  
        legend.justification = c(0, 1)) +
  xlab("Log2(Fold-change)") +
  ylab("-Log10(p.adj)")
ggsave("plots/Figure6J.pdf", plot = p, width = 6.4, height = 5, dpi = 600)


# ============================================================
# Figure 6K: Gene expression comparison between T. denticola and P. gingivalis
# ============================================================

DefaultAssay(seu) <- "RNA"

species1 <- "Treponema denticola"
species2 <- "Porphyromonas gingivalis"

gene_list <- c("DdpA", "AppF", "grdA1", "grdE-2", 
               "OPT",  "PTR2", "SpeA", "AguA")

colors <- c("Treponema denticola" = "#D55E00", "Porphyromonas gingivalis" = "#0072B2")

# -------------------------------
# Loop through each gene
# -------------------------------
for (gene in gene_list) {
  
  # Extract cells for each species
  species1_cells <- WhichCells(seu, expression = species_info == species1)
  species2_cells <- WhichCells(seu, expression = species_info == species2)
  
  # Extract expression values
  species1_expr <- GetAssayData(seu, slot = "data")[gene, species1_cells]
  species2_expr <- GetAssayData(seu, slot = "data")[gene, species2_cells]
  
  # Create combined data frame
  expression_data <- data.frame(
    Expression = c(as.vector(species1_expr), as.vector(species2_expr)),
    Species = factor(rep(c(species1, species2),
                         times = c(length(species1_expr), length(species2_expr)))),
    Sample = c(seu@meta.data[species1_cells, "orig.ident"],
               seu@meta.data[species2_cells, "orig.ident"])
  )
  
  # Wilcoxon test
  wilcox_test_result <- wilcox.test(species1_expr, species2_expr, alternative = "two.sided")
  
  # Sample-wise mean expression
  sample_means <- expression_data %>%
    group_by(Species, Sample) %>%
    summarise(mean_expr = mean(Expression), .groups = "drop")
  
  # -------------------------------
  # Plot
  # -------------------------------
  p <- ggplot(expression_data, aes(x = Species, y = Expression, fill = Species)) +
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
    )  +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
    scale_x_discrete(labels = c("P.g", "T.d")) +
    scale_fill_manual(values = colors) +

    labs(title = substr(gene, 1, 4), y = "Expression Level", x = "") +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c(species1, species2)),
      label = "p.signif",
      vjust = 0.5,
      tip.length = 0,
      size = 6,
      label.y = max(expression_data$Expression) * 1.1,
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns"))
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(size = 14, face = "italic"),
      plot.title = element_text(size = 18, face = "italic", hjust = 0.5),
      axis.ticks.length = unit(0.24, "cm")
    )
  ggsave(paste0("plots/Figure6K_", substr(gene,1,4), ".pdf"), width = 5, height = 5)
  
  # Print Wilcoxon result
  print(gene)
  print(wilcox_test_result)
}
