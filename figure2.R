library(Seurat)
library(SCP)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plot1cell)
library(cowplot)
library(Hmisc)

# load the data
seu <- readRDS("data/seurat_object_celltype.rds")
Idents(seu) <- "seurat_clusters"

# run standard pipeline of normalization and scaling
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
seu@meta.data$sample <- seu@meta.data$orig.ident

# set color palettes
sample_color <- c(rep("#3D98D3FF", 8), rep("#96281BFF", 8))
names(sample_color) <- c(paste0("HC", 1:8), paste0("Pd", 1:8))
group_color <- c("HC" = "#3D98D3FF", "Pd" = "#96281BFF")
cluster_color <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
  "#dda0dd", "#8c564b", "#e377c2", "#7f7f7f", 
  "#636363", "#17becf", "#ffd700", "#bcbd22", 
  "#dcdcdc", "#ff9896", "#c49c94", "#ffbb78"
)
tmp <- fread("data/species_color.csv", sep = ",", header = TRUE)
species_color <- setNames(tmp$color, tmp$species)

## Figure 2A -------------------------------------------------------------------
## Determine the optimal dimension and resolution value 
## (Recommend running on the server)

for (dim in c(5, 10, 15, 20, 30)){
  seu <- seu  %>%
    RunUMAP(reduction = "harmony", dims = 1:dim) %>%
    FindNeighbors(reduction = "harmony", dims = 1:dim) %>%
    FindClusters(resolution = c(seq(0.1,1,0.1))) %>%
    identity()
  P <- clustree(seu@meta.data, prefix = "SCT_snn_res.")
  ggsave(P, filename = paste0("plots/Figure2A_Cutree_", dim, ".pdf"), width = 10, height = 10)
  
  P <- DimPlot(seu, reduction = "umap", label = F, pt.size = 0.005, group.by = "seurat_clusters",
               cols = rep("#2C6339", length(unique(seu@meta.data$seurat_clusters))), raster = FALSE) +
    labs(x = "UMAP1", y = "UMAP2",title = "") +
    theme(axis.text.y = element_blank(),
          text = element_text(size = 15, face="bold"),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank())   
  ggsave(P, filename = paste0("plots/Figure2A_UMAP_dims_",dim,".pdf"), width = 5, height = 5)}


## Figure 2B -------------------------------------------------------------------
## Generate the circos plot of the atlas

# only show most abundant species in the plot
seu@meta.data <- seu@meta.data %>%
  dplyr::mutate(
    species = dplyr::case_when(
      species_info %in% tmp$species ~ species_info,
      TRUE ~ "Others"
    ),
    species = factor(species, levels = tmp$species)
  )

# prepare the data for the circos plot
circ_data <- plot1cell::prepare_circlize_data(seu, scale = 0.62) %>%
  dplyr::mutate(
    group = factor(group, levels = c("HC", "Pd")),
    species = factor(species, levels = as.character(tmp$species)),
    cluster_color = cluster_color[match(seurat_clusters, levels(seu$seurat_clusters))],
    seurat_clusters = factor(seurat_clusters, levels = levels(seu$seurat_clusters)),
    Cluster = factor(Cluster, levels = levels(seu$seurat_clusters))
  ) %>% 
  dplyr::left_join(., tmp, by = "species") %>%
  tibble::column_to_rownames("BC")

# set the random seed for reproducibility
set.seed(1234)

#' 
#' Modify the plot_circlize function from the plot1cell package.
#' Original code: https://github.com/TheHumphreysLab/plot1cell
#' 
plot_circlize2 <- function(data_plot, contour.levels = c(0.2,0.3), pt.size = 0.5,
                           kde2d.n = 1000, contour.nlevels = 100, bg.color="#F9F2E4",
                           track_col = NULL, pt_color_column = NULL, 
                           track_name = "Clusters", label.cex = 0.5, track_lwd = 3) {
  # contour trajectory
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  
  # assign colors to the track 
  track_labels <- names(table(data_plot$Cluster))
  track_colors <- alpha(track_col, 0.8) %||% scales::hue_pal()(length(track_labels))
  
  # assign colors to the points
  if ((!is.null(pt_color_column)) && (pt_color_column %in% colnames(data_plot))) {
    data_plot$pcolor <- data_plot[[pt_color_column]]
  } else {
    col_df <- data.frame(Cluster = track_labels, pcolor = track_colors)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order,]
  }
  
  # obtain the cell count proportions
  cluster_counts <- table(data_plot$Cluster)
  total_cells <- nrow(data_plot)
  cluster_props <- as.numeric(cluster_counts / total_cells)
  
  # set the gap parameters
  n_clusters <- length(cluster_counts)
  gap_degrees <- c(rep(2, n_clusters - 1), 18)
  total_gap_degrees <- sum(gap_degrees)
  
  # calculate the total radians available
  total_available_degrees <- 360 - total_gap_degrees
  total_available_radians <- total_available_degrees * pi / 180
  
  # calculate the radian for each cluster
  cluster_radians <- cluster_props * total_available_radians
  
  # recalculate the polar coordinates considering the gap
  data_plot$x_polar2_new <- numeric(nrow(data_plot))
  current_angle <- 0
  for (i in 1:n_clusters) {
    cluster <- names(cluster_counts)[i]
    cluster_idx <- which(data_plot$Cluster == cluster)
    n_cells <- length(cluster_idx)
    
    # evenly distribute cells
    cell_angles <- seq(current_angle, current_angle + cluster_radians[i], 
                       length.out = n_cells + 1)[1:n_cells]
    
    data_plot$x_polar2_new[cluster_idx] <- cell_angles
    
    # move to the start position of the next cluster
    current_angle <- current_angle + cluster_radians[i] + (gap_degrees[i] * pi / 180)
  }
  
  # replace the polar coordinates
  data_plot$x_polar2 <- data_plot$x_polar2_new
  
  # set the circular layout parameters
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0,0,0,0), track.margin = c(0.01,0), "track.height" = 0.01, 
             gap.degree = gap_degrees,
             points.overflow.warning = FALSE)
  
  # plot the track
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter,
                             CELL_META$cell.ylim[2] + mm_y(4),
                             CELL_META$sector.index,
                             cex=0.5, col = "black", facing = "bending.inside", niceFacing = T)
               })
  
  # plot the arc
  for(i in 1:length(track_labels)){
    dd <- data_plot[data_plot$Cluster == track_labels[i],]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, 
                    col = track_colors[i],  lwd = track_lwd, 
                    sector.index = track_labels[i])
  }
  
  # add the track name
  if (!is.null(track_name)) {
    text(x = 1, y = 0.1, labels = track_name, cex = 0.4, col = "black",srt = -90)
  }
  
  # plot the points (generating the UMAP plot)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$pcolor, 0.2), cex = pt.size);
  
  # # add the contour
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, col = "#ae9c76", add =TRUE)
  
  return(data_plot)
}

#' 
#' Modify the add_track function from the plot1cell package.
#' 
add_track2 <- function(data_plot, group, track_num, track_lwd = 3, colors = NULL,
                       track_margin = c(0.01, 0.01)) {
  if (track_num < 2) {
    stop("The first track is the cluster track. Please change the track_num to a value greater than 1")
  }
  
  # add the track margin
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, track.margin = track_margin)
  celltypes <- names(table(data_plot$Cluster))
  group_names <- names(table(data_plot[,group]))
  
  # set the colors
  if(is.null(colors)){
    col_group <- scales::hue_pal()(length(group_names))
    names(col_group) <- group_names
  } else {
    col_group<- colors
  }
  
  for(i in 1:length(celltypes)) {
    cluster <- celltypes[i]
    data_plot_cl <- data_plot[data_plot$Cluster == cluster, ]
    if (nrow(data_plot_cl) == 0) next
    
    # obtain the polar coordinates of the current cluster
    cluster_min <- min(data_plot_cl$x_polar2)
    cluster_max <- max(data_plot_cl$x_polar2)
    cluster_length <- cluster_max - cluster_min
    
    # calculate the cell counts and proportions of each group in the cluster
    group_counts <- table(data_plot_cl[, group])
    if (length(group_counts) == 0) next
    group_props <- as.numeric(group_counts / sum(group_counts))
    
    # calculate the start and end position
    group_starts <- cluster_min + c(0, cumsum(group_props * cluster_length)[-length(group_props)])
    group_ends <- group_starts + group_props * cluster_length
    
    # plot the arc for each group
    for (j in 1:length(group_counts)) {
      group_name <- names(group_counts)[j]
      circos.segments(
        x0 = group_starts[j], 
        y0 = 0, 
        x1 = group_ends[j], 
        y1 = 0, 
        col = col_group[group_name], 
        sector.index = cluster, 
        lwd = track_lwd
      )
    }
  }
}

# generate the circos plot
pdf("plots/Figure2B.pdf", width = 9, height = 9)
circ_data <- plot_circlize2(circ_data, track_col = cluster_color, pt_color_column = "cluster_color",
                            track_name = NULL, track_lwd = 8, pt.size = 0.2, bg.color = "white",
                            kde2d.n = 200, label.cex = 0.6)
add_track2(circ_data, group = "species", colors = species_color, track_num = 2, track_lwd = 8, track_margin = c(0.03, 0.03))
add_track2(circ_data, group = "group", colors = group_color, track_num = 2, track_lwd = 8, track_margin = c(0.015, 0.015))
dev.off()

# save the plot legend
p <- DimPlot(seu, group.by = "species", cols = species_color, raster = FALSE) +
  theme_void() +
  guides(color = guide_legend(override.aes = list(label = "", size = 3))) +
  theme(legend.text = element_text(face = "italic"))
legend_species <- cowplot::get_legend(p)
pdf("plots/Figure2B_species_legend.pdf", width = 6, height = 4)
cowplot::plot_grid(legend_species, ncol = 1, align = "v")
dev.off()

## Figure 2C -------------------------------------------------------------------
## Plot gene count and cell count distribution for each sample

# retrieve the data from the Seurat object
data_plot <- seu@meta.data %>% 
  dplyr::select(sample, nFeature_RNA, nCount_RNA)
sample_levels <- levels(data_plot$sample)

# calculate cell counts and median gene counts
cell_stats <- data_plot %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(cell_n = n(),
            median_gene = median(nFeature_RNA, na.rm = TRUE),
            .groups = "drop")

# join the stats back
data_plot <- left_join(data_plot, cell_stats, by = "sample")

#  generate the box plot for gene and cell counts
p <- ggplot(data_plot, aes(x = sample, y = nFeature_RNA, fill = cell_n)) +
  geom_boxplot() +
  geom_text(data = cell_stats,
            aes(x = sample,
                y = max(data_plot$nFeature_RNA)*1.05,
                label = round(median_gene, 0)),
            inherit.aes = FALSE,
            size = 3, hjust = 0.5) +
  scale_fill_gradient(low = "#d1e5f0", high = "#2166ac") +
  labs(x = "Sample", y = "Gene Count", fill = "Cell Count") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black")) +
  scale_x_discrete(limits = sample_levels)
ggsave("plots/Figure2C.pdf", plot = p, width = 8, height = 4)

# output a table of the average and median gene counts for each sample
sample_stats <- data_plot %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    mean_nFeature_RNA = round(mean(nFeature_RNA, na.rm = TRUE), 3),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE)
  )
fwrite(sample_stats, "data/Figure2C_median_gene_count.csv", sep = ",", 
       row.names = FALSE, col.names = TRUE, quote = FALSE)

## Figure 2D -------------------------------------------------------------------
## Peform the correlation analysis of gene expression across samples

# extract the data for each group
DefaultAssay(seu) <- "RNA"
seu_HC <- subset(seu, group == "HC")
seu_Pd <- subset(seu, group == "Pd")

# compute the pseudobulk expression across all samples for each gene
smean_all <- AggregateExpression(seu, 
                                 group.by = "orig.ident", 
                                 assays = "RNA", 
                                 slot = "data")$RNA %>% as.matrix()
idx <- rowSums(smean_all) > 0
smean_all <- smean_all[idx,]

# compute correlations
scor_all <- Hmisc::rcorr(smean_all, type = "spearman")
fwrite(scor_all$r, file = "data/Figure2D_corr.csv", sep = ",", row.names = T, col.names = T, quote = F)

## visualization
heatmap_color <- RColorBrewer::brewer.pal(name = "RdBu", n = 11)
pal <- rev(colorRampPalette(heatmap_color)(100))
pdf(file = "plots/Figure2D.pdf", width = 6, height = 5.5)
pp <- pheatmap(scor_all$r, color = pal, border_color = NA, cluster_cols = F, 
               show_colnames = T, cluster_rows = F, treeheight_row = 0, treeheight_col = 0)
print(pp)
dev.off()

## Figure 2E -------------------------------------------------------------------
## generate MDS plots based on gene expression and species-level abundance

# obtain the pseudobulk expression value for each gene across all samples
DefaultAssay(seu) <- "RNA"
smean_all <- AggregateExpression(seu, 
                                 group.by = "orig.ident", 
                                 assays = "RNA", 
                                 slot = "data")$RNA %>% as.matrix()
idx <- rowSums(smean_all) > 0
smean_all <- smean_all[idx,]

# compute the MDS
nmds <- vegan::metaMDS(t(smean_all), distance = "bray", k = 2, trymax = 100)
message("Stress: ", nmds$stress)

# run ANOSIM analysis
dist_mat <- vegan::vegdist(t(smean_all), method = "bray")
groups <- c(rep("HC", 8), rep("Pd", 8))
anosim_res <- vegan::anosim(dist_mat, grouping = groups, permutations = 999)
message("ANOSIM R: ", anosim_res$statistic)
message("ANOSIM p: ", anosim_res$signif)

# format the data for plotting
df_nmds <- nmds$points %>%
  as.data.frame() %>% 
  dplyr::mutate(Sample = rownames(.), 
                Group = dplyr::case_when(grepl("HC", Sample) ~ "HC", TRUE ~ "Pd"),
                Sample = factor(Sample, levels = Sample),
                Group = factor(Group, levels = c("HC", "Pd")))
fwrite(df_nmds, file = "data/Figure2E_MDS_expr.csv", sep = ",", row.names = F, col.names = T, quote = F)

# generate the MDS plot
p <- ggplot(df_nmds, aes(x = MDS1, y = MDS2, color = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Group, fill = Group), 
               geom = "polygon", alpha = 0.1, 
               level = 0.95, show.legend = FALSE) +
  geom_text_repel(aes(label = Sample),
                  size = 3,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sample_color, 
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = group_color) +
  labs(title = "", x = "MDS1", y = "MDS2") +
  theme_pubr() +
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.position = "right",
        plot.title = element_blank(),
        legend.text = element_text(size = 10))
ggsave(filename = "plots/Figure2E_left.pdf", plot = p, width = 7, height = 5)

# identify the functional cluster composition across samples
df_species <- seu@meta.data %>% 
  dplyr::select(orig.ident, species_info) %>%
  dplyr::group_by(orig.ident, species_info) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = species_info, values_from = n, values_fill = list(n = 0)) %>%
  tibble::column_to_rownames(var = "orig.ident")
df_species_norm <- sweep(df_species, 1, rowSums(df_species), "/")

# compute the MDS
nmds <- vegan::metaMDS(df_species_norm, distance = "bray", k = 2, trymax = 100)
message("Stress: ", nmds$stress)

# run ANOSIM analysis
dist_mat <- vegan::vegdist(df_species_norm, method = "bray")
groups <- c(rep("HC", 8), rep("Pd", 8))
anosim_res <- vegan::anosim(dist_mat, grouping = groups, permutations = 999)
message("ANOSIM R: ", anosim_res$statistic)
message("ANOSIM p: ", anosim_res$signif)

# format the data for plotting
df_nmds <- nmds$points %>%
  as.data.frame() %>% 
  dplyr::mutate(Sample = rownames(.), 
                Group = dplyr::case_when(grepl("HC", Sample) ~ "HC", TRUE ~ "Pd"),
                Sample = factor(Sample, levels = Sample),
                Group = factor(Group, levels = c("HC", "Pd")))
fwrite(df_nmds, file = "data/Figure2E_MDS_species.csv", sep = ",", row.names = F, col.names = T, quote = F)

# generate the MDS plot
p <- ggplot(df_nmds, aes(x = MDS1, y = MDS2, color = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Group, fill = Group), 
               geom = "polygon", alpha = 0.1, 
               level = 0.95, show.legend = FALSE) +
  geom_text_repel(aes(label = Sample),
                  size = 3,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sample_color, 
                     guide = guide_legend(ncol = 2)) +
  scale_fill_manual(values = group_color) +
  labs(title = "", x = "MDS1", y = "MDS2") +
  theme_pubr() +
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.position = "right",
        plot.title = element_blank(),
        legend.text = element_text(size = 10))
ggsave(filename = "plots/Figure2E_right.pdf", plot = p, width = 7, height = 5)
