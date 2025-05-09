# ====== Load Required Libraries ======
library(SingleCellExperiment)
library(Matrix.utils)

# ====== Subset for T. denticola and P. gingivalis ======
sub = subset(Combined_seurat_object, species_info %in% c('Treponema denticola', 'Porphyromonas gingivalis'))

# Extract raw counts and metadata
pesudo_counts <- sub@assays$RNA@counts
pesudo_metadata <- sub@meta.data
pesudo_metadata$cluster_id <- factor(rep('same', nrow(pesudo_metadata)))  # Dummy cluster
pesudo_metadata$sample_id <- factor(sub@meta.data$orig.ident)  # Use sample ID

# Create SingleCellExperiment object
sce_bulk <- SingleCellExperiment(
  assay = list(counts = pesudo_counts),
  colData = pesudo_metadata
)

# ====== Pseudo-bulk Aggregation by Species and Sample ======
group <- colData(sce_bulk)[, c("species_info", "sample_id")]
pb <- aggregate.Matrix(t(counts(sce_bulk)), groupings = group, fun = "sum")

# Create sample and species metadata for DESeq2
group <- data.frame(Info = rownames(pb))
group$species_info <- sapply(stringr::str_split(rownames(pb), "_", n = 2), `[`, 1)
group$species_info <- factor(group$species_info, levels = c('Treponema denticola', 'Porphyromonas gingivalis'))
group$sample <- sapply(stringr::str_split(rownames(pb), "_", n = 2), `[`, 2)

# Transpose matrix for DESeq2 input
pb_t <- t(pb)

# ====== Differential Expression Analysis (DESeq2) ======
dds <- DESeqDataSetFromMatrix(pb_t, colData = group, design = ~ species_info)
dds <- DESeq(dds)

# Get DE results: T. denticola vs P. gingivalis
res <- results(dds, 
               contrast = c('species_info', 'Treponema denticola', 'Porphyromonas gingivalis'),
               alpha = 0.05) %>% 
  as.data.frame() %>% 
  arrange(desc(log2FoldChange)) %>%
  mutate(gene = rownames(.))

# ====== Functional Annotation with COG ======
cog.20.def <- read.delim("cog-20.def.tab", header = FALSE)
colnames(cog.20.def) = c('COG_ID','COG_class','Func_anno','gene','Func','time','Anno')
cog.20.def = cog.20.def %>% mutate(gene = ifelse(gene == '', COG_ID, gene))
res = left_join(res, cog.20.def, by = 'gene')

# Export results
write.csv(res, paste0(substr(Sys.Date(), 3, 10), '_Td_Pg_DEG_Pseudolbulk.csv'))

# ====== Volcano Plot with Highlighted Genes ======
highlight_genes <- c("grdB-2.1", "grdE-2", "grdB-2", "grdA1", "grdA", "DdpA", "AppF", "DppD", "DppC", "DppB", 
                     "PotB", "fimA", "mfa2.1", "mfa1", "mfa3", "AguA", "PotD", "PotC", "fimC", "FimA", "eMpr", 
                     "fimE", "COG1913", "OPT", "SpeA", "PepO", "PTR2", "LpxC", "PepC", "PepD2")
highlight_data <- subset(res, gene %in% highlight_genes)

# Generate volcano plot
p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(padj), color = log2FoldChange)) +
  geom_point(data = highlight_data, aes(size = -log10(padj), color = log2FoldChange),
             stroke = 1.0, shape = 21, color = "black") +
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 4)) +
  geom_text(data = highlight_data, aes(label = gene), vjust = -1, hjust = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.justification = c(0, 1)) +
  xlab("Log2FC") + ylab("-Log10(FDR q-value)")

# Save plot
ggsave(filename = paste0(substr(Sys.Date(), 3, 10), "_Pg_Td_volplot.pdf"),
       plot = p, width = 6.4, height = 5, dpi = 600)
