
# CRC  GSE144259, GSE50760, GSE87096
# step1_Libraries 
library(data.table)
library(tidyverse)
library(DESeq2)
library(sva)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(BiocParallel)
library(grid)
library(scales)

options(stringsAsFactors = FALSE)
register(SerialParam())

# step2_Working directory folder structure 
setwd("C:/Users/akhil/comet_downlaod/colorectal_bulk_rna/data")

dir.create("results_CRC/01_objects",                   recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/02_tables/metadata",           recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/02_tables/counts",             recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/02_tables/DEG",                recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/02_tables/enrichment",         recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/QC_PCA",              recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/DEG",                 recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/Heatmap",             recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/enrichment/png",      recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/enrichment/svg",      recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/enrichment/tiff",     recursive = TRUE, showWarnings = FALSE)

# step3_Build Metadata manually 

# GSE144259: 3T + 3N | EXCLUDE 3 metastasis
meta_144259 <- data.frame(
  GSM_ID     = c("GSM4284531","GSM4284532",
                 "GSM4284534","GSM4284535",
                 "GSM4284537","GSM4284538"),
  SampleName = c("CRC1N","CRC1T","CRC2N","CRC2T","CRC3N","CRC3T"),
  Condition  = c("Normal","Tumor","Normal","Tumor","Normal","Tumor"),
  GSE        = "GSE144259",
  stringsAsFactors = FALSE
)

#  GSE50760: 18T + 18N | EXCLUDE 18 metastasis (suffix -3) ---
gsm_50760_tumor  <- paste0("GSM", 1228184:1228201)   # primary CRC, suffix -1
gsm_50760_normal <- paste0("GSM", 1228202:1228219)   # normal colon, suffix -2
# GSM1228220:1228237 = metastasis → NOT included

meta_50760 <- data.frame(
  GSM_ID     = c(gsm_50760_tumor, gsm_50760_normal),
  SampleName = c(
    paste0("AMC_", c(2,3,5,6,7,8,9,10,12,13,17,18,19,20,21,22,23,24), "-T"),
    paste0("AMC_", c(2,3,5,6,7,8,9,10,12,13,17,18,19,20,21,22,23,24), "-N")
  ),
  Condition  = c(rep("Tumor", 18), rep("Normal", 18)),
  GSE        = "GSE50760",
  stringsAsFactors = FALSE
)

# --- GSE87096: RNA-seq ONLY — 6T + 6N | EXCLUDE 24 MeDIP/hMeDIP samples ---
meta_87096 <- data.frame(
  GSM_ID     = c(
    # Normal RNA-seq
    "GSM2322216","GSM2322222","GSM2322228",
    "GSM2322234","GSM2322240","GSM2322246",
    # Tumor RNA-seq
    "GSM2322219","GSM2322225","GSM2322231",
    "GSM2322237","GSM2322243","GSM2322249"
  ),
  SampleName = c(
    "2504588_N","2512618_N","2539382_N","2551349_N","2553763_N","S12_N",
    "2504588_T","2512618_T","2539382_T","2551349_T","2553763_T","S12_T"
  ),
  Condition  = c(rep("Normal", 6), rep("Tumor", 6)),
  GSE        = "GSE87096",
  stringsAsFactors = FALSE
)

#  Combine all metadata ---
final_metadata             <- bind_rows(meta_144259, meta_50760, meta_87096)
final_metadata$condition   <- factor(final_metadata$Condition, levels = c("Normal", "Tumor"))
final_metadata$batch       <- factor(final_metadata$GSE)
rownames(final_metadata)   <- final_metadata$GSM_ID

# Save combined metadata
write.csv(final_metadata,
          "results_CRC/02_tables/metadata/CRC_combined_metadata.csv",
          row.names = FALSE)

message(" Metadata built:")
message("Total samples: ", nrow(final_metadata))
print(table(final_metadata$GSE, final_metadata$Condition))

# step4_Load count matrices 
load_counts <- function(file) {
  df           <- fread(file, header = TRUE) %>% as.data.frame()
  rownames(df) <- df[[1]]
  df[[1]]      <- NULL
  as.matrix(df)
}

message("Loading count matrices...")
mat_144259 <- load_counts("GSE144259_raw_counts_GRCh38.p13_NCBI.tsv.gz")
mat_50760  <- load_counts("GSE50760_raw_counts_GRCh38.p13_NCBI.tsv.gz")
mat_87096  <- load_counts("GSE87096_raw_counts_GRCh38.p13_NCBI.tsv.gz")

message("GSE144259: ", nrow(mat_144259), " genes x ", ncol(mat_144259), " samples")
message("GSE50760:  ", nrow(mat_50760),  " genes x ", ncol(mat_50760),  " samples")
message("GSE87096:  ", nrow(mat_87096),  " genes x ", ncol(mat_87096),  " samples")

# step5_Filter each matrix to keep ONLY selected samples 
mat_144259 <- mat_144259[, intersect(colnames(mat_144259), meta_144259$GSM_ID), drop = FALSE]
mat_50760  <- mat_50760[,  intersect(colnames(mat_50760),  meta_50760$GSM_ID),  drop = FALSE]
mat_87096  <- mat_87096[,  intersect(colnames(mat_87096),  meta_87096$GSM_ID),  drop = FALSE]

message("After sample filtering:")
message("GSE144259: ", ncol(mat_144259), " | GSE50760: ", ncol(mat_50760), " | GSE87096: ", ncol(mat_87096))

# step6_Intersect common genes & merge 
common_genes <- Reduce(intersect, list(
  rownames(mat_144259), rownames(mat_50760), rownames(mat_87096)
))
message("Common genes across all 3 datasets: ", length(common_genes))

mat_144259   <- mat_144259[common_genes, ]
mat_50760    <- mat_50760[common_genes, ]
mat_87096    <- mat_87096[common_genes, ]

raw_counts   <- cbind(mat_144259, mat_50760, mat_87096)
raw_counts   <- raw_counts[rowSums(raw_counts) > 0, ]
message("Merged raw matrix: ", nrow(raw_counts), " genes x ", ncol(raw_counts), " samples")

# step7_Align metadata & 
common_samples <- intersect(colnames(raw_counts), rownames(final_metadata))
raw_counts     <- raw_counts[, common_samples]
meta           <- final_metadata[common_samples, , drop = FALSE]

stopifnot(all(colnames(raw_counts) == rownames(meta)))
stopifnot(!any(is.na(meta$condition)))

message(" Aligned: ", ncol(raw_counts), " samples | ",
        sum(meta$condition == "Tumor"), " Tumor | ",
        sum(meta$condition == "Normal"), " Normal")

# Save aligned metadata and raw count matrix 
write.csv(meta,
          "results_CRC/02_tables/metadata/CRC_metadata_final_used.csv",
          row.names = FALSE)

fwrite(as.data.frame(raw_counts) %>% rownames_to_column("GeneID"),
       "results_CRC/02_tables/counts/CRC_merged_raw_counts.tsv",
       sep = "\t")

# step8_ComBat-seq Batch Correction 
message("Running ComBat-seq batch correction...")
adjusted_counts <- ComBat_seq(
  counts = raw_counts,
  batch  = meta$batch,
  group  = meta$condition
)

fwrite(as.data.frame(adjusted_counts) %>% rownames_to_column("GeneID"),
       "results_CRC/02_tables/counts/CRC_ComBat_adjusted_counts.tsv",
       sep = "\t")

message(" ComBat-seq done. Adjusted counts saved.")

# See which genes/samples are affected
which_exceed <- which(adjusted_counts > .Machine$integer.max, arr.ind = TRUE)
cat("Affected genes:\n")
print(rownames(adjusted_counts)[unique(which_exceed[,"row"])])
cat("Affected samples:\n")
print(colnames(adjusted_counts)[unique(which_exceed[,"col"])])

# Cap those 37 values at integer max
adjusted_counts[adjusted_counts > .Machine$integer.max] <- .Machine$integer.max

# Convert storage mode to integer (required by DESeq2)
storage.mode(adjusted_counts) <- "integer"

# Verify clean
cat("\nAfter fix:\n")
cat("  NA values:                   ", sum(is.na(adjusted_counts)), "\n")
cat("  Values > integer max:        ", sum(adjusted_counts > .Machine$integer.max, na.rm=TRUE), "\n")
cat("  Storage mode:                ", storage.mode(adjusted_counts), "\n")
# Expected output:
# NA values:            0
# Values > integer max: 0
# Storage mode:         integer


# step9_DESeq2 object & pre-filter 
dds <- DESeqDataSetFromMatrix(
  countData = adjusted_counts,
  colData   = meta,
  design    = ~ condition
)

dds$condition <- relevel(dds$condition, ref = "Normal")

keep <- rowSums(counts(dds) >= 10) >= 5
dds  <- dds[keep, ]
message("Genes after pre-filtering: ", nrow(dds))

# step10_Run DESeq2 & VST transformation 
dds      <- DESeq(dds, parallel = FALSE)
vst_mat  <- vst(dds, blind = FALSE)
vst_data <- assay(vst_mat)

# step11_Extract results & map gene symbols 
res <- results(dds,
               contrast = c("condition", "Tumor", "Normal"),
               alpha    = 0.05)

message("DESeq2 result summary:")
print(summary(res))

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(
    symbol = unname(mapIds(org.Hs.eg.db,
                           keys      = gene,
                           keytype   = "ENTREZID",
                           column    = "SYMBOL",
                           multiVals = "first")),
    regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>% arrange(padj)

# step12_Non-coding filter 
noncoding_patterns <- c(
  "^LOC[0-9]","^MIR[0-9]","^LINC[0-9]","^MT-","^SNORD[0-9]","^SNORA[0-9]",
  "^RNU[0-9]","^RN7SL","^SNHG[0-9]","^SCARNA[0-9]","^MALAT","^NEAT",
  "^H19$","^XIST$","^MEG[0-9]","^PEG[0-9]"
)

is_noncoding <- function(symbols) {
  pattern <- paste(noncoding_patterns, collapse = "|")
  result  <- rep(FALSE, length(symbols))
  valid   <- !is.na(symbols)
  result[valid] <- grepl(pattern, symbols[valid], ignore.case = FALSE)
  result
}

res_df$is_noncoding <- is_noncoding(res_df$symbol)

sig_deg         <- res_df %>% filter(regulation != "NS") %>% arrange(padj)
sig_deg_protein <- sig_deg %>% filter(!is_noncoding, !is.na(symbol))

message("Significant DEGs (all):            ", nrow(sig_deg))
message("Significant DEGs (protein-coding): ", nrow(sig_deg_protein))
message("  Upregulated:   ", sum(sig_deg_protein$regulation == "Up"))
message("  Downregulated: ", sum(sig_deg_protein$regulation == "Down"))

write.csv(res_df,          "results_CRC/02_tables/DEG/DESeq2_all_genes_with_symbols.csv",         row.names = FALSE)
write.csv(sig_deg,         "results_CRC/02_tables/DEG/DESeq2_significant_DEGs_all.csv",           row.names = FALSE)
write.csv(sig_deg_protein, "results_CRC/02_tables/DEG/DESeq2_significant_DEGs_proteinCoding.csv", row.names = FALSE)

# step13_VST with gene symbols 
entrez_to_symbol <- mapIds(org.Hs.eg.db,
                           keys      = rownames(vst_data),
                           keytype   = "ENTREZID",
                           column    = "SYMBOL",
                           multiVals = "first")

vst_symbol <- as.data.frame(vst_data) %>%
  rownames_to_column("entrez_id") %>%
  mutate(gene_id = unname(entrez_to_symbol[entrez_id])) %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, everything(), -entrez_id)

write.csv(vst_symbol,
          "results_CRC/02_tables/DEG/VST_normalized_expression_geneSymbol.csv",
          row.names = FALSE)

# step14_Save R objects 
saveRDS(dds,             "results_CRC/01_objects/CRC_dds.rds")
saveRDS(res_df,          "results_CRC/01_objects/CRC_res_df.rds")
saveRDS(sig_deg_protein, "results_CRC/01_objects/CRC_sig_deg_protein.rds")
saveRDS(vst_mat,         "results_CRC/01_objects/CRC_vst_mat.rds")
saveRDS(vst_data,        "results_CRC/01_objects/CRC_vst_data.rds")

# step15_Plot helper 
save_plot <- function(plot_obj, base_path, w, h) {
  ggsave(paste0(base_path, ".png"),  plot_obj, width = w, height = h, dpi = 600)
  ggsave(paste0(base_path, ".tiff"), plot_obj, width = w, height = h, dpi = 600,
         compression = "none")
}
pal <- c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey75")

# step16_MA Plot 
ma_plot <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, colour = regulation)) +
  geom_point(alpha = 0.4, size = 1.2, stroke = 0) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0,         linetype = "solid",  colour = "black", linewidth = 0.3) +
  scale_colour_manual(name = "Regulation", values = pal) +
  labs(x = "Mean expression (log\u2081\u2080)",
       y = "Log\u2082 fold change (Tumor vs Normal)") +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
    legend.position = "right", legend.justification = "top",
    legend.background = element_rect(colour = "black", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )
save_plot(ma_plot, "results_CRC/03_plots/DEG/MA_plot", 8, 6)

# step17_Volcano Plot 
top10_up <- res_df %>% filter(regulation == "Up",   !is.na(symbol), !is_noncoding) %>%
  arrange(padj) %>% slice_head(n = 10)
top10_dn <- res_df %>% filter(regulation == "Down", !is.na(symbol), !is_noncoding) %>%
  arrange(padj) %>% slice_head(n = 10)
top20_labeled <- bind_rows(top10_up, top10_dn)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = regulation)) +
  geom_point(alpha = 0.5, size = 1.5, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = top20_labeled, aes(label = symbol),
    size = 3.5, fontface = "bold.italic",
    max.overlaps = 20, box.padding = 0.5, point.padding = 0.3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    name   = "Regulation",
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70"),
    labels = c("Up" = "Up-regulated", "Down" = "Down-regulated", "NS" = "NS")
  ) +
  labs(
    x = expression(bold(Log[2]~Fold~Change~(Tumor~vs~Normal))),
    y = expression(bold(-log[10](italic(p)[adj])))
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
    legend.position = "right", legend.justification = "top",
    legend.background = element_rect(colour = "black", linewidth = 0.5, fill = "white"),
    panel.grid.minor = element_blank()
  )
save_plot(volcano_plot, "results_CRC/03_plots/DEG/Volcano_plot", 8, 7)

# step18_PCA Plots 
pca_df     <- plotPCA(vst_mat, intgroup = c("condition", "batch"), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

pca_cond <- ggplot(pca_df, aes(PC1, PC2, colour = condition, fill = condition)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.12, colour = NA) +
  stat_ellipse(level = 0.95, linewidth = 1.1) +
  scale_colour_manual(name = "Condition", values = c("Normal" = "#2166AC", "Tumor" = "#D6604D")) +
  scale_fill_manual(  name = "Condition", values = c("Normal" = "#2166AC", "Tumor" = "#D6604D")) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(colour = "black", linewidth = 0.5, fill = "white"),
    panel.grid.minor = element_blank()
  )
save_plot(pca_cond, "results_CRC/03_plots/QC_PCA/PCA_by_Condition", 10, 7)

pca_batch <- ggplot(pca_df, aes(PC1, PC2, colour = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       colour = "GEO Study") +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"),
    legend.position = "right",
    legend.background = element_rect(colour = "black", linewidth = 0.5, fill = "white"),
    panel.grid.minor = element_blank()
  )
save_plot(pca_batch, "results_CRC/03_plots/QC_PCA/PCA_by_Batch", 10, 7)

# step19_ComplexHeatmap Top 25 Up + Top 25 Down
top25_up <- res_df %>%
  filter(!is.na(padj), !is.na(symbol), !is_noncoding, regulation == "Up")   %>%
  arrange(padj) %>% slice_head(n = 25)
top25_dn <- res_df %>%
  filter(!is.na(padj), !is.na(symbol), !is_noncoding, regulation == "Down") %>%
  arrange(padj) %>% slice_head(n = 25)
top_degs <- bind_rows(top25_dn, top25_up)

in_vst   <- top_degs$gene %in% rownames(vst_data)
top_filt <- top_degs[in_vst, ]

col_order <- c(rownames(meta)[meta$condition == "Normal"],
               rownames(meta)[meta$condition == "Tumor"])

hmat          <- vst_data[top_filt$gene, col_order, drop = FALSE]
rownames(hmat) <- top_filt$symbol
hmat_z        <- t(scale(t(hmat)))
valid         <- apply(hmat_z, 1, function(x) !all(is.na(x)) & !any(is.infinite(x)))
hmat_z        <- hmat_z[valid, , drop = FALSE]

reg_vec <- top_filt$regulation[match(rownames(hmat_z), top_filt$symbol)]
row_ord <- c(which(reg_vec == "Down"), which(reg_vec == "Up"))
hmat_z  <- hmat_z[row_ord, ]
reg_vec <- reg_vec[row_ord]

col_split        <- factor(meta[col_order, "condition"], levels = c("Normal", "Tumor"))
condition_colors <- c("Normal" = "#4575B4", "Tumor" = "#FFD700")
n_batch          <- length(levels(meta$batch))
batch_colors     <- setNames(brewer.pal(max(3, n_batch), "Set2")[1:n_batch], levels(meta$batch))
col_fun          <- colorRamp2(c(-2, 0, 2), c("green3", "black", "red"))

ann_df <- data.frame(
  Condition = meta[col_order, "condition"],
  Dataset   = meta[col_order, "batch"],
  row.names = col_order
)

ha_top <- HeatmapAnnotation(
  Condition = ann_df$Condition, Dataset = ann_df$Dataset,
  col = list(Condition = condition_colors, Dataset = batch_colors),
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 10),
  show_annotation_name = TRUE, show_legend = FALSE
)

ha_left <- rowAnnotation(
  Regulation = reg_vec,
  col = list(Regulation = c("Up" = "#8BC34A", "Down" = "#CE93D8")),
  width = unit(0.4, "cm"),
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 10),
  annotation_name_side = "top",
  show_annotation_name = TRUE, show_legend = FALSE
)

ht <- Heatmap(
  hmat_z, name = "Z-score", col = col_fun,
  top_annotation = ha_top, left_annotation = ha_left,
  cluster_rows = FALSE, cluster_columns = TRUE,
  cluster_column_slices = FALSE, column_split = col_split,
  column_gap = unit(0, "mm"), rect_gp = gpar(col = NA),
  show_row_names = TRUE, show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp   = gpar(fontface = "bold.italic", fontsize = 9),
  show_heatmap_legend = FALSE,
  row_split = factor(reg_vec, levels = c("Down","Up")),
  row_gap = unit(0, "mm"), border = TRUE
)

lgd_condition  <- Legend(title = "Condition",  at = c("Normal","Tumor"),
                         legend_gp = gpar(fill = c("#4575B4","#FFD700")),
                         title_gp = gpar(fontface="bold",fontsize=10), labels_gp = gpar(fontface="bold",fontsize=9))
lgd_regulation <- Legend(title = "Regulation", at = c("Down","Up"),
                         labels = c("Down-regulated","Up-regulated"),
                         legend_gp = gpar(fill = c("#CE93D8","#8BC34A")),
                         title_gp = gpar(fontface="bold",fontsize=10), labels_gp = gpar(fontface="bold",fontsize=9))
lgd_dataset    <- Legend(title = "Dataset", at = levels(meta$batch),
                         legend_gp = gpar(fill = batch_colors),
                         title_gp = gpar(fontface="bold",fontsize=10), labels_gp = gpar(fontface="bold",fontsize=9))
lgd_zscore     <- Legend(title = "Z-score", col_fun = col_fun, at = c(-2,0,2),
                         labels = c("-2","0","2"), legend_height = unit(3.5,"cm"), direction = "vertical",
                         title_gp = gpar(fontface="bold",fontsize=10), labels_gp = gpar(fontface="bold",fontsize=9))
all_legends    <- packLegend(lgd_condition, lgd_regulation, lgd_dataset, lgd_zscore,
                             direction = "vertical", gap = unit(4,"mm"))

for (ext in c("png","tiff")) {
  if (ext == "png") png( "results_CRC/03_plots/Heatmap/Heatmap_top25up_top25dn.png",
                         width=14, height=11, units="in", res=600)
  else               tiff("results_CRC/03_plots/Heatmap/Heatmap_top25up_top25dn.tiff",
                          width=14, height=11, units="in", res=600, compression="none")
  draw(ht, annotation_legend_list = all_legends,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       align_annotation_legend = "heatmap_top", merge_legend = FALSE,
       padding = unit(c(2,2,2,2),"mm"))
  dev.off()
}

# step20_GO & KEGG Enrichment 
entrez_all <- unique(as.character(sig_deg_protein$gene[!is.na(sig_deg_protein$gene)]))
cat("Entrez IDs for enrichment:", length(entrez_all), "\n")

table_dir <- "results_CRC/02_tables/enrichment"
plot_dir  <- "results_CRC/03_plots/enrichment"

save_3fmt <- function(plot, fname_base, width, height) {
  ggsave(file.path(plot_dir,"png", paste0(fname_base,".png")),
         plot, width=width, height=height, dpi=600, units="in")
  ggsave(file.path(plot_dir,"svg", paste0(fname_base,".svg")),
         plot, width=width, height=height)
  ggsave(file.path(plot_dir,"tiff",paste0(fname_base,".tiff")),
         plot, width=width, height=height, dpi=600, units="in", compression="none")
  cat(sprintf("  [SAVED] %s (.png / .svg / .tiff)\n", fname_base))
}

make_dotplot <- function(enrich_obj, title_str, top_n = 20) {
  df <- as.data.frame(enrich_obj) %>%
    arrange(p.adjust) %>% slice_head(n = top_n) %>%
    arrange(Count) %>%
    mutate(Description = str_wrap(Description, width = 45),
           Description = factor(Description, levels = unique(Description)))
  if (nrow(df) == 0) { message("No terms: ", title_str); return(NULL) }
  p_min <- min(df$p.adjust, na.rm=TRUE); p_max <- max(df$p.adjust, na.rm=TRUE)
  if (p_min == p_max) p_max <- p_min * 1.01
  ggplot(df, aes(x = Count, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.95) +
    scale_color_gradientn(
      name = "p.adjust",
      colors = c("#2D004B","#6A0573","#C2185B","#F06292","#FFCCBC"),
      values = scales::rescale(c(p_min, p_min+(p_max-p_min)*0.25,
                                 p_min+(p_max-p_min)*0.50,
                                 p_min+(p_max-p_min)*0.75, p_max)),
      limits = c(p_min, p_max),
      guide  = guide_colorbar(barheight=unit(5,"cm"), barwidth=unit(0.55,"cm"),
                              title.position="top", title.hjust=0.5,
                              frame.colour="black", ticks.colour="black")
    ) +
    scale_size_continuous(name = "GeneCount", range = c(3,14), breaks = c(5,10,15,20,25)) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.25))) +
    labs(x = "Count", y = NULL, title = title_str) +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face="bold", hjust=0.5, size=16, color="black"),
      axis.text.x      = element_text(face="bold", size=13, color="black"),
      axis.text.y      = element_text(face="bold", size=12, color="black"),
      axis.title.x     = element_text(face="bold", size=14, color="black"),
      legend.title     = element_text(face="bold", size=11, color="black"),
      legend.text      = element_text(face="bold", size=10, color="black"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="grey90", linewidth=0.3),
      panel.border     = element_rect(color="black", linewidth=0.8, fill=NA)
    )
}

analyses <- list(
  list(type="GO",   ont="BP", label="GO_BP", title="GO Biological Process"),
  list(type="GO",   ont="CC", label="GO_CC", title="GO Cellular Component"),
  list(type="GO",   ont="MF", label="GO_MF", title="GO Molecular Function"),
  list(type="KEGG", ont=NULL, label="KEGG",  title="KEGG Pathway Enrichment")
)

for (a in analyses) {
  cat(sprintf("\n Running %s...\n", a$label))
  if (a$type == "GO") {
    enrich_res <- tryCatch(
      enrichGO(gene=entrez_all, OrgDb=org.Hs.eg.db, ont=a$ont,
               pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE),
      error = function(e) NULL)
  } else {
    enrich_res <- tryCatch({
      res_k <- enrichKEGG(gene=entrez_all, organism="hsa",
                          pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2)
      if (!is.null(res_k) && nrow(as.data.frame(res_k)) > 0)
        setReadable(res_k, OrgDb=org.Hs.eg.db, keyType="ENTREZID") else NULL
    }, error = function(e) NULL)
  }
  if (is.null(enrich_res) || nrow(as.data.frame(enrich_res)) == 0) {
    cat(sprintf("   No significant terms for %s\n", a$label)); next
  }
  cat(sprintf("   %s: %d terms\n", a$label, nrow(as.data.frame(enrich_res))))
  write.csv(as.data.frame(enrich_res),
            file.path(table_dir, sprintf("%s_CRC_MetaAnalysis.csv", a$label)),
            row.names = FALSE)
  p <- make_dotplot(enrich_res, title_str = a$title, top_n = 20)
  if (!is.null(p))
    save_3fmt(p, fname_base = sprintf("%s_CRC", a$label),
              width = 11, height = ifelse(a$label == "KEGG", 9, 10))
}


setwd("C:/Users/akhil/comet_downlaod/colorectal_bulk_rna/data")

# --- GSE144259: 3 control + 3 diseased (exclude metastasis M) ----------------
meta_144259 <- data.frame(
  sample_id    = c("GSM4284531","GSM4284532",
                   "GSM4284534","GSM4284535",
                   "GSM4284537","GSM4284538"),
  geo_accession = c("GSM4284531","GSM4284532",
                    "GSM4284534","GSM4284535",
                    "GSM4284537","GSM4284538"),
  title        = c("CRC1N","CRC1T","CRC2N","CRC2T","CRC3N","CRC3T"),
  patient_id   = c(1, 1, 2, 2, 3, 3),
  batch        = "GSE144259",
  platform     = "GPL11154",
  condition    = c("control","diseased","control","diseased",
                   "control","diseased"),
  stringsAsFactors = FALSE
)

# --- GSE50760: 18 control + 18 diseased (exclude metastasis suffix -3) -------
patient_ids_50760 <- c(2,3,5,6,7,8,9,10,12,13,17,18,19,20,21,22,23,24)

meta_50760 <- data.frame(
  sample_id     = c(paste0("GSM", 1228184:1228201),   # tumor   -1
                    paste0("GSM", 1228202:1228219)),   # normal  -2
  geo_accession = c(paste0("GSM", 1228184:1228201),
                    paste0("GSM", 1228202:1228219)),
  title         = c(paste0("AMC_", patient_ids_50760, "-T"),
                    paste0("AMC_", patient_ids_50760, "-N")),
  patient_id    = c(patient_ids_50760, patient_ids_50760),
  batch         = "GSE50760",
  platform      = "GPL11154",
  condition     = c(rep("diseased", 18), rep("control", 18)),
  stringsAsFactors = FALSE
)

# --- GSE87096: 6 control + 6 diseased (RNA-seq only — exclude MeDIP/hMeDIP) --
meta_87096 <- data.frame(
  sample_id     = c(
    "GSM2322216","GSM2322222","GSM2322228",   # normal rna
    "GSM2322234","GSM2322240","GSM2322246",
    "GSM2322219","GSM2322225","GSM2322231",   # tumor  rna
    "GSM2322237","GSM2322243","GSM2322249"
  ),
  geo_accession = c(
    "GSM2322216","GSM2322222","GSM2322228",
    "GSM2322234","GSM2322240","GSM2322246",
    "GSM2322219","GSM2322225","GSM2322231",
    "GSM2322237","GSM2322243","GSM2322249"
  ),
  title         = c(
    "2504588_N","2512618_N","2539382_N","2551349_N","2553763_N","S12_N",
    "2504588_T","2512618_T","2539382_T","2551349_T","2553763_T","S12_T"
  ),
  patient_id    = c(1, 2, 3, 4, 5, 6,
                    1, 2, 3, 4, 5, 6),
  batch         = "GSE87096",
  platform      = "GPL11154",
  condition     = c(rep("control",  6), rep("diseased", 6)),
  stringsAsFactors = FALSE
)

# --- Combine all 3 datasets 
meta_all3_combined <- bind_rows(meta_144259, meta_50760, meta_87096)

# Verify pattern
message("Total samples: ", nrow(meta_all3_combined))
print(table(meta_all3_combined$batch, meta_all3_combined$condition))

# --- Save 
write.csv(meta_all3_combined,
          "results_CRC/02_tables/metadata/meta_all3_combined.csv",
          row.names = FALSE)

message(" Saved: results_CRC/02_tables/metadata/meta_all3_combined.csv")

# STRING ENSP → Gene Symbols | BioMart primary + manual fallback
library(biomaRt)
library(dplyr)

setwd("C:/Users/akhil/comet_downlaod/colorectal_bulk_rna/data")

string_ids <- c(
  "9606.ENSP00000287598",
  "9606.ENSP00000481380",
  "9606.ENSP00000302530",
  "9606.ENSP00000260731",
  "9606.ENSP00000358813",
  "9606.ENSP00000260363",
  "9606.ENSP00000247191",
  "9606.ENSP00000256442"
)
ensp_ids <- gsub("^9606\\.", "", string_ids)


mart <- NULL
for (mirror in c("useast", "uswest", "asia", "www")) {
  mart <- tryCatch(
    useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror),
    error = function(e) NULL
  )
  if (!is.null(mart)) { cat("Connected via mirror:", mirror, "\n"); break }
}

if (is.null(mart)) stop(" All BioMart mirrors failed — check internet connection")

# --- Query BioMart 
bm_result <- getBM(
  attributes = c("ensembl_peptide_id",
                 "ensembl_gene_id",
                 "external_gene_name",
                 "entrezgene_id",
                 "description"),
  filters    = "ensembl_peptide_id",
  values     = ensp_ids,
  mart       = mart
)

# --- Merge with original order 
result <- data.frame(string_id = string_ids,
                     ensp_id   = ensp_ids,
                     stringsAsFactors = FALSE) %>%
  left_join(bm_result, by = c("ensp_id" = "ensembl_peptide_id")) %>%
  rename(gene_symbol      = external_gene_name,
         ensg_id          = ensembl_gene_id,
         entrez_id        = entrezgene_id,
         gene_description = description)


manual_map <- list(
  "ENSP00000481380" = list(
    gene_symbol      = "CCNA2",
    ensg_id          = "ENSG00000145386",
    entrez_id        = 890L,
    gene_description = "cyclin A2 [Source:HGNC Symbol;Acc:HGNC:1578]"
  )
)

for (ensp in names(manual_map)) {
  idx <- which(result$ensp_id == ensp & is.na(result$gene_symbol))
  if (length(idx) > 0) {
    result$gene_symbol[idx]      <- manual_map[[ensp]]$gene_symbol
    result$ensg_id[idx]          <- manual_map[[ensp]]$ensg_id
    result$entrez_id[idx]        <- manual_map[[ensp]]$entrez_id
    result$gene_description[idx] <- manual_map[[ensp]]$gene_description
    cat("Manual fix applied for:", ensp, "→", manual_map[[ensp]]$gene_symbol, "\n")
  }
}


# CRC WGCNA ANALYSIS 


library(WGCNA)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggExtra)
library(patchwork)   
library(ggplotify)   
library(svglite)     

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Working directory & folder structure
setwd("C:/Users/akhil/comet_downlaod/colorectal_bulk_rna/data")

dir.create("results_CRC/01_objects",            recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/02_tables/WGCNA",       recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/WGCNA/png",    recursive = TRUE, showWarnings = FALSE)
dir.create("results_CRC/03_plots/WGCNA/svg",    recursive = TRUE, showWarnings = FALSE)

# Load Metadata and Expression Matrix
meta <- read.csv("results_CRC/02_tables/metadata/meta_all3_combined.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

vst_symbol <- read.csv(
  "results_CRC/02_tables/DEG/VST_normalized_expression_geneSymbol.csv",
  row.names   = 1,
  check.names = FALSE
)

datExpr_raw <- as.matrix(vst_symbol)
datExpr     <- t(datExpr_raw)

# Align metadata and counts perfectly
common_samples <- intersect(rownames(datExpr), rownames(meta))
datExpr        <- datExpr[common_samples, ]
meta_wgcna     <- meta[common_samples, ]

message("Matched samples : ", nrow(datExpr))
message("Genes for WGCNA : ", ncol(datExpr))

# Filter low-variance genes (top 5000)
gene_vars  <- apply(datExpr, 2, var)
top_genes  <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
datExpr    <- datExpr[, top_genes]
message("Genes after variance filter: ", ncol(datExpr))

# Goodness of sample check
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message("Removed bad samples/genes. Remaining: ", nrow(datExpr), " x ", ncol(datExpr))
}

# Scale-Free Topology Analysis
powers <- c(1:10, seq(12, 50, by = 2))
sft    <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
sft_df <- sft$fitIndices

# MANUAL SOFT POWER SELECTION
# Power set to 14 (WGCNA recommendation for signed networks of this size when SFT curve is flat)
soft_power <- 14
message("Manually selected Soft Power: ", soft_power)

# Plot 1A: SFT Fit
p_sft <- ggplot(sft_df, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 2.5, color = "black") +
  geom_text(vjust = -0.8, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.7) +
  labs(x = "Soft Threshold (Power)", y = "Scale-Free Topology Model Fit (signed R²)", title = "Scale independence") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold"))

# Plot 1B: Mean Connectivity
p_mean_conn <- ggplot(sft_df, aes(Power, mean.k., label = Power)) +
  geom_point(size = 2.5, color = "black") +
  geom_text(vjust = -0.8, size = 3.5, fontface = "bold") +
  labs(x = "Soft Threshold (Power)", y = "Mean Connectivity", title = "Mean connectivity") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold"))

# Build Network & Identify Modules
temp_cor <- cor
cor      <- WGCNA::cor

bwnet <- blockwiseModules(
  datExpr, maxBlockSize = 10000, TOMType = "signed", power = soft_power,
  mergeCutHeight = 0.25, numericLabels = FALSE, randomSeed = 1234, verbose = 3
)
cor <- temp_cor 

moduleColors <- bwnet$colors
MEs          <- bwnet$MEs

# Save Module Assignments
write.csv(data.frame(gene = colnames(datExpr), module = moduleColors),
          "results_CRC/02_tables/WGCNA/Module_Gene_Assignments.csv", row.names = FALSE)

# Plot 1C: Dendrogram 
p_dendro <- as.ggplot(expression({
  par(font.main = 2, font.lab = 2, font.axis = 2, cex.main = 1.2, cex.lab = 1.1)
  plotDendroAndColors(
    bwnet$dendrograms[[1]], moduleColors[bwnet$blockGenes[[1]]], "Module Colors",
    main = "Cluster Dendrogram", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
  )
}))

# Plot 1D: Module Correlation Heatmap 
module_cor  <- cor(MEs, use = "p")
module_dist <- 1 - module_cor
rownames(module_dist) <- gsub("ME", "", rownames(module_dist))
colnames(module_dist) <- gsub("ME", "", colnames(module_dist))
col_fun_mod <- colorRamp2(c(0, 0.75, 1.5), c("red", "white", "blue"))

ht_mod <- Heatmap(
  module_dist, name = "Distance", col = col_fun_mod,
  cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_distance_rows = as.dist(module_dist), clustering_distance_columns = as.dist(module_dist),
  row_names_gp = gpar(fontsize = 9, fontface = "bold"), column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param = list(title = "Distance", title_gp = gpar(fontface = "bold")), border = TRUE
)
p_mod_heat <- as.ggplot(expression(draw(ht_mod)))

# MERGED FIGURE 1 (A, B, C, D) 
layout_fig1 <- "
AABB
CCCC
DDDD
"
fig1 <- p_sft + p_mean_conn + p_dendro + p_mod_heat + 
  plot_layout(design = layout_fig1) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("results_CRC/03_plots/WGCNA/png/Fig1_Network_Topology.png", fig1, width = 14, height = 18, dpi = 600)
ggsave("results_CRC/03_plots/WGCNA/svg/Fig1_Network_Topology.svg", fig1, width = 14, height = 18)
message("Merged Figure 1 (Network Topology) saved.")

# Module-Trait Relationships
traits <- data.frame(
  Diseased = as.numeric(meta_wgcna$condition == "diseased"),
  Control  = as.numeric(meta_wgcna$condition == "control"),
  row.names = rownames(meta_wgcna)
)

moduleTraitCor    <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

rownames(moduleTraitCor)    <- gsub("ME", "", rownames(moduleTraitCor))
rownames(moduleTraitPvalue) <- gsub("ME", "", rownames(moduleTraitPvalue))

ord               <- order(moduleTraitCor[, "Diseased"], decreasing = TRUE)
moduleTraitCor    <- moduleTraitCor[ord, ]
moduleTraitPvalue <- moduleTraitPvalue[ord, ]

# Formatting cell text for Heatmap
cell_text <- matrix("", nrow(moduleTraitCor), ncol(moduleTraitCor))
for (i in seq_len(nrow(moduleTraitCor))) {
  for (j in seq_len(ncol(moduleTraitCor))) {
    r_val <- round(moduleTraitCor[i, j], 2)
    p_fmt <- formatC(moduleTraitPvalue[i, j], format = "e", digits = 0)
    p_fmt <- sub("e\\+0*(\\d+)", "e+\\1", sub("e-0*(\\d+)", "e-\\1", p_fmt))
    cell_text[i, j] <- paste0(r_val, "\n(", p_fmt, ")")
  }
}

# Plot 2A: Trait Heatmap
mod_colors <- rownames(moduleTraitCor)
row_anno   <- rowAnnotation(Module = mod_colors, col = list(Module = setNames(mod_colors, mod_colors)), show_legend = FALSE)
col_fun_trait <- colorRamp2(c(-1, 0, 1), c("#33A02C", "white", "#E31A1C"))

ht_trait <- Heatmap(
  moduleTraitCor, name = "Correlation", col = col_fun_trait,
  cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", left_annotation = row_anno,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"), column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(cell_text[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))
  },
  heatmap_legend_param = list(title = "Correlation", direction = "vertical"), border = TRUE
)
p_trait_heat <- as.ggplot(expression(draw(ht_trait)))

# GS vs MM Scatter Plots
geneTraitSignificance <- as.data.frame(cor(datExpr, traits$Diseased, use = "p"))
names(geneTraitSignificance) <- "GS.Trait"
rownames(geneTraitSignificance) <- colnames(datExpr)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
rownames(geneModuleMembership) <- colnames(datExpr)

# Extract top 6 modules highly correlated with 'Diseased'
top_modules <- rownames(moduleTraitCor)[order(abs(moduleTraitCor[, "Diseased"]), decreasing = TRUE)]
top_modules <- top_modules[1:min(6, length(top_modules))] 

scatter_list <- list()
for (mod in top_modules) {
  inModule   <- moduleColors == mod
  col_idx    <- match(paste0("ME", mod), colnames(MEs))
  
  plot_data  <- data.frame(MM = geneModuleMembership[inModule, col_idx], 
                           GS = geneTraitSignificance$GS.Trait[inModule]) %>% 
    filter(complete.cases(.))
  
  if(nrow(plot_data) < 3) next
  
  r_val <- cor(plot_data$MM, plot_data$GS)
  p_val <- cor.test(plot_data$MM, plot_data$GS)$p.value
  
  p_base <- ggplot(plot_data, aes(x = MM, y = GS)) +
    geom_point(color = "navy", alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +
    annotate("text", x = min(plot_data$MM), y = max(plot_data$GS),
             label = paste0("r = ", round(r_val, 2), "\np = ", format(p_val, scientific = TRUE, digits = 2)),
             hjust = 0, vjust = 1, size = 3.5, fontface = "bold") +
    labs(x = paste0("Module Membership (", mod, ")"), y = "Gene Significance (Diseased)") +
    theme_classic() + 
    theme(axis.title = element_text(face = "bold"), panel.border = element_rect(color = "black", fill = NA))
  
  p_marg <- ggMarginal(p_base, type = "histogram", fill = "orange", xparams = list(fill = "orange"), yparams = list(fill = "red"))
  
  scatter_list[[mod]] <- as.ggplot(p_marg)
}

# MERGED FIGURE 2 (A, B-G) 
p_scatters_grid <- wrap_plots(scatter_list, ncol = 2)

fig2 <- p_trait_heat + p_scatters_grid + 
  plot_layout(widths = c(1, 2)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave("results_CRC/03_plots/WGCNA/png/Fig2_Module_Trait_Scatters.png", fig2, width = 16, height = 12, dpi = 600)
ggsave("results_CRC/03_plots/WGCNA/svg/Fig2_Module_Trait_Scatters.svg", fig2, width = 16, height = 12)
message("Merged Figure 2 (Trait Heatmap & Scatters) saved.")

# Extract Hub Modules
top_pos_mod <- top_modules[moduleTraitCor[top_modules, "Diseased"] > 0][1]
top_neg_mod <- top_modules[moduleTraitCor[top_modules, "Diseased"] < 0][1]

if(!is.na(top_pos_mod)) {
  write.table(colnames(datExpr)[moduleColors == top_pos_mod], 
              paste0("results_CRC/02_tables/WGCNA/Genes_", top_pos_mod, "_Upregulated.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("Extracted top upregulated module genes: ", top_pos_mod)
}

if(!is.na(top_neg_mod)) {
  write.table(colnames(datExpr)[moduleColors == top_neg_mod], 
              paste0("results_CRC/02_tables/WGCNA/Genes_", top_neg_mod, "_Downregulated.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("Extracted top downregulated module genes: ", top_neg_mod)
}

# Save R Objects
saveRDS(list(bwnet = bwnet, MEs = MEs, moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue),
        "results_CRC/01_objects/WGCNA_objects.rds")

message("CRC WGCNA PIPELINE COMPLETE.")



