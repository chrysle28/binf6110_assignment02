### ===============================
## BINF110 Assignment 2
## Differential expression and ORA analysis of yeast transcriptome dataset

## Kenneth Gamueda
## 2026-03-01

### === PACKAGES USED ===========
library("ggplot2")
library("tidyverse")
library("readr")
library("DESeq2")
library("tximport")
library("org.Sc.sgd.db")
library("rtracklayer")
library("txdbmaker")
library("apeglm")
library("ggpubr")
library("EnhancedVolcano")
library("pheatmap")
library("clusterProfiler")
library("enrichplot")
library("cowplot")
library("styler")

### === DIFFERENTIAL EXPRESSION ANALYSIS ========
## Importing files and setting up objects for downstream analysis
# Creating metadata table
df_yeast_metadata <- data.frame(
  sampleID = c("IL20", "IL21", "IL22", "IL23", "IL24", "IL25", "IL29", "IL30", "IL31"),
  stage = c("early", "early", "early", "thin", "thin", "thin", "mature", "mature", "mature"),
  SRA = c("SRR10551665", "SRR10551664", "SRR10551663", "SRR10551662", "SRR10551661", "SRR10551660", "SRR10551659", "SRR10551658", "SRR10551657")
)

# Importing quant files from Salmon
files <- file.path("quants", df_yeast_metadata$SRA, "quant.sf")
names(files) <- df_yeast_metadata$SRA
files

# Creating tx2gene table from GTF file
gtf_yeast <- import("genomic.gtf")
txdb <- makeTxDbFromGRanges(gtf_yeast)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

# Creating vcount table from Salmon quant files and tx2gene table
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# Creating DESeqDataSet object with stage as factor
dds_TXI <- DESeqDataSetFromTximport(txi,
  colData = df_yeast_metadata,
  design = ~stage
)
dds_TXI$stage <- factor(dds_TXI$stage, levels = c("early", "thin", "mature"))

## Visualizing overall data structure
vsd_overall <- vst(dds_TXI, blind = F)

# PCA plot of quantified data
PCA_data <- plotPCA(vsd_overall, intgroup = "stage", returnData = T)
percentVar <- round(100 * attr(PCA_data, "percentVar"))

ggplot(PCA_data, aes(PC1, PC2, color = stage)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of all samples (overall data structure)") +
  theme(plot.title = element_text(face = "bold", size = 18)) +
  coord_fixed()

# Sample-sample distance heatmap of data
sample_dists <- dist(t(assay(vsd_overall)))
sample_mat <- as.matrix(sample_dists)
rownames(sample_mat) <- colnames(dds_TXI)
pheatmap(sample_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample distance heatmap (overall structure)",
)

## Conducting differential expression analysis (pairwise comparisons for each stage)
dds_e <- DESeq(dds_TXI)
resultsNames(dds_e)

evt_res <- results(dds_e, alpha = 0.05, name = "stage_thin_vs_early")
evt_res # early (ref) vs thin

evm_res <- results(dds_e, alpha = 0.05, contrast = c("stage", "mature", "early"))
evm_res # early (ref) vs mature

# Conducting LFC shrinkage
evt_LFC <- lfcShrink(dds_e,
  coef = "stage_thin_vs_early",
  type = "apeglm"
)

evm_LFC <- lfcShrink(dds_e,
  coef = "stage_mature_vs_early",
  type = "apeglm"
)

# Repeating previous steps with thin as reference level
dds_TXI$stage <- relevel(dds_TXI$stage, ref = "thin")
dds_t <- DESeq(dds_TXI)
resultsNames(dds_t)

tvm_res <- results(dds_t, alpha = 0.05, contrast = c("stage", "mature", "thin"))
tvm_res # thin (ref) vs mature

tvm_LFC <- lfcShrink(dds_t,
  coef = "stage_mature_vs_thin",
  type = "apeglm"
)

## Summary of DESeq2 results
summary(evt_res)
sum(evt_res$padj < 0.05, na.rm = TRUE)

summary(evm_res)
sum(evm_res$padj < 0.05, na.rm = TRUE)

summary(tvm_res)
sum(tvm_res$padj < 0.05, na.rm = TRUE)

### === VISUALIZATION OF DE RESULTS ========
## Creating MA plots for each pairwise comparison
# Early vs thin (top selected by P-adj)
p1 <- ggmaplot(evt_LFC,
  fdr = 0.05,
  fc = 2,
  main = "Thin (test) vs Early (ref) - top genes selected via P-adj",
  top = 5,
  size = 2,
  alpha = 0.7,
  select.top.method = "padj",
  legend = "top",
  font.label = c("bold"),
  label.rectangle = T,
  font.legend = c(10, "bold"),
  font.main = c(14, "bold"),
  ggtheme = ggplot2::theme_minimal()
)

# Early vs thin (top selected by FC)
p2 <- ggmaplot(evt_LFC,
  fdr = 0.05,
  fc = 2,
  main = "Thin (test) vs Early (ref) - top genes selected via FC",
  top = 5,
  size = 2,
  alpha = 0.7,
  select.top.method = "fc",
  legend = "top",
  font.label = c("bold"),
  label.rectangle = T,
  font.legend = c(10, "bold"),
  font.main = c(14, "bold"),
  ggtheme = ggplot2::theme_minimal()
)

pA <- plot_grid(p1, p2)
pA

# Thin vs mature (top selected by P-adj)
p3 <- ggmaplot(tvm_LFC,
  fdr = 0.05,
  fc = 2,
  main = "Mature (test) vs Thin (ref) - top genes selected via P-adj",
  top = 5,
  size = 2,
  alpha = 0.7,
  select.top.method = "padj",
  legend = "top",
  font.label = c("bold"),
  label.rectangle = T,
  font.legend = c(10, "bold"),
  font.main = c(14, "bold"),
  ggtheme = ggplot2::theme_minimal()
)

# Thin vs mature (top selected by FC)
p4 <- ggmaplot(tvm_LFC,
  fdr = 0.05,
  fc = 2,
  main = "Mature (test) vs Thin (ref) - top genes selected via FC",
  top = 5,
  size = 2,
  alpha = 0.7,
  select.top.method = "fc",
  legend = "top",
  font.label = c("bold"),
  label.rectangle = T,
  font.legend = c(10, "bold"),
  font.main = c(14, "bold"),
  ggtheme = ggplot2::theme_minimal()
)

pB <- plot_grid(p3, p4)
pB

# Early vs mature (top selected by P-adj)
p5 <- ggmaplot(evm_LFC,
  fdr = 0.05,
  fc = 2,
  main = "Mature (test) vs Early (ref) - (top genes selected via P-adj)",
  top = 5,
  size = 2,
  alpha = 0.7,
  select.top.method = "padj",
  legend = "top",
  font.label = c("bold"),
  label.rectangle = T,
  font.legend = c(10, "bold"),
  font.main = c(14, "bold"),
  ggtheme = ggplot2::theme_minimal()
)

# Early vs mature (top selected by P-adj)
p6 <- ggmaplot(evm_LFC,
  fdr = 0.05,
  fc = 2,
  main = "Mature (test) vs Early (ref) -  top genes selected via FC",
  top = 5,
  size = 2,
  alpha = 0.7,
  select.top.method = "fc",
  legend = "top",
  font.label = c("bold"),
  label.rectangle = T,
  font.legend = c(10, "bold"),
  font.main = c(14, "bold"),
  ggtheme = ggplot2::theme_minimal()
)

pC <- plot_grid(p5, p6)
pC

plot_grid(pA, pB, pC, nrow = 3)
rm(p1, p2, p3, p4, p5, p6, pA, pB, pC)

## VOLCANO PLOTS
# Creating custom key-value pairs for plotting in Enhanced Volcano
create_keyvals <- function(LFC, padj_cutoff = 0.05, lfc_cutoff = 1) {
  keyvals <- ifelse(LFC$padj < padj_cutoff & abs(LFC$log2FoldChange) > lfc_cutoff,
    ifelse(LFC$log2FoldChange > 0, "darkred", "royalblue"),
    "gray"
  )
  keyvals[is.na(keyvals)] <- "gray"
  names(keyvals)[keyvals == "darkred"] <- "Up"
  names(keyvals)[keyvals == "royalblue"] <- "Down"
  names(keyvals)[keyvals == "gray"] <- "NS"
  return(keyvals)
}

keyvals_evt <- create_keyvals(evt_LFC)
keyvals_tvm <- create_keyvals(tvm_LFC)
keyvals_evm <- create_keyvals(evm_LFC)

# Early vs thin
p1 <- EnhancedVolcano(evt_LFC,
  lab = rownames(evt_LFC),
  title = "Early vs thin",
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 2,
  pointSize = 2,
  labSize = 4,
  cutoffLineType = "twodash",
  cutoffLineWidth = 0.8,
  colCustom = keyvals_evt,
  legendPosition = "top",
  boxedLabels = T,
  drawConnectors = T,
  gridlines.major = F,
  gridlines.minor = T,
  selectLab = c(
    "YGR087C", "YGR088W", "YHR094C", "YKR097W", "YCR105W",
    "YGR296W", "YEL077C", "YHR055C", "YHR052C-B", "YEL069C"
  )
)
p1

# Thin vs mature
p2 <- EnhancedVolcano(tvm_LFC,
  lab = rownames(tvm_LFC),
  title = "Thin vs mature",
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 2,
  pointSize = 2,
  labSize = 4,
  cutoffLineType = "twodash",
  cutoffLineWidth = 0.8,
  colCustom = keyvals_tvm,
  legendPosition = "top",
  boxedLabels = T,
  drawConnectors = T,
  gridlines.major = F,
  gridlines.minor = F,
  selectLab = c("YKR075C", "YPR127W", "YDR085C", "YEL070W", "YBR117C", "YEL072W", "YCR106W", "YCR107W", "YEL071W", "YCR105")
)
p2

# Early vs mature
p3 <- EnhancedVolcano(evm_LFC,
  lab = rownames(evm_LFC),
  title = "Early vs mature",
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 2,
  pointSize = 2,
  labSize = 4,
  cutoffLineType = "twodash",
  cutoffLineWidth = 0.8,
  colCustom = keyvals_evm,
  legendPosition = "top",
  boxedLabels = T,
  drawConnectors = T,
  gridlines.major = F,
  gridlines.minor = F,
  selectLab = c("YGR087C", "YHR094C", "YJL052W", "YGL055W", "YIR019C", "YGR296W", "YFL068W", "YDR545W", "YEL071W", "YEl070W")
)
p3

plot_grid(p1, p2, p3, ncol = 3)
rm(p1, p2, p3)

## HEATMAP
# Creating dataframes of LFC objects from each comparison
evt_df <- as.data.frame(evt_LFC)
tvm_df <- as.data.frame(tvm_LFC)
evm_df <- as.data.frame(evm_LFC)

# Combining dataframes to get list of all genes
all_genes <- union(union(rownames(evt_df), rownames(tvm_df)), rownames(evm_df))

# Creating combined dataframe of genes with their padj values
all_df <- data.frame(
  evt = evt_df[all_genes, "padj"],
  tvm = tvm_df[all_genes, "padj"],
  evm = evm_df[all_genes, "padj"]
)
rownames(all_df) <- all_genes

# Ordering genes by their padj value (smallest first)
all_df$min_padj <- apply(all_df, 1, function(x) min(x, na.rm = TRUE))
all_df <- all_df[order(all_df$min_padj), ]

# Extracting the 50 most significant genes
top_genes <- head(rownames(all_df), 30)
top_genes

# Extracting counts with variance stabilizing transformation
vsd <- vst(dds_e)
mat <- assay(vsd)[top_genes, ]

# Creating heatmap
annotation_colors <- list(stage = c(early = "#548292", thin = "#ca5152", mature = "#648f67"))

pheatmap(mat,
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  annotation_col = as.data.frame(colData(dds_e)[, "stage", drop = F]),
  annotation_colors = annotation_colors,
  annotation_names_col = F,
  show_colnames = F,
  main = "Heatmap of the top 30 most significant genes (based on P-adj)"
)

### === OVER-REPRESENTATION ANALYSIS ========
# Function to run workflow
run_enrichment_workflow <- function(df,
                                    padj_cutoff = 0.05,
                                    lfc_cutoff = 1,
                                    ontology = "BP",
                                    showCategories = 15,
                                    enrichment = "GO",
                                    comparison_name = "Comparison") {
  # Adding ORF as column from rownames
  df$ORF <- rownames(df)

  # Obtaining significant genes (based on padj and lfc)
  sig <- df %>%
    filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()

  # Obtaining background list of all genes
  all <- df %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()

  # Obtaining upregulated and downregulated genes
  up <- df %>%
    filter(padj < padj_cutoff & log2FoldChange > lfc_cutoff) %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()

  down <- df %>%
    filter(padj < padj_cutoff & log2FoldChange < -lfc_cutoff) %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()

  # GO enrichment
  if (enrichment == "GO") {
    # For significant genes
    res <- enrichGO(
      gene = sig,
      universe = all,
      OrgDb = org.Sc.sgd.db,
      keyType = "ORF",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )

    # For up/downregulated genes
    compare_res <- compareCluster(
      geneClusters = list(Upregulated = up, Downregulated = down),
      fun = "enrichGO",
      OrgDb = org.Sc.sgd.db,
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      keyType = "ORF"
    )

    # KEGG enrichment
  } else if (enrichment == "KEGG") {
    # For significant genes
    res <- enrichKEGG(
      gene = sig,
      organism = "sce",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )

    # For up/downregulated genes
    compare_res <- compareCluster(
      geneClusters = list(Upregulated = up, Downregulated = down),
      fun = "enrichKEGG",
      organism = "sce",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }

  # Results and plots
  list(
    enrich_res = res,
    compare_res = compare_res,
    e_bar = barplot(res, showCategory = showCategories) + ggtitle(paste(enrichment, ontology, "Enrichment", paste0("(", comparison_name, ")"))) +
      theme(plot.title = element_text(face = "bold", size = 16)),
    e_dot = dotplot(res, showCategory = showCategories) + ggtitle(paste(enrichment, ontology, "Enrichment", paste0("(", comparison_name, ")"))) +
      theme(plot.title = element_text(face = "bold", size = 16)),
    c_dot = dotplot(compare_res, showCategory = 10) + ggtitle(paste(enrichment, ontology, "Comparison", paste0("(", comparison_name, ")"))) +
      theme(plot.title = element_text(face = "bold", size = 16))
  )
}

# Early vs thin GO
evt_GO_results <- run_enrichment_workflow(evt_df, enrichment = "GO", ont = "BP", comparison_name = "Thin vs Early")
evt_GO_results$e_bar
evt_GO_results$e_dot
evt_GO_results$c_dot
plot_grid(evt_GO_results$e_bar, evt_GO_results$c_dot)

# Early vs thin KEGG
evt_KEGG_results <- run_enrichment_workflow(evt_df, enrichment = "KEGG", comparison_name = "Thin vs Early")
evt_KEGG_results$e_bar
evt_KEGG_results$e_dot
evt_KEGG_results$c_dot
plot_grid(evt_KEGG_results$e_bar, evt_KEGG_results$c_dot)

# Thin vs mature GO
tvm_GO_results <- run_enrichment_workflow(tvm_df, enrichment = "GO", ont = "BP", comparison_name = "Mature vs Thin")
tvm_GO_results$e_bar
tvm_GO_results$e_dot
tvm_GO_results$c_dot
plot_grid(tvm_GO_results$e_bar, tvm_GO_results$c_dot)

# Thin vs mature KEGG
tvm_KEGG_results <- run_enrichment_workflow(tvm_df, enrichment = "KEGG", comparison_name = "Mature vs Thin")
tvm_KEGG_results$e_bar
tvm_KEGG_results$e_dot
tvm_KEGG_results$c_dot
plot_grid(tvm_KEGG_results$e_bar, tvm_KEGG_results$c_dot)

# Early vs mature GO
evm_GO_results <- run_enrichment_workflow(evm_df, enrichment = "GO", ont = "BP", comparison_name = "Mature vs Early")
evm_GO_results$e_bar
evm_GO_results$e_dot
evm_GO_results$c_dot
plot_grid(evm_GO_results$e_bar, evm_GO_results$c_dot)

# Early vs mature KEGG
evm_KEGG_results <- run_enrichment_workflow(evm_df, enrichment = "KEGG", comparison_name = "Mature vs Early")
evm_KEGG_results$e_bar
evm_KEGG_results$e_dot
evm_KEGG_results$c_dot
plot_grid(evm_KEGG_results$e_bar, evm_KEGG_results$c_dot)

# Comparison of all stages with up/downregulation
all_stages <- list("Early vs Thin" = evt_df, "Thin vs Mature" = tvm_df, "Early vs Mature" = evm_df)

# Function to split into upregulated and downregulated genes
get_up_down <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1) {
  df$ORF <- rownames(df)
  up <- df %>%
    filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(ORF)
  down <- df %>%
    filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(ORF)
  return(list(Up = na.omit(up), Down = na.omit(down)))
}

# Creating list of genes to pass to compareCluster
gene_clusters <- unlist(
  lapply(all_stages, get_up_down,
    padj_cutoff = 0.05,
    lfc_cutoff = 2
  ),
  recursive = F
)

names(gene_clusters) <- gsub("\\.Up", " (Up)", names(gene_clusters))
names(gene_clusters) <- gsub("\\.Down", " (Down)", names(gene_clusters))

# Conducting GO enrichment analysis with all stages
all_GO <- compareCluster(
  geneCluster = gene_clusters,
  fun = "enrichGO",
  OrgDb = org.Sc.sgd.db,
  ont = "BP",
  pAdjustMethod = "BH",
  keyType = "ORF"
)

# Conducting KEGG analysis with all stages
all_KEGG <- compareCluster(
  geneClusters = gene_clusters,
  fun = "enrichKEGG",
  organism = "sce",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# Visualizing compareCluster results
dotplot(all_GO, showCategory = 5) +
  theme(plot.title = element_text(face = "bold", size = 18)) +
  ggtitle("Comparison of GO Enrichment Across Stages")

dotplot(all_KEGG) +
  theme(plot.title = element_text(face = "bold", size = 18)) +
  ggtitle("Comparison of KEGG Pathways Across Stages")
