# 05_deseq2.R
# Import Salmon quantifications (tximport), run DESeq2 with a paired design
# (donor cell line + dexamethasone treatment), make QC/result plots, and
# validate against the published airway count table.
# Run from the repo root (open the RStudio project), e.g. source("scripts/05_deseq2.R")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(pheatmap)
})
set.seed(1)
outdir <- "results/deseq2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Sample sheet -------------------------------------------------------
samples <- read.delim("samples.tsv")
samples$cell <- factor(samples$cell)
samples$dex  <- factor(samples$dex, levels = c("untrt", "trt"))  # untreated = reference
rownames(samples) <- samples$run

# ---- 2. Import Salmon transcript quantifications, summarise to genes -------
tx2gene <- read.delim("data/ref/tx2gene.tsv")
files <- file.path("results/salmon", samples$run, "quant.sf")
names(files) <- samples$run
stopifnot(all(file.exists(files)))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")])

# ---- 3. DESeq2: paired design (controls for donor), test dex effect --------
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ cell + dex)
keep <- rowSums(counts(dds) >= 10) >= 4   # >=10 reads in at least 4 samples (smallest group)
dds <- dds[keep, ]
cat("Genes kept after filtering:", nrow(dds), "\n")
dds <- DESeq(dds)

res     <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
res_shr <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")
summary(res)

sym <- unique(tx2gene[, c("GENEID", "SYMBOL")])
out <- data.frame(
  gene_id       = rownames(res),
  symbol        = sym$SYMBOL[match(rownames(res), sym$GENEID)],
  baseMean      = res$baseMean,
  log2FC        = res$log2FoldChange,
  log2FC_shrunk = res_shr$log2FoldChange,
  pvalue        = res$pvalue,
  padj          = res$padj
)
out <- out[order(out$padj), ]
write.csv(out, file.path(outdir, "dex_vs_untrt_results.csv"), row.names = FALSE)

# ---- 4. Plots ---------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = dex, shape = cell)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", pv[1], "% variance"),
       y = paste0("PC2: ", pv[2], "% variance"),
       title = "PCA of VST-transformed counts") +
  theme_bw()
ggsave(file.path(outdir, "pca.png"), p, width = 6, height = 4.5, dpi = 150)

png(file.path(outdir, "ma_plot_shrunk.png"), width = 1000, height = 800, res = 150)
plotMA(res_shr, ylim = c(-5, 5), main = "Dexamethasone vs untreated (apeglm-shrunken LFC)")
dev.off()

top <- head(out$gene_id[!is.na(out$padj)], 30)
mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- out$symbol[match(top, out$gene_id)]
ann <- as.data.frame(colData(vsd)[, c("dex", "cell")])
pheatmap(mat, annotation_col = ann, main = "Top 30 DE genes (row-centred VST)",
         filename = file.path(outdir, "top30_heatmap.png"), width = 7, height = 8)

# ---- 5. Sanity check: known glucocorticoid-responsive genes -----------------
markers <- c("FKBP5", "TSC22D3", "PER1", "KLF15", "DUSP1", "ZBTB16", "CRISPLD2")
cat("\nKnown dexamethasone-responsive genes:\n")
print(out[out$symbol %in% markers, c("symbol", "log2FC_shrunk", "padj")], row.names = FALSE)

# ---- 6. Validation against the published airway count table -----------------
# The Bioconductor 'airway' package holds the authors' gene counts (full depth,
# GRCh37 / Ensembl 75, genome alignment). Differences are expected: this analysis
# uses 2M read pairs per sample, GENCODE v50 protein-coding transcripts, and Salmon.
if (requireNamespace("airway", quietly = TRUE)) {
  suppressPackageStartupMessages(library(airway))
  data("airway")
  airway$dex <- factor(airway$dex, levels = c("untrt", "trt"))
  dds_pub <- DESeqDataSet(airway, design = ~ cell + dex)
  dds_pub <- dds_pub[rowSums(counts(dds_pub) >= 10) >= 4, ]
  dds_pub <- DESeq(dds_pub, quiet = TRUE)
  res_pub <- results(dds_pub, contrast = c("dex", "trt", "untrt"), alpha = 0.05)

  ours <- data.frame(id = sub("\\..*$", "", out$gene_id), lfc_ours = out$log2FC, padj_ours = out$padj)
  ours <- ours[!duplicated(ours$id), ]
  pub  <- data.frame(id = rownames(res_pub), lfc_pub = res_pub$log2FoldChange, padj_pub = res_pub$padj)
  m <- merge(ours, pub, by = "id")

  sig_ours <- m$id[which(m$padj_ours < 0.05)]
  sig_pub  <- m$id[which(m$padj_pub  < 0.05)]
  both     <- intersect(sig_ours, sig_pub)
  either   <- m$id %in% union(sig_ours, sig_pub)
  rho <- cor(m$lfc_ours[either], m$lfc_pub[either], method = "spearman", use = "complete.obs")

  val <- sprintf(paste0(
    "Genes present in both analyses: %d\n",
    "Significant (padj < 0.05): this pipeline %d | published counts %d | overlap %d\n",
    "Share of this pipeline's DE genes also DE in published counts: %.1f%%\n",
    "Spearman correlation of log2FC (genes DE in either): %.3f\n"),
    nrow(m), length(sig_ours), length(sig_pub), length(both),
    100 * length(both) / max(length(sig_ours), 1), rho)
  cat("\n", val, sep = "")
  writeLines(val, file.path(outdir, "validation_vs_published.txt"))

  pv2 <- ggplot(m[either, ], aes(lfc_pub, lfc_ours)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    labs(x = "log2FC, published airway counts", y = "log2FC, this pipeline",
         title = sprintf("Agreement with published counts (Spearman rho = %.2f)", rho)) +
    theme_bw()
  ggsave(file.path(outdir, "validation_scatter.png"), pv2, width = 5.5, height = 5, dpi = 150)
} else {
  message("Package 'airway' not installed - skipping validation. BiocManager::install('airway')")
}

# ---- 7. Record software versions --------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("\nDone. Outputs in", outdir, "\n")
