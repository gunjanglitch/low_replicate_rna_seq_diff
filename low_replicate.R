# ==========================================================
# 1. LOAD LIBRARIES
# ==========================================================

library(airway)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(gprofiler2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggrepel)


# ==========================================================
# 2. LOAD AIRWAY DATASET
# ==========================================================

data(airway)

counts <- assay(airway)
metadata <- colData(airway)

dim(counts)
dim(metadata)

metadata[, c("SampleName", "cell", "dex")]

table(metadata$cell, metadata$dex)


# ==========================================================
#                 BEFORE: FULL DATASET
#              8 SAMPLES / 4 CELL LINES
# ==========================================================


# ==========================================================
# 3. CREATE DESEQ2 OBJECT
# ==========================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ cell + dex
)


# ==========================================================
# 4. FILTER LOW-COUNT GENES
# ==========================================================

keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep, ]

dim(dds)


# ==========================================================
# 5. NORMALIZATION
# ==========================================================

dds <- estimateSizeFactors(dds)

sizeFactors(dds)


# ==========================================================
# 6. VST TRANSFORMATION
# ==========================================================

vsd <- vst(dds)


# ==========================================================
# 7. PCA
# ==========================================================

pca <- prcomp(t(assay(vsd)))

summary(pca)

plotPCA(
  vsd,
  intgroup = c("dex", "cell")
)


# ==========================================================
# 8. DIFFERENTIAL EXPRESSION ANALYSIS
# ==========================================================

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("dex", "trt", "untrt")
)

summary(res)


# ==========================================================
# 9. DEFINE DEGs
#    padj < 0.05
#    |log2FC| > 1
# ==========================================================

deg <- res[
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1,
]

nrow(deg)

sum(deg$log2FoldChange > 1)
sum(deg$log2FoldChange < -1)


# ==========================================================
# 10. GENE ID → GENE SYMBOL
# ==========================================================

res$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

deg$gene_symbol <- res$gene_symbol[
  match(rownames(deg), rownames(res))
]


# ==========================================================
# 11. DEG TABLE
# ==========================================================

deg_table <- as.data.frame(deg)

deg_table$ensembl_id <- rownames(deg_table)

rownames(deg_table) <- NULL

deg_table <- deg_table[, c(
  "ensembl_id",
  "gene_symbol",
  "baseMean",
  "log2FoldChange",
  "pvalue",
  "padj"
)]

head(deg_table)


# ==========================================================
# 12. TOP 10 GENES
# ==========================================================

top10 <- deg_table[
  order(deg_table$padj),
][1:10, ]

top10[, c(
  "gene_symbol",
  "log2FoldChange",
  "padj"
)]


# ==========================================================
# 13. VOLCANO PLOT — BEFORE
# ==========================================================

top10_ids <- top10$ensembl_id

top10_symbols <- res$gene_symbol[
  match(top10_ids, rownames(res))
]

top10_symbols <- top10_symbols[
  !is.na(top10_symbols)
]

EnhancedVolcano(
  res,
  lab = res$gene_symbol,
  selectLab = top10_symbols,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "BEFORE: Dexamethasone vs Untreated",
  subtitle = "8 samples / 4 cell lines",
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  drawConnectors = TRUE,
  boxedLabels = TRUE,
  max.overlaps = 20
)


# ==========================================================
# 14. MA PLOT — BEFORE
# ==========================================================

ma_data_before <- as.data.frame(res)

ma_data_before$ensembl_id <- rownames(ma_data_before)

ma_data_before$gene_symbol <- res$gene_symbol

ma_data_before$group <- "Not significant"

ma_data_before$group[
  ma_data_before$padj < 0.05 &
    ma_data_before$log2FoldChange > 1
] <- "Upregulated"

ma_data_before$group[
  ma_data_before$padj < 0.05 &
    ma_data_before$log2FoldChange < -1
] <- "Downregulated"


ggplot(
  ma_data_before,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange,
    color = group
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_hline(
    yintercept = c(-1, 0, 1),
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_text_repel(
    data = ma_data_before[
      ma_data_before$ensembl_id %in% top10$ensembl_id &
        !is.na(ma_data_before$gene_symbol),
    ],
    aes(label = gene_symbol),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "blue",
      "Not significant" = "grey",
      "Upregulated" = "red"
    )
  ) +
  labs(
    title = "BEFORE: Dexamethasone vs Untreated",
    subtitle = "8 samples / 4 cell lines",
    x = "Log10 Mean Expression",
    y = "Log2 Fold Change",
    color = "Gene status"
  ) +
  theme_classic()


# ==========================================================
# 15. HEATMAP — BEFORE
# ==========================================================

top50 <- deg_table[
  order(deg_table$padj),
][1:50, ]

heatmap_data <- assay(vsd)[
  top50$ensembl_id,
]

rownames(heatmap_data) <- top50$gene_symbol

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = as.data.frame(
    colData(vsd)[, c("dex", "cell")]
  ),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  main = "Top 50 DEGs — BEFORE"
)


# ==========================================================
# 16. SEPARATE UPREGULATED / DOWNREGULATED GENES
# ==========================================================

up_genes <- rownames(
  deg[deg$log2FoldChange > 1, ]
)

down_genes <- rownames(
  deg[deg$log2FoldChange < -1, ]
)

length(up_genes)
length(down_genes)


# ==========================================================
# 17. GO + KEGG — UPREGULATED
# ==========================================================

go_kegg_up <- gost(
  query = up_genes,
  organism = "hsapiens",
  sources = c(
    "GO:BP",
    "GO:MF",
    "GO:CC",
    "KEGG"
  )
)

table(go_kegg_up$result$source)


# ==========================================================
# 18. GO + KEGG — DOWNREGULATED
# ==========================================================

go_kegg_down <- gost(
  query = down_genes,
  organism = "hsapiens",
  sources = c(
    "GO:BP",
    "GO:MF",
    "GO:CC",
    "KEGG"
  )
)

table(go_kegg_down$result$source)


# ==========================================================
# 19. GO PLOT — UPREGULATED
#    BP + CC + MF
# ==========================================================

go_up <- go_kegg_up$result[
  go_kegg_up$result$source %in%
    c("GO:BP", "GO:CC", "GO:MF"),
]

go_up_top <- do.call(
  rbind,
  lapply(
    c("GO:BP", "GO:CC", "GO:MF"),
    function(x) {

      tmp <- go_up[
        go_up$source == x,
      ]

      tmp[
        order(tmp$p_value),
      ][1:min(10, nrow(tmp)), ]

    }
  )
)


ggplot(
  go_up_top,
  aes(
    x = -log10(p_value),
    y = reorder(term_name, -log10(p_value)),
    color = source,
    size = intersection_size
  )
) +
  geom_point() +
  labs(
    title = "GO Enrichment — Upregulated Genes",
    x = "-log10(P-value)",
    y = "GO Term",
    color = "GO Category",
    size = "Gene Count"
  ) +
  theme_classic()


# ==========================================================
# 20. KEGG PLOT — UPREGULATED
# ==========================================================

kegg_up <- go_kegg_up$result[
  go_kegg_up$result$source == "KEGG",
]

kegg_up <- kegg_up[
  order(kegg_up$p_value),
]


ggplot(
  kegg_up,
  aes(
    x = intersection_size,
    y = reorder(term_name, intersection_size)
  )
) +
  geom_col(fill = "steelblue") +
  geom_text(
    aes(label = intersection_size),
    hjust = -0.2
  ) +
  labs(
    title = "KEGG Enrichment — Upregulated Genes",
    x = "Number of DEGs",
    y = "KEGG Pathway"
  ) +
  theme_classic()


# ==========================================================
# 21. GO PLOT — DOWNREGULATED
# ==========================================================

go_down <- go_kegg_down$result[
  go_kegg_down$result$source %in%
    c("GO:BP", "GO:CC", "GO:MF"),
]

go_down_top <- do.call(
  rbind,
  lapply(
    c("GO:BP", "GO:CC", "GO:MF"),
    function(x) {

      tmp <- go_down[
        go_down$source == x,
      ]

      tmp[
        order(tmp$p_value),
      ][1:min(10, nrow(tmp)), ]

    }
  )
)


ggplot(
  go_down_top,
  aes(
    x = -log10(p_value),
    y = reorder(term_name, -log10(p_value)),
    color = source,
    size = intersection_size
  )
) +
  geom_point() +
  labs(
    title = "GO Enrichment — Downregulated Genes",
    x = "-log10(P-value)",
    y = "GO Term",
    color = "GO Category",
    size = "Gene Count"
  ) +
  theme_classic()


# ==========================================================
# 22. KEGG PLOT — DOWNREGULATED
# ==========================================================

kegg_down <- go_kegg_down$result[
  go_kegg_down$result$source == "KEGG",
]

kegg_down <- kegg_down[
  order(kegg_down$p_value),
]


ggplot(
  kegg_down,
  aes(
    x = -log10(p_value),
    y = reorder(term_name, -log10(p_value)),
    size = intersection_size,
    color = -log10(p_value)
  )
) +
  geom_point() +
  labs(
    title = "KEGG Enrichment — Downregulated Genes",
    x = "-log10(P-value)",
    y = "KEGG Pathway",
    size = "Gene Count",
    color = "-log10(P-value)"
  ) +
  theme_classic()



############################################################
#                 AFTER: REDUCED DATASET
#              4 SAMPLES / 2 CELL LINES
############################################################


# ==========================================================
# 23. SELECT 2 CELL LINES
# ==========================================================

selected_samples <- c(
  "SRR1039508",
  "SRR1039509",
  "SRR1039512",
  "SRR1039513"
)

counts_after <- counts[
  ,
  selected_samples
]

metadata_after <- metadata[
  selected_samples,
]

metadata_after[
  ,
  c("SampleName", "cell", "dex")
]


# ==========================================================
# 24. CREATE DESEQ2 OBJECT — AFTER
# ==========================================================

dds_after <- DESeqDataSetFromMatrix(
  countData = counts_after,
  colData = metadata_after,
  design = ~ cell + dex
)


# ==========================================================
# 25. FILTER LOW-COUNT GENES — AFTER
# ==========================================================

keep_after <- rowSums(
  counts(dds_after) >= 10
) >= 2

dds_after <- dds_after[
  keep_after,
]

dim(dds_after)


# ==========================================================
# 26. NORMALIZATION — AFTER
# ==========================================================

dds_after <- estimateSizeFactors(
  dds_after
)

sizeFactors(dds_after)


# ==========================================================
# 27. VST — AFTER
# ==========================================================

vsd_after <- vst(
  dds_after
)


# ==========================================================
# 28. PCA — AFTER
# ==========================================================

plotPCA(
  vsd_after,
  intgroup = c("dex", "cell")
)


# ==========================================================
# 29. DIFFERENTIAL EXPRESSION — AFTER
# ==========================================================

dds_after <- DESeq(
  dds_after
)

res_after <- results(
  dds_after,
  contrast = c("dex", "trt", "untrt")
)


# ==========================================================
# 30. DEFINE AFTER DEGs
# ==========================================================

deg_after <- res_after[
  !is.na(res_after$padj) &
    res_after$padj < 0.05 &
    abs(res_after$log2FoldChange) > 1,
]

nrow(deg_after)

sum(deg_after$log2FoldChange > 1)

sum(deg_after$log2FoldChange < -1)


# ==========================================================
# 31. GENE ID → GENE SYMBOL — AFTER
# ==========================================================

res_after$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_after),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

deg_after$gene_symbol <- res_after$gene_symbol[
  match(
    rownames(deg_after),
    rownames(res_after)
  )
]


# ==========================================================
# 32. AFTER DEG TABLE
# ==========================================================

deg_after_table <- as.data.frame(
  deg_after
)

deg_after_table$ensembl_id <- rownames(
  deg_after_table
)

rownames(deg_after_table) <- NULL

deg_after_table <- deg_after_table[
  ,
  c(
    "ensembl_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )
]

head(deg_after_table)


# ==========================================================
# 33. TOP 10 AFTER GENES
# ==========================================================

top10_after <- deg_after_table[
  order(deg_after_table$padj),
][1:10, ]

top10_after[
  ,
  c(
    "gene_symbol",
    "log2FoldChange",
    "padj"
  )
]


# ==========================================================
# 34. VOLCANO PLOT — AFTER
# ==========================================================

top10_ids_after <- top10_after$ensembl_id

top10_symbols_after <- res_after$gene_symbol[
  match(
    top10_ids_after,
    rownames(res_after)
  )
]

top10_symbols_after <- top10_symbols_after[
  !is.na(top10_symbols_after)
]


EnhancedVolcano(
  res_after,
  lab = res_after$gene_symbol,
  selectLab = top10_symbols_after,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "AFTER: Dexamethasone vs Untreated",
  subtitle = "4 samples / 2 cell lines",
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  drawConnectors = TRUE,
  boxedLabels = TRUE,
  max.overlaps = 20
)


# ==========================================================
# 35. MA PLOT — AFTER
# ==========================================================

ma_data_after <- as.data.frame(
  res_after
)

ma_data_after$ensembl_id <- rownames(
  ma_data_after
)

ma_data_after$gene_symbol <- res_after$gene_symbol

ma_data_after$group <- "Not significant"

ma_data_after$group[
  ma_data_after$padj < 0.05 &
    ma_data_after$log2FoldChange > 1
] <- "Upregulated"

ma_data_after$group[
  ma_data_after$padj < 0.05 &
    ma_data_after$log2FoldChange < -1
] <- "Downregulated"


ggplot(
  ma_data_after,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange,
    color = group
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_hline(
    yintercept = c(-1, 0, 1),
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_text_repel(
    data = ma_data_after[
      ma_data_after$ensembl_id %in%
        top10_after$ensembl_id &
        !is.na(ma_data_after$gene_symbol),
    ],
    aes(label = gene_symbol),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "blue",
      "Not significant" = "grey",
      "Upregulated" = "red"
    )
  ) +
  labs(
    title = "AFTER: Dexamethasone vs Untreated",
    subtitle = "4 samples / 2 cell lines",
    x = "Log10 Mean Expression",
    y = "Log2 Fold Change",
    color = "Gene status"
  ) +
  theme_classic()


# ==========================================================
# 36. HEATMAP — AFTER
# ==========================================================

top50_after <- deg_after_table[
  order(deg_after_table$padj),
][1:50, ]

heatmap_data_after <- assay(
  vsd_after
)[
  top50_after$ensembl_id,
]

rownames(heatmap_data_after) <-
  top50_after$gene_symbol


pheatmap(
  heatmap_data_after,
  scale = "row",
  annotation_col = as.data.frame(
    colData(vsd_after)[
      ,
      c("dex", "cell")
    ]
  ),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  main = "Top 50 DEGs — AFTER"
)


###############################################
#BEFORE vs AFTER
###############################################
before_genes <- rownames(deg)
after_genes <- rownames(deg_after)
length(before_genes)
length(after_genes)

#find shared and unique DEG
common_genes <- intersect(
  before_genes,
  after_genes
)

before_only <- setdiff(
  before_genes,
  after_genes
)

after_only <- setdiff(
  after_genes,
  before_genes
)
length(common_genes)
length(before_only)
length(after_only)

#simple overlap plot
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    BEFORE = before_genes,
    AFTER = after_genes
  ),
  filename = NULL,
  fill = c("steelblue", "tomato"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "DEG Overlap: BEFORE vs AFTER"
)

grid::grid.draw(venn.plot)

#calculate what fraction of the AFTER DEGs are also present in BEFORE:
length(common_genes) / length(after_genes) * 100
#And what fraction of the BEFORE DEGs remain in AFTER:
length(common_genes) / length(before_genes) * 100

