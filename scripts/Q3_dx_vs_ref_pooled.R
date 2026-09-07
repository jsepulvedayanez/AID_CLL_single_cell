suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
  library(DESeq2)
  library(tibble)
})

set.seed(42)
merged <- readRDS("results/CLL_scRNAseq_merged_annotated.rds")

# Fix proliferation: use MKI67 (gold standard), not Phase (Phase was 1/3-split -> unreliable)
mki67 <- FetchData(merged, vars = "MKI67")
merged$is_proliferating <- mki67$MKI67 > 0 & merged$tumor_final

cat("MKI67+ tumor cells:", sum(merged$is_proliferating),
    "=", round(100*sum(merged$is_proliferating)/sum(merged$tumor_final),2), "% of tumor\n")
cat("MKI67+ by status:\n"); print(table(merged$status[merged$is_proliferating]))
cat("MKI67+ by cluster:\n"); print(table(merged$cell_type[merged$is_proliferating]))

tumor <- subset(merged, subset = tumor_final == TRUE)
cat("\nTumor cells: dx", sum(tumor$status=="diagnosis"),
    "| ref", sum(tumor$status=="refractory"), "\n")
cat("Refractory tumor B by subject:\n")
print(table(tumor$subject[tumor$status=="refractory"]))

# ---------- Q3a: module scores Dx vs Ref (pooled) ----------
scores <- c("IFNa1","IFNg1","KRAS1","AID_program1","Prolif1")
res_list <- lapply(scores, function(s) {
  dx <- tumor@meta.data[[s]][tumor$status=="diagnosis"]
  ref <- tumor@meta.data[[s]][tumor$status=="refractory"]
  wt <- suppressWarnings(wilcox.test(dx, ref))
  tibble(score = s,
         median_dx = median(dx), median_ref = median(ref),
         mean_dx = mean(dx), mean_ref = mean(ref),
         p_value = wt$p.value)
})
score_res <- bind_rows(res_list)
cat("\n========== Q3: PATHWAY SCORES, tumor Dx vs Ref (pooled) ==========\n")
print(as.data.frame(score_res))
write.csv(score_res, "results/Q3_pathway_scores_dx_vs_ref.csv", row.names = FALSE)

# ---------- Q3b: fgsea Dx vs Ref (pooled), hallmark + KEGG ----------
hallmark_sets <- msigdbr(species="Homo sapiens", collection="H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
kegg_sets <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:KEGG_LEGACY") %>%
  split(x = .$gene_symbol, f = .$gs_name)

mat <- GetAssayData(tumor, layer = "data")
cells_dx <- colnames(tumor)[tumor$status=="diagnosis"]
cells_ref <- colnames(tumor)[tumor$status=="refractory"]
lfc <- rowMeans(mat[, cells_dx, drop=FALSE]) - rowMeans(mat[, cells_ref, drop=FALSE])
lfc <- sort(lfc[is.finite(lfc)], decreasing = TRUE)

gsea_h <- fgsea(pathways=hallmark_sets, stats=lfc, minSize=15, maxSize=500, eps=0)
gsea_h <- gsea_h[order(gsea_h$pval), ]
gsea_k <- fgsea(pathways=kegg_sets, stats=lfc, minSize=15, maxSize=500, eps=0)
gsea_k <- gsea_k[order(gsea_k$pval), ]

cat("\n========== Q3: GSEA Hallmark, tumor Dx vs Ref (pooled, NES>0 = up in Dx) ==========\n")
top_h <- gsea_h %>% filter(padj < 0.25) %>% select(pathway, NES, pval, padj) %>% as.data.frame()
print(head(top_h, 30))
write.csv(as.data.frame(gsea_h)[, setdiff(names(gsea_h), "leadingEdge")], "results/GSEA_Hallmark_Dx_vs_Ref_pooled.csv", row.names=FALSE)
write.csv(as.data.frame(gsea_k)[, setdiff(names(gsea_k), "leadingEdge")], "results/GSEA_KEGG_Dx_vs_Ref_pooled.csv", row.names=FALSE)

# Focus on the 3 named pathways
cat("\n--- Named pathways of interest ---\n")
focus <- gsea_h %>% filter(grepl("INTERFERON_ALPHA|INTERFERON_GAMMA|KRAS", pathway)) %>%
  select(pathway, NES, pval, padj) %>% as.data.frame()
print(focus)
focus_k <- gsea_k %>% filter(grepl("KRAS|JAK_STAT|INTERFERON", pathway)) %>%
  select(pathway, NES, pval, padj) %>% as.data.frame()
cat("\nKEGG related:\n"); print(head(focus_k, 15))

# ---------- Q3c: pseudobulk DESeq2 (pooled, design ~ subject + status) ----------
cat("\n========== Q3: pseudobulk DESeq2 (design ~ subject + status) ==========\n")
run_pb <- function(obj, cell_subset, comp_col, a, b, sample_col="sample_id") {
  sub <- subset(obj, subset = cell_type %in% cell_subset)
  pb <- AggregateExpression(sub, assays="RNA", slot="counts",
                            group.by=sample_col, return.seurat=FALSE)$RNA
  colnames(pb) <- gsub("-", "_", colnames(pb))   # Seurat replaces '_' in idents
  # Build colData in base R (dplyr pipeline trips on rownames)
  md <- unique(sub@meta.data[, c(sample_col, comp_col, "subject")])
  md <- md[md[[comp_col]] %in% c(a, b), , drop = FALSE]
  rownames(md) <- md[[sample_col]]
  common <- intersect(colnames(pb), rownames(md))
  pb <- pb[, common, drop = FALSE]
  md <- md[common, , drop = FALSE]
  md[[comp_col]] <- factor(md[[comp_col]], levels = c(b, a))
  dds <- DESeqDataSetFromMatrix(countData = round(pb), colData = md,
                                design = as.formula(paste0("~ subject + ", comp_col)))
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c(comp_col, a, b)) %>% as.data.frame() %>%
    rownames_to_column("gene") %>% arrange(padj)
  list(dds = dds, results = res, n_samples = ncol(pb))
}

pb_res <- tryCatch(
  run_pb(merged, c("CLL_B","CLL_B_plasma","CLL_B_prolif"),
         "status","diagnosis","refractory"),
  error = function(e) { cat("Pseudobulk DESeq2 error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(pb_res)) {
  cat("Pseudobulk samples:", pb_res$n_samples, "\n")
  write.csv(pb_res$results, "results/DESeq2_pseudobulk_Dx_vs_Ref_pooled.csv", row.names=FALSE)
  sig <- pb_res$results %>% filter(!is.na(padj), padj < 0.05)
  cat("DE genes (padj<0.05):", nrow(sig), "\n")
  print(head(sig %>% select(gene, log2FoldChange, padj), 30))
}

cat("\n========== DONE ==========\n")
