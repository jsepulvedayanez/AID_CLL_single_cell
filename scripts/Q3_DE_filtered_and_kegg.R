suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(msigdbr)
  library(fgsea)
  library(DESeq2)
})

merged <- readRDS("results/CLL_scRNAseq_merged_annotated.rds")
bmark <- FetchData(merged, vars = c("MS4A1","CD79A","CD19"))
merged$is_tumorB <- merged$celltype_demux == "Bcell" &
  (bmark$MS4A1 > 0 | bmark$CD79A > 0 | bmark$CD19 > 0)
tumor <- subset(merged, subset = is_tumorB == TRUE)

# ---------- 1. Filtrar genes no informativos para DE ----------
drop_re <- "^(IGH|IGK|IGL|IGJ|JCHAIN|TRA|TRB|TRD|TRG|MT-|MT|RPS|RPL|HBA|HBB|HBD|HBG|HBE|HBZ|MALAT1)"
keep_genes <- rownames(tumor)[!grepl(drop_re, rownames(tumor))]
cat("Genes totales:", nrow(tumor), "| tras filtrar IG/TCR/mito/ribo/hb:",
    length(keep_genes), "\n")
tumor_f <- tumor[keep_genes, ]

# ---------- 2. Pseudobulk DESeq2 (genes filtrados) ----------
pb <- AggregateExpression(tumor_f, assays = "RNA", slot = "counts",
                          group.by = "sample_id", return.seurat = FALSE)$RNA
colnames(pb) <- gsub("-", "_", colnames(pb))

md <- unique(tumor_f@meta.data[, c("sample_id","status","subject")])
rownames(md) <- md$sample_id
md <- md[colnames(pb), ]
md$status <- factor(md$status, levels = c("refractory","diagnosis"))

dds <- DESeqDataSetFromMatrix(countData = round(pb), colData = md,
                              design = ~ subject + status)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("status","diagnosis","refractory")) %>%
  as.data.frame() %>% rownames_to_column("gene") %>% arrange(padj)

sig <- res %>% filter(!is.na(padj), padj < 0.05)
cat("\nGenes DE (padj<0.05) tras filtrar:", nrow(sig), "\n")
cat("\nTop 30 genes DE:\n")
print(head(res %>% select(gene, log2FoldChange, padj), 30))

write.csv(res, "results/DESeq2_pseudobulk_Dx_vs_Ref_pooled_FILTERED.csv", row.names = FALSE)

# ---------- 3. GSEA KEGG + Hallmark (tumorB, sin genes filtrados) ----------
hallmark_sets <- msigdbr(species="Homo sapiens", collection="H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
kegg_sets <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:KEGG_LEGACY") %>%
  split(x = .$gene_symbol, f = .$gs_name)

mat <- GetAssayData(tumor, layer = "data")
lfc <- rowMeans(mat[, tumor$status == "diagnosis", drop=FALSE]) -
       rowMeans(mat[, tumor$status == "refractory", drop=FALSE])
lfc <- sort(lfc[is.finite(lfc)], decreasing = TRUE)

gsea_k <- fgsea(pathways=kegg_sets, stats=lfc, minSize=15, maxSize=500, eps=0)
gsea_k <- gsea_k[order(gsea_k$pval), ]

cat("\n========== TOP 20 KEGG (padj<0.25, NES>0 = up en Dx) ==========\n")
print(gsea_k %>% filter(padj < 0.25) %>%
        select(pathway, NES, pval, padj) %>% head(20))

write.csv(as.data.frame(gsea_k)[, setdiff(names(gsea_k),"leadingEdge")],
          "results/GSEA_KEGG_Dx_vs_Ref_pooled.csv", row.names=FALSE)

# Hallmark top 20 (para referencia completa)
gsea_h <- fgsea(pathways=hallmark_sets, stats=lfc, minSize=15, maxSize=500, eps=0)
gsea_h <- gsea_h[order(gsea_h$pval), ]
cat("\n========== TOP 20 HALLMARK (padj<0.25) ==========\n")
print(gsea_h %>% filter(padj < 0.25) %>%
        select(pathway, NES, pval, padj) %>% head(20))

cat("\nDONE\n")
