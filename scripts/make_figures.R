suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

merged <- readRDS("results/CLL_scRNAseq_merged_annotated.rds")
mki67 <- FetchData(merged, vars = "MKI67")
merged$is_proliferating <- mki67$MKI67 > 0 & merged$is_tumorB

# 1. Corrected annotation UMAP
p1 <- DimPlot(merged, reduction = "umap", group.by = "cell_type",
              label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle("Corrected annotation (12 clusters)") + theme_classic()

# 2. Status + subject
p2 <- DimPlot(merged, reduction = "umap", group.by = "status", pt.size = 0.4) +
  ggtitle("Diagnosis vs Refractory") + theme_classic()

# 3. AICDA expression (rare AID+ subpopulation)
p3 <- FeaturePlot(merged, features = "AICDA", reduction = "umap",
                  order = TRUE, min.cutoff = 0, max.cutoff = 1) &
  scale_color_viridis(option = "magma")
p3 <- p3 + ggtitle("AICDA (AID) expression")

# 4. MKI67 (proliferation)
p4 <- FeaturePlot(merged, features = "MKI67", reduction = "umap",
                  order = TRUE, min.cutoff = 0, max.cutoff = 1) &
  scale_color_viridis(option = "magma")
p4 <- p4 + ggtitle("MKI67 (proliferation)")

ggsave("figures/01_annotation_umap.png", p1, width = 7, height = 6, dpi = 150)
ggsave("figures/02_status_umap.png", p2, width = 7, height = 6, dpi = 150)
ggsave("figures/03_AICDA_feature.png", p3, width = 7, height = 6, dpi = 150)
ggsave("figures/04_MKI67_feature.png", p4, width = 7, height = 6, dpi = 150)

# 5. Pathway scores Dx vs Ref (tumor B cells)
tumor <- subset(merged, subset = is_tumorB == TRUE)
p5 <- VlnPlot(tumor, features = c("IFNa1","IFNg1","KRAS1"),
              group.by = "status", ncol = 3, pt.size = 0) +
  plot_annotation(title = "Pathway module scores: tumor, Dx vs Ref (pooled)")
ggsave("figures/05_pathway_scores_dx_vs_ref.png", p5, width = 11, height = 5, dpi = 150)

cat("Figures saved to figures/\n")
