suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

merged <- readRDS("results/CLL_scRNAseq_merged_annotated.rds")
bmark <- FetchData(merged, vars = c("MS4A1","CD79A","CD19"))
merged$is_tumorB <- merged$celltype_demux == "Bcell" &
  (bmark$MS4A1 > 0 | bmark$CD79A > 0 | bmark$CD19 > 0)
tumor <- subset(merged, subset = is_tumorB == TRUE)

pal <- c("Diagnóstico" = "#3B82BF", "Refractario" = "#E64B35")

# Función de plot reutilizable
plot_pathways <- function(df, title, subtitle, n_dx, n_ref) {
  long <- df %>%
    pivot_longer(-status, names_to = "pathway", values_to = "score") %>%
    mutate(pathway = recode(pathway,
      IFNa1 = "IFN-α respuesta", IFNg1 = "IFN-γ respuesta",
      KRAS1 = "KRAS signaling"))
  long$pathway <- factor(long$pathway,
    levels = c("IFN-α respuesta", "IFN-γ respuesta", "KRAS signaling"))
  long$status <- factor(long$status, levels = c("Diagnóstico", "Refractario"))

  meds <- long %>% group_by(pathway, status) %>%
    summarise(med = median(score), .groups = "drop") %>%
    group_by(pathway) %>%
    mutate(rng = diff(range(long$score[long$pathway == pathway])),
           y_label = med + 0.18 * rng) %>% ungroup()

  pvals <- long %>% group_by(pathway) %>%
    summarise(p = wilcox.test(score[status=="Diagnóstico"],
                              score[status=="Refractario"])$p.value,
              .groups = "drop") %>%
    mutate(label = paste0("p = ", formatC(p, format = "e", digits = 2)))

  ggplot(long, aes(x = status, y = score, fill = status)) +
    geom_violin(trim = TRUE, alpha = 0.35, color = NA, scale = "width") +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9,
                 color = "grey20", fill = "white") +
    geom_point(data = meds, aes(x = status, y = med, fill = status),
               shape = 23, size = 3, color = "black", stroke = 0.7) +
    geom_label(data = meds, aes(x = status, y = y_label, label = sprintf("%.3f", med)),
               fill = "white", color = "black", size = 3.5, fontface = "bold",
               label.padding = unit(0.15, "lines"), inherit.aes = FALSE) +
    geom_text(data = pvals, aes(x = 1.5, y = Inf, label = label),
              vjust = 1.2, size = 3.4, inherit.aes = FALSE, fontface = "italic") +
    facet_wrap(~ pathway, scales = "free_y") +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = "Module score (z)", fill = NULL,
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 12),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40", size = 10))
}

# ---- 1. COMPARACIÓN COMPLETA ----
df_full <- FetchData(tumor, vars = c("status","IFNa1","IFNg1","KRAS1"))
df_full$status <- ifelse(df_full$status == "diagnosis", "Diagnóstico", "Refractario")
g1 <- plot_pathways(df_full,
  "Vías en B tumorales: Diagnóstico vs Refractario — COMPARACIÓN COMPLETA",
  paste0("Todas las células. n = ", sum(df_full$status=="Diagnóstico"),
         " Dx / ", sum(df_full$status=="Refractario"), " Ref"))
ggsave("figures/05_pathway_scores_dx_vs_ref.png", g1, width = 11, height = 4.5, dpi = 160)

# ---- 2. SUBSET JUSTO (submuestreo estratificado por sujeto, seed 1) ----
ref_cells <- colnames(tumor)[tumor$status == "refractory"]
dx_cells  <- colnames(tumor)[tumor$status == "diagnosis"]
ref_n <- table(tumor$subject[ref_cells])
set.seed(1)
sel <- unlist(lapply(names(ref_n), function(sb) {
  pool <- dx_cells[tumor$subject[dx_cells] == sb]
  sample(pool, size = min(as.integer(ref_n[sb]), length(pool)))
}))

df_sub <- FetchData(tumor, vars = c("status","IFNa1","IFNg1","KRAS1"))
keep <- c(sel, ref_cells)
df_sub <- df_sub[keep, ]
df_sub$status <- ifelse(df_sub$status == "diagnosis", "Diagnóstico", "Refractario")

g2 <- plot_pathways(df_sub,
  "Vías en B tumorales: Diagnóstico vs Refractario — SUBSET JUSTO (submuestreo)",
  paste0("Diagnóstico submuestreado por sujeto hasta igualar refractario. n = ",
         sum(df_sub$status=="Diagnóstico"), " Dx / ",
         sum(df_sub$status=="Refractario"), " Ref (seed 1 de 20)"))
ggsave("figures/06_pathway_scores_subset_433vs433.png", g2, width = 11, height = 4.5, dpi = 160)

cat("Figuras guardadas: 05 (completo) y 06 (subset).\n")
