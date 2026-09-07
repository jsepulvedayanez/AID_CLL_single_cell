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

df <- FetchData(tumor, vars = c("status", "IFNa1", "IFNg1", "KRAS1"))
df$status <- ifelse(df$status == "diagnosis", "Diagnóstico", "Refractario")
long <- df %>%
  pivot_longer(-status, names_to = "pathway", values_to = "score") %>%
  mutate(pathway = recode(pathway,
    IFNa1 = "IFN-α respuesta", IFNg1 = "IFN-γ respuesta",
    KRAS1 = "KRAS signaling"))
long$pathway <- factor(long$pathway,
  levels = c("IFN-α respuesta", "IFN-γ respuesta", "KRAS signaling"))
long$status <- factor(long$status, levels = c("Diagnóstico", "Refractario"))

# Medianas + etiqueta desplazada POR ENCIMA de la caja (fuera del violín)
meds <- long %>%
  group_by(pathway, status) %>%
  summarise(med = median(score), .groups = "drop") %>%
  group_by(pathway) %>%
  mutate(rng = diff(range(long$score[long$pathway == pathway])),
         y_label = med + 0.18 * rng) %>%
  ungroup()

pvals <- long %>% group_by(pathway) %>%
  summarise(p = wilcox.test(score[status=="Diagnóstico"],
                            score[status=="Refractario"])$p.value, .groups = "drop") %>%
  mutate(label = paste0("p = ", formatC(p, format = "e", digits = 2)))

pal <- c("Diagnóstico" = "#3B82BF", "Refractario" = "#E64B35")

g <- ggplot(long, aes(x = status, y = score, fill = status)) +
  geom_violin(trim = TRUE, alpha = 0.35, color = NA, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9,
               color = "grey20", fill = "white") +
  # mediana: rombo sobre la caja
  geom_point(data = meds, aes(x = status, y = med, fill = status),
             shape = 23, size = 3, color = "black", stroke = 0.7) +
  # etiqueta de mediana POR ENCIMA, con fondo blanco para que no se confunda
  geom_label(data = meds, aes(x = status, y = y_label, label = sprintf("%.3f", med)),
             fill = "white", color = "black", size = 3.5, fontface = "bold",
             label.size = 0.25, label.padding = unit(0.15, "lines"),
             inherit.aes = FALSE) +
  # p-valor arriba del todo
  geom_text(data = pvals, aes(x = 1.5, y = Inf, label = label),
            vjust = 1.2, size = 3.4, inherit.aes = FALSE, fontface = "italic") +
  facet_wrap(~ pathway, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Module score (z)", fill = NULL,
       title = "Vías en células B tumorales: Diagnóstico vs Refractario",
       subtitle = "Violín + caja. Mediana: rombo ▾ con valor encima. n = 8.521 Dx / 433 Ref") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10))

ggsave("figures/05_pathway_scores_dx_vs_ref.png", g, width = 11, height = 4.5, dpi = 160)
cat("Figura guardada.\n")
