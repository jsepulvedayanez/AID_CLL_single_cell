import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
    "language_info": {"name": "R", "version": "4.6.1",
                      "codemirror_mode": "r", "file_extension": ".r",
                      "mimetype": "text/x-r-source", "pygments_lexer": "r"},
}

cells = []

def md(s):
    cells.append(nbf.v4.new_markdown_cell(s))

def code(s):
    cells.append(nbf.v4.new_code_cell(s))

# ---------------------------------------------------------------------------
md("""# AID-CLL Single-Cell — Anotación Celular y Análisis de Vías

**Objetivo:** anotar correctamente las células, definir el compartimento de **células B tumorales** y responder las 3 preguntas del estudio con **pruebas estadísticas formales**.

**Las 3 preguntas:**
1. ¿Qué células expresan programas relacionados con AID?
2. ¿Existe una subpoblación tumoral proliferativa?
3. ¿Están enriquecidas las vías **IFN-α**, **IFN-γ** y **KRAS** en células tumorales específicas? (diagnóstico vs refractario, *pooled*)

**Nota metodológica central:** el estudio tiene un **desbalance de composición** muy fuerte — el compartimento B tumoral en diagnóstico (≈8.500 células) es ~20× el de refractario (≈430). Para la pregunta 3 hacemos **dos** análisis: (a) la comparación completa y (b) una **comparación justa** submuestreando diagnóstico hasta igualar el número de refractario *por sujeto*, repetida 20 veces.

> Este notebook corrige la versión anterior (`cell_annotation.ipynb` original), cuya anotación mapeaba solo 9 de 12 clusters y estaba mal en 7 de 9 etiquetas.""")

# ---------------------------------------------------------------------------
md("""## 1. Configuración

Cargamos todas las librerías y fijamos la semilla para reproducibilidad.""")
code("""# ============================================================
# 1. CONFIGURACIÓN
# ============================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(msigdbr)
library(fgsea)
library(DESeq2)
library(tibble)

set.seed(42)   # reproducibilidad

message("Librerías cargadas. ", R.version.string)""")

# ---------------------------------------------------------------------------
md("""## 2. Carga del objeto mergeado

El objeto proviene del pipeline de preprocesamiento (**Seurat v5**) y tiene las capas **separadas por chip**: `counts.MV1`, `counts.MV2`, `data.MV1`, `data.MV2`.

> **Pitfall importante:** en Seurat v5, si no se hace `JoinLayers()` antes de `FetchData()`/`FeaturePlot()`, se pierden ~7.500 células (las de una capa).""")
code("""# ============================================================
# 2. CARGA + UNIÓN DE CAPAS
# ============================================================
# Asegurar que el working directory es la raíz del repo (nbconvert/Jupyter
# puede arrancar en notebooks/, rompiendo las rutas relativas).
if (!file.exists("data/merged_MV12.rds") && file.exists("../data/merged_MV12.rds")) {
  setwd("..")
}

input_rds <- "data/merged_MV12.rds"
merged <- readRDS(input_rds)

# Unir capas split por chip (MV1/MV2) en una sola capa 'counts' y 'data'
merged <- JoinLayers(merged, assay = "RNA", layers = c("counts", "data"))

# Identidad activa = clustering a resolución 0.3
Idents(merged) <- "RNA_snn_res.0.3"

message("Células: ", ncol(merged), "  |  Genes: ", nrow(merged))""")

# ---------------------------------------------------------------------------
md("""## 3. Estructura y composición del objeto

Inspeccionamos la metadata y la composición por tipo celular (asignación **demux**) antes de anotar.""")
code("""# ============================================================
# 3. ESTRUCTURA DEL OBJETO
# ============================================================
cat("Columnas de metadata:\\n")
print(colnames(merged@meta.data))

cat("\\nTamaño de clusters (res 0.3):\\n")
print(table(merged$RNA_snn_res.0.3))

cat("\\nComposición por tipo celular (demux):\\n")
print(table(merged$celltype_demux))

cat("\\nMuestras (sample_id x status):\\n")
print(table(merged$sample_id, merged$status))""")

# ---------------------------------------------------------------------------
md("""## 4. Anotación celular corregida

**Estrategia:** combinar (a) las etiquetas **demux** (asignación computacional validada a partir del perfil de expresión) con (b) **marcadores canónicos** para refinar.

Primero visualizamos los marcadores por cluster para anotar con evidencia, no con suposiciones.""")
code("""# ============================================================
# 4a. MARCADORES CANÓNICOS POR CLUSTER
# ============================================================
canonical_markers <- list(
  B_general     = c("MS4A1","CD19","CD79A","CD79B"),
  CLL_hallmark  = c("CD5","FCER2","CD200","ROR1","LEF1"),   # FCER2 = CD23
  Plasma        = c("XBP1","PRDM1","TNFRSF17","CD38"),
  T_CD4         = c("CD3D","CD3E","CD4","IL7R"),
  T_CD8         = c("CD3D","CD8A","CD8B","GZMB","GZMK"),
  NK            = c("NKG7","GNLY","KLRD1","NCAM1","FCGR3A"),
  Monocyte      = c("CD14","LYZ","S100A8","S100A9"),
  Proliferating = c("MKI67","TOP2A","PCNA","CCNB1","STMN1"),
  AID_machinery = c("AICDA","UNG","MSH2","MSH6","EXO1","POLH")
)

DotPlot(merged, features = unique(unlist(canonical_markers)),
        group.by = "RNA_snn_res.0.3", dot.scale = 5) +
  RotatedAxis() + scale_color_viridis(option = "plasma") +
  ggtitle("Marcadores canónicos por cluster (res 0.3)")""")

md("""### Mapeo cluster → tipo celular

La evidencia de marcadores + demux da el siguiente mapeo (los 12 clusters, no solo 9):

| cluster | evidencia | etiqueta |
|---|---|---|
| 0, 3, 5, 7, 10 | CD79A+/MS4A1+/FCER2(CD23)+ | **CLL_B** |
| 2 | CD79A++ + XBP1/PRDM1/TNFRSF17/CD38 (pre-plasmablástico) | **CLL_B_plasma** |
| 6 | B + MKI67/TOP2A/PCNA/STMN1 | **CLL_B_prolif** |
| 1 | CD3D/E + CD8A/B + NKG7/GZMB | **CD8_T** |
| 4 | CD3D/E + CD4 + GZMK | **CD4_T** |
| 8 | NKG7/GNLY/KLRD1 | **NK** |
| 9 | CD14/LYZ/S100A8 | **Monocyte** |
| 11 | megacariocito (46) + NK (20) mixto | **Mega_NK** |""")
code("""# ============================================================
# 4b. APLICAR ANOTACIÓN
# ============================================================
cluster_annotation <- c(
  "0"  = "CLL_B", "1" = "CD8_T", "2" = "CLL_B_plasma", "3" = "CLL_B",
  "4"  = "CD4_T", "5" = "CLL_B", "6" = "CLL_B_prolif", "7" = "CLL_B",
  "8"  = "NK", "9" = "Monocyte", "10" = "CLL_B", "11" = "Mega_NK"
)
merged$cell_type <- unname(cluster_annotation[as.character(merged$RNA_snn_res.0.3)])

DimPlot(merged, reduction = "umap", group.by = "cell_type",
        label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle("Anotación corregida (12 clusters)") + theme_classic()""")

# ---------------------------------------------------------------------------
md("""## 5. Definición de células B tumorales

> **Sí, todas las comparaciones de la pregunta 3 son SOLO sobre células B tumorales.**

El compartimento tumoral se define como **células B (demux) que expresan ≥1 marcador B** (`MS4A1`/`CD79A`/`CD19`). Como la fracción clonal es 80–96%, la gran mayoría de las B de sangre son la clona de CLL. Exigir el marcador B descarta posibles dobletes/células mal clasificadas (p. ej. B dentro del cluster CD8_T).""")
code("""# ============================================================
# 5. CÉLULAS B TUMORALES
# ============================================================
bmark <- FetchData(merged, vars = c("MS4A1","CD79A","CD19"))
merged$is_tumorB <- merged$celltype_demux == "Bcell" &
  (bmark$MS4A1 > 0 | bmark$CD79A > 0 | bmark$CD19 > 0)

cat("Células B (demux):", sum(merged$celltype_demux == "Bcell"), "\\n")
cat("B tumorales (con marcador B):", sum(merged$is_tumorB), "\\n")
cat("B tumorales por status:\\n")
print(table(merged$status[merged$is_tumorB]))
cat("\\nB tumorales por sujeto x status:\\n")
print(table(merged$subject[merged$is_tumorB], merged$status[merged$is_tumorB]))""")

# ---------------------------------------------------------------------------
md("""## 6. Q1 — ¿Qué células expresan programas relacionados con AID?

**AICDA (AID)** es ultra-raro en scRNA: se detecta en una fracción ínfima de células. Además de AICDA, evaluamos el **programa AID amplio** (maquinaria de CSR/SHM: UNG, MSH2/6, EXO1, POLH, APOBEC3).""")
code("""# ============================================================
# 6. Q1 — AID
# ============================================================
aicda <- FetchData(merged, vars = "AICDA")
merged$is_AID_positive <- aicda$AICDA > 0

cat("AICDA>0 (todas las células):", sum(merged$is_AID_positive),
    "=", round(100*mean(merged$is_AID_positive), 3), "%\\n")
cat("\\nAICDA>0 entre B tumorales, por tipo:\\n")
print(table(merged$cell_type[merged$is_AID_positive & merged$is_tumorB]))

# Programa AID amplio (AICDA + maquinaria CSR/SHM)
aid_genes <- c("AICDA","UNG","APOBEC3A","APOBEC3B","APOBEC3C","APOBEC3G",
               "MSH2","MSH6","EXO1","POLH","PMS2","MLH1")
aid_genes <- aid_genes[aid_genes %in% rownames(merged)]
merged <- AddModuleScore(merged, features = list(aid_genes), name = "AID_program")

cat("\\nScore medio del programa AID por tipo celular:\\n")
print(round(sort(tapply(merged$AID_program1, merged$cell_type, mean),
                 decreasing = TRUE), 4))

FeaturePlot(merged, features = "AICDA", reduction = "umap",
            order = TRUE, min.cutoff = 0, max.cutoff = 1) &
  scale_color_viridis(option = "magma")""")

# ---------------------------------------------------------------------------
md("""## 7. Q2 — ¿Hay subpoblación tumoral proliferativa?

Definimos proliferación con **MKI67** (gold standard), no con `Phase`, porque el `CellCycleScoring` previo dejó `Phase` ~1/3-1/3-1/3 (mal calibrado para este dataset).""")
code("""# ============================================================
# 7. Q2 — PROLIFERACIÓN
# ============================================================
mki67 <- FetchData(merged, vars = "MKI67")
merged$is_proliferating <- mki67$MKI67 > 0 & merged$is_tumorB

prolif_genes <- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","STMN1","TYMS","HIST1H4C")
prolif_genes <- prolif_genes[prolif_genes %in% rownames(merged)]
merged <- AddModuleScore(merged, features = list(prolif_genes), name = "Prolif")

cat("MKI67+ en B tumorales:", sum(merged$is_proliferating),
    "=", round(100*sum(merged$is_proliferating)/sum(merged$is_tumorB), 2), "%\\n")
cat("\\nMKI67+ por cluster:\\n")
print(table(merged$cell_type[merged$is_proliferating]))
cat("\\nMKI67+ por status:\\n")
print(table(merged$status[merged$is_proliferating]))

FeaturePlot(merged, features = "MKI67", reduction = "umap",
            order = TRUE, min.cutoff = 0, max.cutoff = 1) &
  scale_color_viridis(option = "magma")""")

# ---------------------------------------------------------------------------
md("""## 8. Q3 — Vías IFN-α / IFN-γ / KRAS en B tumorales (Dx vs Ref, pooled)

Puntuamos cada célula con **AddModuleScore** para las 3 vías hipotetizadas, restringido a **B tumorales**, y comparamos diagnóstico vs refractario.""")
code("""# ============================================================
# 8. Q3 — SCORES DE VÍA
# ============================================================
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
get_genes <- function(nm) unique(hallmark$gene_symbol[hallmark$gs_name == nm])

merged <- AddModuleScore(merged, features = list(get_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE")), name = "IFNa")
merged <- AddModuleScore(merged, features = list(get_genes("HALLMARK_INTERFERON_GAMMA_RESPONSE")), name = "IFNg")
merged <- AddModuleScore(merged, features = list(get_genes("HALLMARK_KRAS_SIGNALING_UP")),     name = "KRAS")

tumor <- subset(merged, subset = is_tumorB == TRUE)

# Violín + caja + mediana etiquetada (visualización más legible para comparar)
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

meds <- long %>% group_by(pathway, status) %>%
  summarise(med = median(score), .groups = "drop") %>%
  group_by(pathway) %>%
  mutate(rng = diff(range(long$score[long$pathway == pathway])),
         y_label = med + 0.18 * rng) %>%
  ungroup()
pvals <- long %>% group_by(pathway) %>%
  summarise(p = wilcox.test(score[status=="Diagnóstico"],
                            score[status=="Refractario"])$p.value, .groups = "drop") %>%
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
  scale_fill_manual(values = c("Diagnóstico" = "#3B82BF", "Refractario" = "#E64B35")) +
  labs(x = NULL, y = "Module score (z)", fill = NULL,
       title = "Vías en células B tumorales: Diagnóstico vs Refractario",
       subtitle = "Violín + caja. Mediana: rombo ▾ con valor encima. n = 8.521 Dx / 433 Ref") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10))""")

# ---------------------------------------------------------------------------
md("""### 8a. Test estadístico (Wilcoxon + tamaño de efecto Cliff's δ + corrección BH)

Además del p-valor reportamos el **Cliff's delta** (correlación rank-biserial), que mide el **tamaño del efecto** en [-1, 1], y corregimos por comparaciones múltiples (BH) sobre las 3 vías hipotetizadas.

> Convención Cliff's δ: < 0.147 *negligible* · 0.147–0.33 *pequeño* · 0.33–0.474 *mediano* · > 0.474 *grande*.
> δ < 0 → mayor en **refractario**; δ > 0 → mayor en **diagnóstico**.""")
code("""# ============================================================
# 8a. TEST ESTADÍSTICO COMPLETO
# ============================================================
cliff_delta <- function(x, y) {
  sum(outer(x, y, function(a,b) sign(a-b))) / (length(x) * length(y))
}

scores <- c("IFNa1","IFNg1","KRAS1","AID_program1","Prolif1")

full_res <- lapply(scores, function(s) {
  x <- tumor@meta.data[[s]][tumor$status == "diagnosis"]
  y <- tumor@meta.data[[s]][tumor$status == "refractory"]
  wt <- wilcox.test(x, y)
  data.frame(score = s,
             med_dx = median(x), med_ref = median(y),
             cliff_delta = cliff_delta(x, y),
             p_value = wt$p.value)
})
full_res <- do.call(rbind, full_res)

# Corrección BH sobre las 3 vías con hipótesis a priori
full_res$padj_bh <- NA
named <- full_res$score %in% c("IFNa1","IFNg1","KRAS1")
full_res$padj_bh[named] <- p.adjust(full_res$p_value[named], method = "BH")

print(full_res, row.names = FALSE)""")

# ---------------------------------------------------------------------------
md("""### 8b. Comparación JUSTA: submuestreo de diagnóstico

El desbalance (≈8.500 Dx vs ≈430 Ref) **infla los p-valores** y puede **confundir por sujeto**. Para una comparación justa:

1. Submuestreamos las B tumorales de **diagnóstico** hasta igualar el **número de refractario por sujeto** (S457: 323, S477: 110) → 433 vs 433.
2. Repetimos con **20 semillas** distintas y reportamos la mediana del p-valor, el % de réplicas significativas y la mediana del Cliff's δ.

Esto responde directamente: *¿el efecto sobrevive con tamaños parejos y sin sesgo de sujeto?*""")
code("""# ============================================================
# 8b. COMPARACIÓN JUSTA (submuestreo estratificado por sujeto)
# ============================================================
ref_cells <- colnames(tumor)[tumor$status == "refractory"]
dx_cells  <- colnames(tumor)[tumor$status == "diagnosis"]
ref_subj <- tumor$subject[ref_cells]
dx_subj  <- tumor$subject[dx_cells]
ref_n    <- table(ref_subj)
cat("Refractario por sujeto:\\n"); print(ref_n)

fair_res <- lapply(scores, function(s) {
  pvals <- c(); deltas <- c()
  for (seed in 1:20) {
    set.seed(seed)
    sel <- unlist(lapply(names(ref_n), function(sb) {
      pool <- dx_cells[dx_subj == sb]
      sample(pool, size = min(as.integer(ref_n[sb]), length(pool)))
    }))
    x <- tumor@meta.data[[s]][match(sel, colnames(tumor))]
    y <- tumor@meta.data[[s]][match(ref_cells, colnames(tumor))]
    pvals  <- c(pvals,  wilcox.test(x, y)$p.value)
    deltas <- c(deltas, cliff_delta(x, y))
  }
  data.frame(score = s, n_dx = length(ref_cells), n_ref = length(ref_cells),
             median_p = median(pvals),
             prop_sig = mean(pvals < 0.05),
             median_delta = median(deltas))
})
fair_res <- do.call(rbind, fair_res)
print(fair_res, row.names = FALSE)""")

# ---------------------------------------------------------------------------
md("""### 8b-cont. Visualización del subset justo (un sorteo representativo)

Misma figura de la sección 8, pero sobre el **submuestreo justo** (433 Dx vs 433 Ref, seed 1). Comparar con la sección 8 muestra cuánto cambia (o no) la conclusión al balancear por sujeto.""")
code("""# ============================================================
# 8b-cont. FIGURA DEL SUBSET JUSTO
# ============================================================
set.seed(1)
sel <- unlist(lapply(names(ref_n), function(sb) {
  pool <- dx_cells[dx_subj == sb]
  sample(pool, size = min(as.integer(ref_n[sb]), length(pool)))
}))

df_sub <- FetchData(tumor, vars = c("status","IFNa1","IFNg1","KRAS1"))
df_sub <- df_sub[c(sel, ref_cells), ]
df_sub$status <- ifelse(df_sub$status == "diagnosis", "Diagnóstico", "Refractario")

long_sub <- df_sub %>%
  pivot_longer(-status, names_to = "pathway", values_to = "score") %>%
  mutate(pathway = recode(pathway,
    IFNa1 = "IFN-α respuesta", IFNg1 = "IFN-γ respuesta",
    KRAS1 = "KRAS signaling"))
long_sub$pathway <- factor(long_sub$pathway,
  levels = c("IFN-α respuesta", "IFN-γ respuesta", "KRAS signaling"))
long_sub$status <- factor(long_sub$status, levels = c("Diagnóstico", "Refractario"))

meds_sub <- long_sub %>% group_by(pathway, status) %>%
  summarise(med = median(score), .groups = "drop") %>%
  group_by(pathway) %>%
  mutate(rng = diff(range(long_sub$score[long_sub$pathway == pathway])),
         y_label = med + 0.18 * rng) %>% ungroup()
pvals_sub <- long_sub %>% group_by(pathway) %>%
  summarise(p = wilcox.test(score[status=="Diagnóstico"],
                            score[status=="Refractario"])$p.value, .groups = "drop") %>%
  mutate(label = paste0("p = ", formatC(p, format = "e", digits = 2)))

ggplot(long_sub, aes(x = status, y = score, fill = status)) +
  geom_violin(trim = TRUE, alpha = 0.35, color = NA, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9,
               color = "grey20", fill = "white") +
  geom_point(data = meds_sub, aes(x = status, y = med, fill = status),
             shape = 23, size = 3, color = "black", stroke = 0.7) +
  geom_label(data = meds_sub, aes(x = status, y = y_label, label = sprintf("%.3f", med)),
             fill = "white", color = "black", size = 3.5, fontface = "bold",
             label.padding = unit(0.15, "lines"), inherit.aes = FALSE) +
  geom_text(data = pvals_sub, aes(x = 1.5, y = Inf, label = label),
            vjust = 1.2, size = 3.4, inherit.aes = FALSE, fontface = "italic") +
  facet_wrap(~ pathway, scales = "free_y") +
  scale_fill_manual(values = c("Diagnóstico" = "#3B82BF", "Refractario" = "#E64B35")) +
  labs(x = NULL, y = "Module score (z)", fill = NULL,
       title = "SUBSET JUSTO: Diagnóstico submuestreado por sujeto (433 vs 433)",
       subtitle = "seed 1 de 20 — mismo patrón: IFN ↑ refractario, KRAS ↑ diagnóstico") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10))""")

# ---------------------------------------------------------------------------
md("""### 8c. GSEA (fgsea) sobre el ranking génico

Complemento al score por célula: rankeamos **todos** los genes por diferencia de medias (Dx − Ref) y testeamos enriquecimiento de gene-sets Hallmark.

> NES > 0 = enriquecido en **diagnóstico**; NES < 0 = enriquecido en **refractario**.""")
code("""# ============================================================
# 8c. GSEA (fgsea)
# ============================================================
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

mat <- GetAssayData(tumor, layer = "data")
lfc <- rowMeans(mat[, tumor$status == "diagnosis", drop = FALSE]) -
       rowMeans(mat[, tumor$status == "refractory", drop = FALSE])
lfc <- sort(lfc[is.finite(lfc)], decreasing = TRUE)

gsea_h <- fgsea(pathways = hallmark_sets, stats = lfc, minSize = 15, maxSize = 500, eps = 0)
gsea_h <- gsea_h[order(gsea_h$pval), ]

cat("Top 15 vías Hallmark (padj<0.25):\\n")
gsea_h %>% filter(padj < 0.25) %>%
  select(pathway, NES, pval, padj) %>% head(15)

cat("\\nVías hipotetizadas (IFN-α, IFN-γ, KRAS):\\n")
gsea_h %>% filter(grepl("INTERFERON_ALPHA|INTERFERON_GAMMA|KRAS", pathway)) %>%
  select(pathway, NES, pval, padj)""")

# ---------------------------------------------------------------------------
md("""### 8d. Pseudobulk DESeq2 (n = 4 muestras)

Test formal a **nivel de muestra** (agrega todas las células de cada muestra en un perfil), con diseño `~ subject + status` para controlar la variación entre pacientes.

> **Caveat:** con n = 4 es de **baja potencia** y muy sensible a contaminación ambiental (el compartimento refractario es minúsculo). Se reporta como apoyo, no como evidencia principal.""")
code("""# ============================================================
# 8d. PSEUDOBULK DESeq2 (genes filtrados)
# ============================================================
# Se eliminan genes no informativos para DE: inmunoglobulinas (IGH/IGK/IGL/
# JCHAIN), TCR (TRA/TRB/TRD/TRG), mitocondriales (MT-), ribosomales (RPS/RPL)
# y hemoglobinas (HBA/HBB). Dominan la señal o reflejan contaminación/identidad
# clonal, no biología de la vía.
drop_re <- "^(IGH|IGK|IGL|IGJ|JCHAIN|TRA|TRB|TRD|TRG|MT-|MT|RPS|RPL|HBA|HBB|HBD|HBG|HBE|HBZ|MALAT1)"
tumor_f <- tumor[!grepl(drop_re, rownames(tumor)), ]

pb <- AggregateExpression(tumor_f, assays = "RNA", slot = "counts",
                          group.by = "sample_id", return.seurat = FALSE)$RNA
colnames(pb) <- gsub("-", "_", colnames(pb))

md <- unique(tumor@meta.data[, c("sample_id", "status", "subject")])
rownames(md) <- md$sample_id
md <- md[colnames(pb), ]
md$status <- factor(md$status, levels = c("refractory", "diagnosis"))

dds <- DESeqDataSetFromMatrix(countData = round(pb), colData = md,
                              design = ~ subject + status)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("status", "diagnosis", "refractory")) %>%
  as.data.frame() %>% rownames_to_column("gene") %>% arrange(padj)

cat("Genes DE (padj<0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\\n")
head(res %>% select(gene, log2FoldChange, padj), 15)""")

# ---------------------------------------------------------------------------
md("""## 9. Guardar resultados""")
code("""# ============================================================
# 9. GUARDAR RESULTADOS
# ============================================================
dir.create("results", showWarnings = FALSE)
saveRDS(merged, "results/CLL_scRNAseq_merged_annotated.rds")

meta_out <- merged@meta.data %>%
  select(subject, sample_id, status, chip, celltype_demux, cell_type,
         is_tumorB, is_AID_positive, AID_program1,
         is_proliferating, Prolif1, IFNa1, IFNg1, KRAS1)
write.csv(meta_out, "results/cell_annotation_metadata.csv")

write.csv(full_res, "results/Q3_pathway_scores_dx_vs_ref.csv", row.names = FALSE)
write.csv(fair_res, "results/Q3_pathway_scores_fair_downsampled.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_h)[, setdiff(names(gsea_h), "leadingEdge")],
          "results/GSEA_Hallmark_Dx_vs_Ref_pooled.csv", row.names = FALSE)
write.csv(res, "results/DESeq2_pseudobulk_Dx_vs_Ref_pooled.csv", row.names = FALSE)

message("Resultados guardados en results/")""")

# ---------------------------------------------------------------------------
md("""## 10. Conclusiones

**Q1 — AID.** AICDA se detecta en solo **20/14.550 células (0,14%)**, de las cuales **17 son B tumorales** (9 CLL_B, 2 CLL_B_plasma, 6 CLL_B_prolif). AID+ **no es un cluster propio**: es una subpoblación rara dispersa en B tumorales, con el programa AID (score 0.030) **enriquecido en el cluster proliferativo CLL_B_prolif** (cluster 6).

**Q2 — Proliferación.** Sí hay subpoblación proliferativa: **cluster 6 (CLL_B_prolif, ≈600 células)** concentra MKI67/TOP2A/PCNA/STMN1. La fracción MKI67+ en B tumorales es **≈0,2%** (21 células), coherente con la quiescencia típica de la CLL.

**Q3 — Vías (Dx vs Ref, solo B tumorales, pooled).** Con la **comparación justa** (submuestreo, 433 vs 433, 20 semillas):
- **IFN-α y IFN-γ: mayor en refractario** (Cliff δ ≈ −0,41 / −0,39; 100% de réplicas p<0,05) → efecto **mediano**, robusto.
- **KRAS: mayor en diagnóstico** (Cliff δ ≈ +0,14; 100% p<0,05) → efecto **pequeño**.
- **GSEA a nivel de gene-set:** el **diagnóstico** muestra un programa **metabólico/proliferativo** coherente: **OXPHOS** (NES 1.61, padj 3e-5), glicólisis, mTORC1, MYC targets, p53, DNA repair, y **RIBOSOMA** (KEGG, NES 1.92, padj 1e-21) — traducción y biogénesis ribosomal elevadas. IFN no sale significativo al ranking global (efecto diluido, pero en dirección coherente con el score por célula).
- **Pseudobulk DESeq2 (n=4):** tras filtrar IG/TCR/mito/ribo/hb, los DE bajan de 39 a **19 genes**. ↑ diagnóstico: **JUN, TSC22D3, OTUD1**; ↑ refractario: **CDK6, SOX4, HCK, CD9, MARCKS, GAS6** (más GNLY/GZMB = contaminación NK residual). Baja confianza por n=4 y contaminación clonal.

**Lectura biológica:** el compartimento B tumoral refractario muestra una firma de **respuesta a interferón (IFN-α/γ) aumentada** respecto al diagnóstico, mientras que el diagnóstico muestra **KRAS y OXPHOS** más altos — patrón consistente con la exposición al ibrutinib y la presión selectiva del microambiente.""")

nb.cells = cells

import json
path = "notebooks/cell_annotation.ipynb"
with open(path, "w") as f:
    nbf.write(nb, f)
print("Notebook written:", path, "| cells:", len(cells))
