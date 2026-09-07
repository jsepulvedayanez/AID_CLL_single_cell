suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
  library(DESeq2)
})

set.seed(42)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---------- Load & join ----------
merged <- readRDS("data/merged_MV12.rds")
merged <- JoinLayers(merged, assay = "RNA", layers = c("counts", "data"))
Idents(merged) <- "RNA_snn_res.0.3"

# ---------- 1. Corrected annotation (all 12 clusters) ----------
# Evidence: demux labels + canonical markers (inspect_markers.R output).
cluster_annotation <- c(
  "0"  = "CLL_B",          # CD79A++ FCER2(CD23)+ CD200+  -> S457 dx bulk
  "1"  = "CD8_T",          # CD3D/E+ CD8A/B+ NKG7+ GZMB+   -> S477 ref dominant
  "2"  = "CLL_B_plasma",   # CD79A++ XBP1+ PRDM1+ TNFRSF17+ CD38+ CD200+ (pre-plasmablastic CLL)
  "3"  = "CLL_B",          # CD79A+ CD19+ FCER2+ CD5+      -> S477 dx
  "4"  = "CD4_T",          # CD3D/E+ CD4+ GZMK+            -> ref dominant
  "5"  = "CLL_B",          # MS4A1++ CD79A+ CD79B+         -> S477 dx
  "6"  = "CLL_B_prolif",   # B markers + MKI67/TOP2A/PCNA/STMN1 + highest AICDA -> proliferative CLL
  "7"  = "CLL_B",          # CD79A+ FCER2+ ROR1+           -> S457 dx (low-intensity B)
  "8"  = "NK",             # NKG7++ GNLY++ KLRD1+ GZMB+ FCGR3A+
  "9"  = "Monocyte",       # CD14+ LYZ+ S100A8/A9+ FCGR3A+
  "10" = "CLL_B",          # MS4A1++ CD79A+ CD79B+         -> S477
  "11" = "Mega_NK"         # megakaryocyte (46) + NK (20) mixed
)
merged$cell_type <- unname(cluster_annotation[as.character(merged$RNA_snn_res.0.3)])

# Tumor compartment = B cells (demux) — clone fraction 80-96% so B≈tumor
merged$is_Bcell    <- merged$celltype_demux == "Bcell"
merged$is_tumor    <- merged$cell_type %in% c("CLL_B","CLL_B_plasma","CLL_B_prolif")
merged$tumor_final <- merged$is_tumor & merged$is_Bcell

# ---------- 2. AID definition (data-driven) ----------
# AICDA is ultra-sparse (20 cells). Two tiers:
#  (a) AID+ : any detected AICDA (counts/data > 0)
#  (b) AID-program module: AICDA + CSR/SHM machinery
aicda <- FetchData(merged, vars = "AICDA")
merged$AICDA_expr  <- aicda$AICDA
merged$is_AID_positive <- aicda$AICDA > 0

aid_module <- c("AICDA","UNG","APOBEC3A","APOBEC3B","APOBEC3C","APOBEC3G",
                "MSH2","MSH6","EXO1","POLH","PMS2","MLH1")
aid_module <- aid_module[aid_module %in% rownames(merged)]
merged <- AddModuleScore(merged, features = list(aid_module), name = "AID_program")

# ---------- 3. Proliferation definition (data-driven) ----------
# Combine cell-cycle phase + MKI67 + a proliferation module (no hardcoded cluster).
prolif_genes <- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","STMN1","TYMS","HIST1H4C")
prolif_genes <- prolif_genes[prolif_genes %in% rownames(merged)]
merged <- AddModuleScore(merged, features = list(prolif_genes), name = "Prolif")
merged$prolif_score <- merged$S.Score + merged$G2M.Score
merged$is_proliferating <- merged$Phase %in% c("S", "G2M") & merged$is_tumor

# ---------- 4. Q3 pathway module scores (IFN-a, IFN-g, KRAS) ----------
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
get_genes <- function(nm) unique(hallmark$gene_symbol[hallmark$gs_name == nm])

merged <- AddModuleScore(merged, features = list(get_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE")), name = "IFNa")
merged <- AddModuleScore(merged, features = list(get_genes("HALLMARK_INTERFERON_GAMMA_RESPONSE")), name = "IFNg")
merged <- AddModuleScore(merged, features = list(get_genes("HALLMARK_KRAS_SIGNALING_UP")),     name = "KRAS")

# ---------- Summaries ----------
cat("========== ANNOTATION ==========\n")
print(table(cell_type = merged$cell_type, status = merged$status))
cat("\nTumor (B) cells:", sum(merged$tumor_final), "of", ncol(merged), "\n")
cat("Tumor cells by status:\n"); print(table(merged$status[merged$tumor_final]))

cat("\n========== Q1: AID ==========\n")
cat("AICDA>0 cells (all):", sum(merged$is_AID_positive), "=", round(100*mean(merged$is_AID_positive),3), "%\n")
aid_by_type <- table(cell_type = merged$cell_type, AID = merged$is_AID_positive)
print(aid_by_type)
cat("\nAICDA>0 among tumor B cells:", sum(merged$is_AID_positive & merged$tumor_final), "\n")
cat("AICDA>0 among tumor B by status:\n")
print(table(merged$status[merged$tumor_final & merged$is_AID_positive]))
# AID program score by cell type (top)
aid_score_by_type <- tapply(merged$AID_program1, merged$cell_type, mean)
cat("\nAID_program mean score by cell_type:\n")
print(round(sort(aid_score_by_type, decreasing = TRUE), 4))

cat("\n========== Q2: PROLIFERATION ==========\n")
cat("Phase distribution (all cells):\n"); print(table(merged$Phase))
cat("\nPhase distribution (tumor B):\n"); print(table(merged$Phase[merged$tumor_final]))
cat("\nProliferating tumor cells (Phase S/G2M):", sum(merged$is_proliferating),
    "=", round(100*sum(merged$is_proliferating)/sum(merged$tumor_final),2), "% of tumor\n")
cat("MKI67>0 cells:", sum(merged$is_proliferating | (FetchData(merged, vars='MKI67')$MKI67>0 & merged$tumor_final)), "\n")
cat("Proliferating tumor by status:\n")
print(table(merged$status[merged$is_proliferating]))
cat("\nProliferating tumor by cluster:\n")
print(table(merged$cell_type[merged$is_proliferating]))

# ---------- Save annotation metadata ----------
meta_out <- merged@meta.data %>%
  select(subject, sample_id, status, chip, celltype_demux, cell_type,
         is_tumor, tumor_final, is_AID_positive, AICDA_expr, AID_program1,
         is_proliferating, prolif_score, S.Score, G2M.Score, Phase,
         IFNa1, IFNg1, KRAS1)
write.csv(meta_out, "results/cell_annotation_metadata.csv", row.names = TRUE)

saveRDS(merged, "results/CLL_scRNAseq_merged_annotated.rds")
cat("\nSaved annotated object + metadata CSV.\n")
