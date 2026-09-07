# AID-Positive CLL Single-Cell Analysis

## Overview

This repository contains the analysis code for paired single-cell RNA-seq (scRNA-seq) and single-cell B-cell receptor sequencing (scBCR-seq) data from two chronic lymphocytic leukemia (CLL) patients.

The study focuses on AID-positive disease and compares diagnosis versus refractory samples to characterize cell states, tumor composition, clonality, and pathway activity in peripheral blood.

## Study design

### Patients

This project includes two CLL patients:

* **Patient 457**
* **Patient 477**

Both patients are:

* Male
* Diagnosed with **Binet C** CLL
* **IGHV unmutated**
* AID-positive in a subpopulation of cells
* Sampled from **peripheral blood**

### Clinical and molecular context

| Patient | Diagnosis date | Treatment           | Isotype | IGHV status | IGH                          | IGL               |           % tumor clone | Observations     |
| ------- | -------------- | ------------------- | ------- | ----------- | ---------------------------- | ----------------- | ----------------------: | ---------------- |
| 457     | 2017           | Ibrutinib           | IgM     | Unmutated   | VH1-69*01 / D3-3*01 / J6*06  | Lambda + d (weak) | 96.1 (diagnosis sample) | AID (+), LPL (+) |
| 477     | 08/2018        | Ibrutinib (05/2020) | IgG     | Unmutated   | VH4-39*01 / D6-19*01 / J5*02 | Kappa             |   80 (diagnosis sample) | AID (+), LPL (+) |

## Samples

Each patient has paired samples at two time points:

* **Diagnosis**
* **Refractory**

### 10x Genomics experiment/chip information

| Patient | Time point | ID   | 10x chip / experiment |
| ------- | ---------- | ---- | --------------------- |
| 457     | Diagnosis  | SID2 | MV1                   |
| 457     | Refractory | SID2 | MV2                   |
| 477     | Diagnosis  | SID1 | MV2                   |
| 477     | Refractory | SID1 | MV1                   |


## Data types

* **scRNA-seq**: transcriptomic profiling of individual cells
* **scBCR-seq**: B-cell receptor profiling to assess clonotype structure and clonal relationships

## Main questions
1. Which cells express AID-related programs?
2. Is there a proliferative tumor-cell subpopulation?
3. Are IFN-\u03b1, IFN-\u03b3, and KRAS pathways enriched in specific tumor cells?

## Analysis workflow

1. Quality control and filtering
2. Normalization and integration
3. Batch assessment across chips (MV1 / MV2)
4. Cell clustering and annotation
5. Identification of B cells and tumor cells
6. Detection of AID-positive subpopulations
7. Proliferation scoring and filtering, if applicable
8. Pathway enrichment analysis


## Important notes

* The dataset comes from **blood**, so it includes more than only B cells.
* Diagnosis and refractory samples are expected to have strong **cell-composition imbalance**, with diagnosis enriched for B cells and refractory samples containing a more mixed blood-cell compartment.
* The B-cell fraction is not equivalent to the tumor-cell fraction; tumor cells still need to be defined using clonotype and expression-based annotation.
* The **proliferative fraction still needs to be defined by analysis**.
* The chip/experiment information is important and should be considered during batch correction and interpretation.
* Because samples are paired by patient and time point, comparisons should account for both **biological state** and **technical batch**.

## Repository structure

```text
README.md
data/                            # datos crudos (matrices 10x .h5, Seurat .rds)
notebooks/
  preproccesing.ipynb             # QC, integración, clustering (Seurat v5)
  cell_annotation.ipynb           # anotación + Q1/Q2/Q3 (corregido y ejecutado)
  cell_annotation_original.ipynb  # versión previa (anotación con errores)
scripts/
  optimized_annotation.R          # anotación corregida (12 clusters)
  Q3_dx_vs_ref_pooled.R           # análisis de vías Dx vs Ref (pooled)
  Q3_DE_filtered_and_kegg.R       # DE filtrado (IG/TCR/mito/ribo/hb) + GSEA KEGG
  make_figures.R                  # figuras base (UMAP, AICDA, MKI67)
  make_figures_subset.R           # figuras 05 (completo) y 06 (subset 433 vs 433)
  _build_notebook.py              # construye cell_annotation.ipynb (nbformat)
results/                          # salidas (CSV + RDS anotado, no se trackea .rds)
figures/                          # figuras PNG (05 completo, 06 subset justo)
```

## Results (resumen)

Análisis ejecutado en `notebooks/cell_annotation.ipynb`. Células B tumorales = demux "Bcell" + ≥1 marcador B (MS4A1/CD79A/CD19) → **8.954 células** (8.521 diagnóstico, 433 refractario).

1. **AID:** AICDA detectable en 20/14.550 células (0,14%); 17 son B tumorales, con el programa AID enriquecido en el cluster proliferativo (CLL_B_prolif).
2. **Proliferación:** subpoblación proliferativa (cluster 6, MKI67/TOP2A/PCNA/STMN1); fracción MKI67+ ≈0,2% (CLL quiescente).
3. **Vías (Dx vs Ref, pooled, solo B tumorales):** con comparación justa (submuestreo 433 vs 433 estratificado por sujeto, 20 semillas): **IFN-α/IFN-γ ↑ en refractario** (Cliff δ ≈ −0,41/−0,39), **KRAS ↑ en diagnóstico** (δ ≈ +0,14).
4. **Enriquecimiento (GSEA):** diagnóstico muestra un programa **metabólico/proliferativo** — OXPHOS (NES 1.61), glicólisis, mTORC1, MYC, y **RIBOSOMA** (KEGG, NES 1.92, padj 1e-21).
5. **Pseudobulk DESeq2 (n=4):** tras filtrar IG/TCR/mito/ribo/hb, 19 genes DE (padj<0.05). ↑ diagnóstico: JUN, TSC22D3, OTUD1; ↑ refractario: CDK6, SOX4, HCK, CD9, MARCKS, GAS6 (más GNLY/GZMB = contaminación NK residual). Baja confianza por n=4 y contaminación clonal.

## Suggested output folders

```text
results/qc/
results/integration/
results/clonotypes/
results/proliferation/
results/pathways/
figures/
```

## Status

This project is in progress and is being organized for analysis and interpretation of paired scRNA-seq and scBCR-seq data in AID-positive CLL.
