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
7. Detection of AID-positive subpopulations
8. Proliferation scoring and filtering, if applicable
9. Pathway enrichment analysis

## Important notes

* The dataset comes from **blood**, so it includes more than only B cells.
* The **proliferative fraction still needs to be defined by analysis**.
* The chip/experiment information is important and should be considered during batch correction and interpretation.
* Because samples are paired by patient and time point, comparisons should account for both **biological state** and **technical batch**.

## Suggested repository structure

```text
README.md
data/
metadata/
scripts/
notebooks/
results/
figures/
```

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

