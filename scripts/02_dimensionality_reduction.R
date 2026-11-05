## 1. Packages and paths

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(batchelor)
  library(BiocSingular)
})

data_dir <- "data"
output_dir <- "results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## 2. Load single-cell data

spe_qc <- readRDS(file.path(data_dir, "spe_qc.rds"))

stopifnot(
  is(spe_qc, "SpatialExperiment"),
  "exprs" %in% assayNames(spe_qc),
  all(c("batch", "patient_id") %in% names(colData(spe_qc))),
  !anyNA(spe_qc$batch),
  !anyDuplicated(colnames(spe_qc))
)

use_markers <- rownames(spe_qc)[
  !grepl("DNA|Histone", rownames(spe_qc), ignore.case = TRUE)
]
rowData(spe_qc)$use_channel <- rownames(spe_qc) %in% use_markers
stopifnot(length(use_markers) > 20, length(unique(spe_qc$batch)) > 1)

patient_ids <- sort(unique(as.character(spe_qc$patient_id)))
spe_qc$patient_short <- factor(
  sprintf("P%02d", match(spe_qc$patient_id, patient_ids)),
  levels = sprintf("P%02d", seq_along(patient_ids))
)

## 3. UMAP before batch correction

set.seed(220226)
spe_qc <- runUMAP(
  spe_qc,
  subset_row = use_markers,
  exprs_values = "exprs",
  ncomponents = 2,
  name = "UMAP_raw"
)

## 4. Batch correction with fastMNN

set.seed(220228)
out <- fastMNN(
  spe_qc,
  batch = spe_qc$batch,
  subset.row = use_markers,
  assay.type = "exprs",
  d = 20,
  auto.merge = TRUE,
  BSPARAM = ExactParam()
)

stopifnot(identical(colnames(spe_qc), colnames(out)))
reducedDim(spe_qc, "fastMNN") <- reducedDim(out, "corrected")

## 5. UMAP after batch correction

set.seed(220228)
spe_qc <- runUMAP(
  spe_qc,
  dimred = "fastMNN",
  ncomponents = 2,
  name = "UMAP_mnnCorrected"
)
