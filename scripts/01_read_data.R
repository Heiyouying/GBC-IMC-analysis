## 1. Packages and paths

suppressPackageStartupMessages({
  library(imcRtools)
  library(SpatialExperiment)
})

data_dir <- "data"
steinbock_dir <- file.path(data_dir, "steinbock")
read_images <- FALSE

## 2. Read and combine single-cell measurements

seg_dirs <- list.dirs(steinbock_dir, recursive = FALSE, full.names = TRUE)
seg_dirs <- sort(seg_dirs[grepl("_Segmentation$", basename(seg_dirs))])
stopifnot(length(seg_dirs) > 0)

spe_list <- lapply(seg_dirs, read_steinbock, pattern = "\\.csv$")
marker_order <- lapply(spe_list, rownames)
markers <- marker_order[[1]]

stopifnot(
  !anyDuplicated(markers),
  all(vapply(marker_order, setequal, logical(1), markers))
)

spe_list <- lapply(spe_list, function(spe) spe[markers, ])
spe_qc <- do.call(cbind, spe_list)
colnames(spe_qc) <- paste0(spe_qc$sample_id, "_", spe_qc$ObjectNumber)
stopifnot(!anyDuplicated(colnames(spe_qc)))

## 3. Add sample and patient metadata

sample_metadata <- read.csv(
  file.path(data_dir, "sample_metadata.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

metadata_columns <- c(
  "sample_id", "slide", "acq_index", "core_id", "patient_id",
  "core_label", "image_id", "batch"
)

stopifnot(
  all(metadata_columns %in% names(sample_metadata)),
  !anyDuplicated(sample_metadata$sample_id)
)

sample_metadata <- sample_metadata[, metadata_columns, drop = FALSE]
sample_index <- match(spe_qc$sample_id, sample_metadata$sample_id)
stopifnot(!anyNA(sample_index))

for (column in setdiff(metadata_columns, "sample_id")) {
  colData(spe_qc)[[column]] <- sample_metadata[[column]][sample_index]
}

stopifnot(
  !anyNA(spe_qc$batch),
  !anyNA(spe_qc$patient_id),
  all(nzchar(as.character(spe_qc$batch))),
  all(nzchar(as.character(spe_qc$patient_id)))
)

## 4. Transform intensities and select analysis channels

assay(spe_qc, "exprs") <- asinh(counts(spe_qc) / 1)
rowData(spe_qc)$use_channel <- !grepl(
  "DNA|Histone", rownames(spe_qc), ignore.case = TRUE
)

## 5. Read multichannel images and segmentation masks

if (read_images) {
  suppressPackageStartupMessages(library(cytomapper))

  image_list <- lapply(seq_along(seg_dirs), function(i) {
    current_images <- loadImages(file.path(seg_dirs[i], "img"))
    channelNames(current_images) <- marker_order[[i]]
    getChannels(current_images, markers)
  })

  mask_list <- lapply(seg_dirs, function(directory) {
    loadImages(file.path(directory, "masks"), as.is = TRUE)
  })

  images <- do.call(c, image_list)
  masks <- do.call(c, mask_list)
  sample_ids <- unique(as.character(spe_qc$sample_id))

  stopifnot(
    !anyDuplicated(names(images)),
    !anyDuplicated(names(masks)),
    all(sample_ids %in% names(images)),
    all(sample_ids %in% names(masks))
  )

  images <- getImages(images, sample_ids)
  masks <- getImages(masks, sample_ids)
  image_metadata <- sample_metadata[
    match(sample_ids, sample_metadata$sample_id), , drop = FALSE
  ]
  mcols(images) <- mcols(masks) <- S4Vectors::DataFrame(image_metadata)

  saveRDS(images, file.path(data_dir, "images.rds"), compress = "xz")
  saveRDS(masks, file.path(data_dir, "masks.rds"), compress = "xz")
}

## 6. Save the SpatialExperiment object

saveRDS(spe_qc, file.path(data_dir, "spe_qc.rds"), compress = "xz")
