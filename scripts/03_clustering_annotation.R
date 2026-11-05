## 1. Packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scran)
  library(bluster)
  library(scuttle)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
})

## 2. Marker expression heatmaps

plot_marker_heatmap <- function(spe, group, markers) {
  medians <- aggregateAcrossCells(
    as(spe, "SingleCellExperiment"),
    ids = colData(spe)[, group, drop = FALSE],
    statistics = "median",
    use.assay.type = "exprs",
    subset.row = markers
  )

  marker_matrix <- assay(medians, "exprs")
  colnames(marker_matrix) <- as.character(colData(medians)[[group]])
  marker_zscores <- t(scale(t(marker_matrix)))
  marker_zscores[!is.finite(marker_zscores)] <- 0

  draw(Heatmap(
    marker_zscores,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#3B4CC0", "#F7F7F7", "#B40426")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8)
  ))
}

## 3. Compare clustering parameters

compare_parameters <- function(spe, k_values, label, sample_seed = 1L) {
  clustering_matrix <- reducedDim(spe, "fastMNN")
  sweep <- clusterSweep(
    clustering_matrix,
    BLUSPARAM = SNNGraphParam(),
    k = k_values,
    type = c("rank", "jaccard"),
    cluster.fun = "louvain",
    BPPARAM = BiocParallel::SerialParam(RNGseed = 230214)
  )

  sweep_clusters <- as.data.frame(sweep$clusters)
  sweep_summary <- as.data.frame(sweep$parameters)

  set.seed(sample_seed)
  silhouette_cells <- sample(
    seq_len(nrow(clustering_matrix)),
    min(3000L, nrow(clustering_matrix))
  )
  silhouette_distances <- dist(clustering_matrix[silhouette_cells, , drop = FALSE])

  sweep_summary$avg_silhouette <- vapply(sweep_clusters, function(clusters) {
    sampled_clusters <- as.integer(factor(clusters[silhouette_cells]))
    mean(cluster::silhouette(sampled_clusters, silhouette_distances)[, "sil_width"])
  }, numeric(1))
  sweep_summary$n_clusters <- vapply(
    sweep_clusters, function(clusters) length(unique(clusters)), integer(1)
  )
  print(sweep_summary)

  print(ggplot(sweep_summary, aes(k, avg_silhouette, shape = type)) +
    geom_line(aes(group = type)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = label, x = "k", y = "Mean silhouette width", shape = "SNN weighting"))

  print(ggplot(sweep_summary, aes(k, n_clusters, shape = type)) +
    geom_line(aes(group = type)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = label, x = "k", y = "Number of clusters", shape = "SNN weighting"))

  list(clusters = sweep$clusters, parameters = sweep$parameters, summary = sweep_summary)
}


## 4. Global clustering with the selected parameters

set.seed(220702)
global_sweep <- compare_parameters(
  spe_qc, c(30L, 40L, 50L, 60L), "Global clustering"
)
spe_qc$nn_clusters <- clusterCells(
  spe_qc,
  use.dimred = "fastMNN",
  BLUSPARAM = SNNGraphParam(k = 60L, type = "rank", cluster.fun = "louvain")
)

phenotype_markers <- rownames(spe_qc)[
  !grepl("DNA|Histone", rownames(spe_qc), ignore.case = TRUE) &
    !rownames(spe_qc) %in% c("80ArAr", "190BCKG", "208Pb")
]
plot_marker_heatmap(spe_qc, "nn_clusters", phenotype_markers)

## 5. CAF subclustering with the selected parameters

caf_cluster_ids <- c()
spe_caf <- spe_qc[, spe_qc$nn_clusters %in% caf_cluster_ids]
caf_sweep <- compare_parameters(
  spe_caf, c(10L, 20L, 30L, 40L, 50L, 60L), "CAF subclustering", sample_seed = 220702
)

caf_markers <- c(
  "Vimentin", "Cadherin_11", "CD34", "CD146", "Myosin_11",
  "FSP1", "HLA_DRA", "CD74", "Fibulin", "RGS5", "MME",
  "MMP11", "CA9", "KRT8", "KRT19", "Podoplanin", "FAP",
  "Versican", "COL1A2", "Fibronectin", "Periostin", "Tenascin_C",
  "Complement_C3", "IGF1", "aSMA_115In", "aSMA_113In"
)

set.seed(220702)
spe_caf$caf_clusters <- clusterCells(
  spe_caf,
  use.dimred = "fastMNN",
  BLUSPARAM = SNNGraphParam(k = 50L, type = "rank", cluster.fun = "louvain")
)
plot_marker_heatmap(spe_caf, "caf_clusters", caf_markers)

## 6. Immune subclustering with the selected parameters

immune_cluster_ids <- c()
spe_imm <- spe_qc[, spe_qc$nn_clusters %in% immune_cluster_ids]
immune_sweep <- compare_parameters(
  spe_imm, c(10L, 20L, 30L, 40L, 50L, 60L), "Immune subclustering"
)

immune_markers <- c("CD45", "CD3", "CD8a", "CD68", "HLA_DRA", "CD74", "Ki67")

set.seed(220712)
spe_imm$imm_clusters <- clusterCells(
  spe_imm,
  use.dimred = "fastMNN",
  BLUSPARAM = SNNGraphParam(k = 40L, type = "rank", cluster.fun = "louvain")
)
plot_marker_heatmap(spe_imm, "imm_clusters", immune_markers)

## Final curated cell annotations are provided in data/cell_annotations.csv and data/spe_annotated.rds.
