## 1. Packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## 2. Reading

spe_qc <- readRDS("data/spe_annotated.rds")
tumour_fraction <- tapply(
  spe_qc[["major cell type"]] == "Tumor cells", spe_qc$core_id, mean
)
spe_qc <- spe_qc[, spe_qc$core_id %in% names(tumour_fraction)[tumour_fraction >= 0.03]]

known_cells <- data.frame(
  core_id = as.character(spe_qc$core_id),
  Cell_type = as.character(spe_qc[["all cell type"]]),
  x = spatialCoords(spe_qc)[, 1], y = spatialCoords(spe_qc)[, 2]
) %>% filter(Cell_type != "Undefined cells")

cell_colours <- c(
  "Tumor cells" = "#E03826", "Intratumoral T cells" = "#FAA3B5",
  "Intratumoral macrophages" = "#D4A345", "apCAF" = "#B8E6C2",
  "iCAF" = "#7AB270", "matCAF" = "#4F7038", "tCAF" = "#1E4729",
  "vCAF" = "#26808C", "Antigen-presenting cells" = "#EBB847",
  "CD8+ T cells" = "#3FADC9", "Proliferative T cells" = "#8CC6E0",
  "Other T cells" = "#1F388C", "Classical macrophages" = "#C78C26",
  "Non-classical macrophages" = "#7A4F26", "Epithelial" = "#FFE4B5",
  "Endothelial" = "#E28743"
)
known_types <- names(cell_colours)

## 3. Cell-centered window composition within 20 micrometers

window_radius <- 20
min_neighbours <- 3
window_results <- list()

for (core in unique(known_cells$core_id)) {
  core_cells <- filter(known_cells, core_id == core)
  neighbours <- dbscan::frNN(as.matrix(core_cells[, c("x", "y")]), eps = window_radius)$id
  neighbours <- lapply(seq_along(neighbours), function(i) setdiff(neighbours[[i]], i))
  eligible <- which(lengths(neighbours) >= min_neighbours)
  if (length(eligible) == 0) next
  local_composition <- t(vapply(neighbours[eligible], function(index) {
    tabulate(match(core_cells$Cell_type[index], known_types), nbins = length(known_types)) / length(index)
  }, numeric(length(known_types))))
  colnames(local_composition) <- known_types
  window_results[[core]] <- data.frame(
    Center_type = core_cells$Cell_type[eligible], local_composition, check.names = FALSE
  ) %>%
    group_by(Center_type) %>%
    summarise(across(all_of(known_types), mean), .groups = "drop") %>%
    pivot_longer(all_of(known_types), names_to = "Neighbor_type", values_to = "Proportion") %>%
    mutate(core_id = core)
}
window_by_core <- bind_rows(window_results)
window_composition <- window_by_core %>%
  group_by(Center_type, Neighbor_type) %>%
  summarise(Proportion = mean(Proportion), .groups = "drop")

for (center in known_types) {
  plot_data <- window_composition %>% filter(Center_type == center) %>%
    arrange(desc(Proportion)) %>% mutate(Neighbor_type = factor(Neighbor_type, levels = Neighbor_type))
  if (nrow(plot_data) == 0) next
  print(ggplot(plot_data, aes(Neighbor_type, Proportion)) +
    geom_col(fill = "#3B5BA5", colour = "#254375", width = 0.6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Relative frequency", title = paste(center, "20μm window")))
}

## 4. Nearest-neighbour distances between cell types

distance_limit <- 75
distance_results <- list()
for (core in unique(known_cells$core_id)) {
  core_cells <- filter(known_cells, core_id == core)
  coordinates <- as.matrix(core_cells[, c("x", "y")])
  for (target in unique(core_cells$Cell_type)) {
    target_index <- which(core_cells$Cell_type == target)
    target_coordinates <- coordinates[target_index, , drop = FALSE]
    distances <- FNN::get.knnx(target_coordinates, coordinates, k = 1)$nn.dist[, 1]
    distances[target_index] <- if (length(target_index) > 1)
      FNN::get.knn(target_coordinates, k = 1)$nn.dist[, 1] else NA_real_
    distance_results[[paste(core, target)]] <- data.frame(
      core_id = core, Center_type = core_cells$Cell_type,
      Target_type = target, Distance = distances
    ) %>%
      filter(!is.na(Distance), Distance <= distance_limit) %>%
      group_by(core_id, Center_type, Target_type) %>%
      summarise(Mean_distance = mean(Distance), .groups = "drop")
  }
}
distance_by_core <- bind_rows(distance_results)

for (center in known_types) {
  plot_data <- distance_by_core %>% filter(Center_type == center)
  target_order <- plot_data %>% group_by(Target_type) %>%
    summarise(Median = median(Mean_distance), .groups = "drop") %>%
    arrange(Median) %>% pull(Target_type)
  plot_data$Target_type <- factor(plot_data$Target_type, levels = target_order)
  print(ggplot(plot_data, aes(Target_type, Mean_distance, fill = Target_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.5) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
    scale_fill_manual(values = cell_colours) +
    theme_bw() +
    theme(
      panel.grid = element_blank(), legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = NULL, y = "Mean nearest-neighbour distance (75μm)", title = paste("Distance from", center)))
}
