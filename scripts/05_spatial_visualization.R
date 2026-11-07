## 1. Packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(ggplot2)
})

## 2. Reading

spe_qc <- readRDS("data/spe_annotated.rds")
coordinates <- spatialCoords(spe_qc)
plot_data <- data.frame(
  core_id = as.character(spe_qc$core_id),
  x = coordinates[, 1],
  y = coordinates[, 2],
  "major cell type" = as.character(spe_qc[["major cell type"]]),
  "all cell type" = as.character(spe_qc[["all cell type"]]),
  check.names = FALSE
)

## 3. Cell type colors

major_colours <- c(
  "Tumor cells" = "#B40426",
  "CAFs" = "#3B7A57",
  "Immune cells" = "#4A7BB0",
  "Epithelial" = "#FFE4B5",
  "Endothelial" = "#E28743",
  "Undefined cells" = "#8C969E"
)

all_colours <- c(
  "Tumor cells" = "#B40426",
  "Intratumoral T cells" = "#F59999",
  "Intratumoral macrophages" = "#D4A345",
  "CD8+ T cells" = "#40A6CC",
  "Proliferative T cells" = "#8CC7E0",
  "Other T cells" = "#1F388C",
  "Classical macrophages" = "#C78C26",
  "Non-classical macrophages" = "#7A4F26",
  "Antigen-presenting cells" = "#EBB847",
  "apCAF" = "#B8E6C2",
  "tCAF" = "#1F4729",
  "vCAF" = "#26808C",
  "iCAF" = "#7AB270",
  "matCAF" = "#4F7038",
  "Epithelial" = "#FFE4B5",
  "Endothelial" = "#E28743",
  "Undefined cells" = "#8C969E"
)

## 4. Spatial plots for individual cores

plot_spatial <- function(data, annotation, colours) {
  ggplot(data, aes(x, y, colour = .data[[annotation]])) +
    geom_point(size = 1, alpha = 0.85) +
    scale_colour_manual(values = colours, breaks = names(colours)) +
    coord_equal() +
    scale_y_reverse() +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = unique(data$core_id), subtitle = annotation, colour = NULL)
}

core_ids <- unique(plot_data$core_id)

for (core_id in core_ids) {
  core_data <- plot_data[plot_data$core_id == core_id, ]
  print(plot_spatial(core_data, "major cell type", major_colours))
  print(plot_spatial(core_data, "all cell type", all_colours))
}
