## 1. Packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(imcRtools)
  library(BiocParallel)
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
spe_qc[["all cell type"]] <- droplevels(factor(spe_qc[["all cell type"]]))
cell_types <- levels(spe_qc[["all cell type"]])

## 3. Build the 15-micrometre spatial graph

spe_qc <- buildSpatialGraph(
  spe_qc, img_id = "core_id", type = "expansion", threshold = 15,
  name = "expansion_graph_15um"
)

## 4. Cell-cell interaction and avoidance testing

set.seed(221029)
interaction_results <- as.data.frame(testInteractions(
  spe_qc, group_by = "core_id", label = "all cell type",
  colPairName = "expansion_graph_15um", method = "classic",
  iter = 1000, p_threshold = 0.01,
  BPPARAM = SerialParam(RNGseed = 221029)
))

## 5. Summarise cell type pairs across cores

interaction_summary <- interaction_results %>%
  filter(!is.na(sigval)) %>%
  group_by(from_label, to_label) %>%
  summarise(sum_sigval = sum(sigval), .groups = "drop") %>%
  complete(
    from_label = cell_types, to_label = cell_types,
    fill = list(sum_sigval = 0)
  )

## 6. Interaction and avoidance heatmap

heatmap_data <- interaction_summary %>%
  mutate(
    from_label = factor(from_label, levels = cell_types),
    to_label = factor(to_label, levels = rev(cell_types))
  )

print(ggplot(heatmap_data, aes(from_label, to_label, fill = sum_sigval)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
    midpoint = 0
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9), panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, fill = "Sum sigval", title = "Cell-cell interaction and avoidance"))
