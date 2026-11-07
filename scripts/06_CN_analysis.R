## 1. Packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(imcRtools)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(ggalluvial)
  library(ggrepel)
  library(fmsb)
  library(survival)
  library(survminer)
})

## 2. Reading

spe_qc <- readRDS("data/spe_annotated.rds")
clinical <- read.csv(
  "data/clinical.csv", fileEncoding = "UTF-8-BOM",
  na.strings = c("", "NA", "Rx", "Nx", "Tx", "Mx"), check.names = FALSE
)

tumour_fraction <- tapply(
  spe_qc[["major cell type"]] == "Tumor cells", spe_qc$core_id, mean
)
spe_qc <- spe_qc[, spe_qc$core_id %in% names(tumour_fraction)[tumour_fraction >= 0.03]]
spe_qc[["all cell type"]] <- as.character(spe_qc[["all cell type"]])

cell_types <- c(
  "Tumor cells", "Intratumoral T cells", "Intratumoral macrophages",
  "apCAF", "iCAF", "matCAF", "tCAF", "vCAF", "Antigen-presenting cells",
  "CD8+ T cells", "Proliferative T cells", "Other T cells",
  "Classical macrophages", "Non-classical macrophages",
  "Epithelial", "Endothelial", "Undefined cells"
)
cell_colours <- c(
  "Tumor cells" = "#E03826", "Intratumoral T cells" = "#FAA3B5",
  "Intratumoral macrophages" = "#D4A345", "apCAF" = "#B8E6C2",
  "iCAF" = "#7AB270", "matCAF" = "#4F7038", "tCAF" = "#1E4729",
  "vCAF" = "#26808C", "Antigen-presenting cells" = "#EBB847",
  "CD8+ T cells" = "#3FADC9", "Proliferative T cells" = "#8CC6E0",
  "Other T cells" = "#1F388C", "Classical macrophages" = "#C78C26",
  "Non-classical macrophages" = "#7A4F26", "Epithelial" = "#FFE4B5",
  "Endothelial" = "#E28743", "Undefined cells" = "#8C969E"
)

## 3. Identify CNs using 20 nearest neighbours within 75 micrometers

spe_qc <- buildSpatialGraph(
  spe_qc, img_id = "core_id", type = "knn", k = 20, max_dist = 75,
  name = "knn20_75um"
)
spe_qc <- aggregateNeighbors(
  spe_qc, colPairName = "knn20_75um", aggregate_by = "metadata",
  count_by = "all cell type", proportions = TRUE, name = "neighbour_composition"
)

neighbour_composition <- as.matrix(spe_qc$neighbour_composition)
included <- complete.cases(neighbour_composition)
set.seed(123)
cn_clustering <- kmeans(neighbour_composition[included, , drop = FALSE], centers = 7, nstart = 50)
cn_levels <- paste0("CN", 1:7)
cn_labels <- rep(NA_character_, ncol(spe_qc))
cn_labels[included] <- paste0("CN", cn_clustering$cluster)
spe_qc$CN <- factor(cn_labels, levels = cn_levels)

cn_counts <- table(spe_qc$CN, factor(spe_qc[["all cell type"]], levels = cell_types))
cn_composition <- unclass(prop.table(cn_counts, margin = 1))

## 4. CN composition heatmap with column Z-scores

heatmap_data <- scale(cn_composition)
heatmap_data[!is.finite(heatmap_data)] <- 0
draw(Heatmap(
  heatmap_data, name = "Z-score",
  col = circlize::colorRamp2(c(-3, 0, 3), c("#3B5BA5", "#F7F7F7", "#B40426")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  rect_gp = grid::gpar(col = "white"), border = TRUE,
  column_names_rot = 45, column_names_gp = grid::gpar(fontsize = 8)
))

## 5. CN and cell type Sankey plot
cn_colours <- setNames(
  c("#EDC948", "#59A14F", "#76B7B2", "#E15759", "#9C755F", "#F28E2B", "#4E79A7"),
  cn_levels
)
sankey_data <- as.data.frame(cn_counts, responseName = "Frequency")
names(sankey_data)[1:2] <- c("CN", "Cell_type")
sankey_data <- sankey_data %>% filter(Frequency > 5) %>% mutate(CN_label = CN)
sankey_long <- to_lodes_form(
  sankey_data, axes = c("CN_label", "Cell_type"),
  key = "Axis", value = "Node", id = "Flow"
)

print(ggplot(sankey_long, aes(Axis, Frequency, stratum = Node, alluvium = Flow)) +
  geom_alluvium(aes(fill = CN), width = 0.08, alpha = 0.6) +
  geom_stratum(
    data = filter(sankey_long, Axis == "CN_label"),
    aes(fill = CN), width = 0.04, colour = "grey30"
  ) +
  geom_stratum(
    data = filter(sankey_long, Axis == "Cell_type"),
    fill = "white", width = 0.04, colour = "grey30"
  ) +
  geom_text(
    data = filter(sankey_long, Axis == "CN_label"), stat = "stratum",
    aes(label = after_stat(stratum)), nudge_x = -0.06, hjust = 1, size = 3
  ) +
  geom_text_repel(
    data = filter(sankey_long, Axis == "Cell_type"), stat = "stratum",
    aes(label = after_stat(stratum)), direction = "y", nudge_x = 0.15,
    hjust = 0, size = 3, min.segment.length = 0, max.overlaps = Inf, seed = 123
  ) +
  scale_fill_manual(values = cn_colours) +
  scale_x_discrete(labels = c("Cellular neighbourhoods", "Cell types"), expand = c(0.4, 0.4)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = NULL, y = "Cell count"))

## 6. CN marker-expression radar plots

radar_markers <- c(
  "CD34", "CD146", "Podoplanin", "Periostin", "Versican", "Myosin_11",
  "CD3", "CD8a", "CD68", "CD74", "HLA_DRA", "Complement_C3",
  "KRT8", "KRT19", "CA9", "Ki67"
)
radar_means <- data.frame(
  CN = spe_qc$CN, t(assay(spe_qc, "exprs")[radar_markers, ]), check.names = FALSE
) %>%
  filter(!is.na(CN)) %>%
  group_by(CN) %>%
  summarise(across(all_of(radar_markers), mean), .groups = "drop")
radar_limit <- max(2, ceiling(max(as.matrix(radar_means[, radar_markers])) * 2) / 2)

previous_par <- par(mar = c(2, 2, 4, 2))
for (cn in cn_levels) {
  radar_data <- as.data.frame(rbind(
    rep(radar_limit, length(radar_markers)), rep(0, length(radar_markers)),
    radar_means[radar_means$CN == cn, radar_markers]
  ))
  names(radar_data) <- radar_markers
  radarchart(
    radar_data, axistype = 1, seg = 4,
    caxislabels = round(seq(0, radar_limit, length.out = 5), 2),
    pcol = cn_colours[cn], pfcol = adjustcolor(cn_colours[cn], alpha.f = 0.3),
    pty = 32, plty = 1, plwd = 2, cglty = 1, cglcol = "grey85",
    calcex = 0.75, vlcex = 0.8, title = cn
  )
}
par(previous_par)

## 7. CN proportions and clinical relevance

core_cn_counts <- table(as.character(spe_qc$core_id), spe_qc$CN)
core_cn_proportions <- data.frame(
  core_id = rownames(core_cn_counts), unclass(prop.table(core_cn_counts, margin = 1)),
  row.names = NULL, check.names = FALSE
)
clinical_cn <- inner_join(clinical, core_cn_proportions, by = "core_id")
cn_long <- clinical_cn %>%
  pivot_longer(all_of(cn_levels), names_to = "CN", values_to = "Proportion") %>%
  mutate(CN = factor(CN, levels = cn_levels))
clinical_vars <- c(
  "Gender", "Age", "T_stage", "N_stage", "M_stage", "Grade",
  "Residual_Tumor", "Lymphatic_Invasion", "Venous_Invasion",
  "Perineural_Invasion", "UICC"
)

box_data <- cn_long %>%
  mutate(
    T_stage = recode(T_stage, T1 = "T1-T2", T2 = "T1-T2", T3 = "T3-T4", T4 = "T3-T4"),
    Grade = recode(Grade, G1 = "G1-G2", G2 = "G1-G2"),
    Residual_Tumor = recode(Residual_Tumor, R1 = "R1-R2", R2 = "R1-R2"),
    UICC = recode(UICC, UICC3 = "UICC3-UICC4", UICC4 = "UICC3-UICC4")
  )

for (variable in clinical_vars) {
  plot_data <- box_data %>%
    filter(!is.na(.data[[variable]]), !is.na(Proportion)) %>%
    mutate(Group = factor(.data[[variable]]))
  print(ggplot(plot_data, aes(CN, Proportion, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, position = position_dodge(0.7)) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
      size = 0.6, alpha = 0.7, show.legend = FALSE
    ) +
    stat_compare_means(
      aes(group = Group), method = "wilcox.test", method.args = list(exact = FALSE),
      label = "p.signif", hide.ns = TRUE
    ) +
    scale_fill_manual(values = c("#2159A1", "#E15759")) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "top") +
    labs(x = NULL, y = "CN proportion", fill = variable))
}

for (cn in cn_levels) {
  km_data <- clinical_cn %>% transmute(time, event, Proportion = .data[[cn]]) %>% drop_na()
  cutpoint <- surv_cutpoint(
    km_data, time = "time", event = "event", variables = "Proportion",
    minprop = 0.1, progressbar = FALSE
  )
  km_data$Group <- factor(
    surv_categorize(cutpoint, labels = c("Low", "High"))$Proportion,
    levels = c("High", "Low")
  )
  fit <- survfit(Surv(time, event) ~ Group, data = km_data)
  print(ggsurvplot(
    fit, data = km_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    palette = c("#D73027", "#2B83BA"), legend.title = cn,
    legend.labs = c("High", "Low"), xlim = c(0, 5), break.time.by = 1,
    surv.scale = "percent", xlab = "Time (years)", ylab = "Overall survival (%)"
  ))
}

for (variable in clinical_vars) {
  roe_data <- cn_long %>%
    filter(!is.na(.data[[variable]]), !is.na(Proportion)) %>%
    mutate(Group = .data[[variable]]) %>%
    group_by(CN) %>% mutate(Expected = mean(Proportion)) %>%
    group_by(Group, CN) %>%
    summarise(Observed = mean(Proportion), Expected = first(Expected), .groups = "drop") %>%
    mutate(
      Ro_e = Observed / Expected,
      Group = factor(Group, levels = rev(sort(unique(Group)))),
      Enrichment = cut(
        Ro_e, breaks = c(-Inf, 0.8, 1.3, 1.8, Inf), labels = c("+/-", "+", "++", "+++")
      )
    )
  print(ggplot(roe_data, aes(CN, Group, fill = Ro_e)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = Enrichment, colour = Enrichment), fontface = "bold") +
    scale_fill_gradientn(
      colours = c("#FDF0E6", "#FD8D3C", "#A63603"),
      limits = c(0, 2.5), oob = scales::squish, name = "Ro/e"
    ) +
    scale_colour_manual(
      values = rep("black", 4), drop = FALSE, name = NULL,
      labels = c("+/-: 0-0.8", "+: 0.8-1.3", "++: 1.3-1.8", "+++: >1.8")
    ) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "top") +
    labs(x = NULL, y = variable))
}

## 8. Cellular neighborhood interface protein expression gradient analysis

cell_meta <- data.frame(
  cell_id = colnames(spe_qc), core_id = as.character(spe_qc$core_id),
  Cell_type = spe_qc[["all cell type"]], CN = as.character(spe_qc$CN),
  x = spatialCoords(spe_qc)[, 1], y = spatialCoords(spe_qc)[, 2]
)
tumour_cn <- rownames(cn_composition)[which.max(cn_composition[, "Tumor cells"])]
target_cns <- setdiff(cn_levels, tumour_cn)
interface_range <- 50
interface_results <- list()

for (core in unique(cell_meta$core_id)) {
  core_cells <- filter(cell_meta, core_id == core)
  tumour_cells <- filter(core_cells, CN == tumour_cn)
  if (nrow(tumour_cells) == 0) next
  for (cn in target_cns) {
    target_cells <- filter(core_cells, CN == cn)
    if (nrow(target_cells) == 0) next
    region <- filter(core_cells, CN %in% c(tumour_cn, cn))
    region$Distance <- ifelse(
      region$CN == tumour_cn,
      -FNN::get.knnx(target_cells[, c("x", "y")], region[, c("x", "y")], k = 1)$nn.dist[, 1],
      FNN::get.knnx(tumour_cells[, c("x", "y")], region[, c("x", "y")], k = 1)$nn.dist[, 1]
    )
    interface_results[[paste(core, cn)]] <- region %>%
      filter(abs(Distance) <= interface_range) %>% mutate(Interface = cn)
  }
}
interface_data <- bind_rows(interface_results)
gradient_markers <- setdiff(
  rownames(spe_qc), c("80ArAr", "190BCKG", "208Pb", "DNA1", "DNA2", "Histone_H3")
)

for (cell_group in c("All cells", cell_types)) {
  group_data <- if (cell_group == "All cells") interface_data else
    filter(interface_data, Cell_type == cell_group)
  group_data <- group_data %>%
    group_by(Interface) %>% filter(n_distinct(Distance) >= 4) %>% ungroup()
  if (nrow(group_data) == 0) next
  k_basis <- min(10L, min(tapply(group_data$Distance, group_data$Interface, function(x) length(unique(x)))))
  for (marker in gradient_markers) {
    plot_data <- group_data %>%
      mutate(Expression = as.numeric(assay(spe_qc, "exprs")[marker, cell_id]))
    print(ggplot(plot_data, aes(Distance, Expression, colour = Interface)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_smooth(
        method = "gam", formula = y ~ s(x, bs = "cs", k = k_basis),
        se = TRUE, linewidth = 1.2, alpha = 0.08
      ) +
      scale_colour_manual(values = cn_colours, breaks = target_cns) +
      theme_bw() +
      labs(
        x = paste0("Inside ", tumour_cn, " | Outside ", tumour_cn, " (50μm)"),
        y = paste(marker, "expression"), title = paste(cell_group, marker, sep = ": "),
        colour = paste("Interface with", tumour_cn)
      ))
  }
}