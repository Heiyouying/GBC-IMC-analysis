## 1. Packages

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
})

## 2. Read clinical data

spe_qc <- readRDS("data/spe_annotated.rds")
clinical <- read.csv(
  "data/clinical.csv",
  fileEncoding = "UTF-8-BOM",
  na.strings = c("", "NA", "Rx", "Nx", "Tx", "Mx"),
  check.names = FALSE
)

caf_types <- c("tCAF", "vCAF", "matCAF", "iCAF", "apCAF")
immune_types <- c(
  "Antigen-presenting cells", "CD8+ T cells", "Proliferative T cells",
  "Other T cells", "Classical macrophages", "Non-classical macrophages"
)
clinical_vars <- c(
  "Gender", "Age", "T_stage", "N_stage", "M_stage", "Grade",
  "Residual_Tumor", "Lymphatic_Invasion", "Venous_Invasion",
  "Perineural_Invasion", "UICC"
)

## 3. Calculate cell proportions within each core

proportion_table <- function(cell_counts) {
  proportions <- prop.table(cell_counts, margin = 1)
  proportions[is.nan(proportions)] <- NA_real_
  data.frame(
    core_id = rownames(proportions), unclass(proportions),
    row.names = NULL, check.names = FALSE
  )
}

major_counts <- table(spe_qc$core_id, spe_qc[["major cell type"]])
subtype_counts <- table(spe_qc$core_id, spe_qc[["all cell type"]])
major_proportions <- proportion_table(major_counts)
caf_proportions <- proportion_table(subtype_counts[, caf_types, drop = FALSE])
immune_proportions <- proportion_table(subtype_counts[, immune_types, drop = FALSE])

clinical_combined <- clinical %>%
  left_join(major_proportions, by = "core_id") %>%
  left_join(caf_proportions, by = "core_id") %>%
  left_join(immune_proportions, by = "core_id")

caf_long <- clinical_combined %>%
  pivot_longer(all_of(caf_types), names_to = "CAF_type", values_to = "Proportion") %>%
  mutate(CAF_type = factor(CAF_type, levels = caf_types))

## 4. Box plots of CAF proportions

box_data <- caf_long %>%
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

  print(ggplot(plot_data, aes(CAF_type, Proportion, fill = Group)) +
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
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "top") +
    labs(x = NULL, y = "CAF subtype proportion", fill = variable))
}

## 5. Kaplan-Meier curves

survival_vars <- c("CAFs", caf_types)
for (cell_type in survival_vars) {
  km_data <- clinical_combined %>%
    transmute(time, event, Proportion = .data[[cell_type]]) %>%
    drop_na()
  cutpoint <- surv_cutpoint(
    km_data, time = "time", event = "event",
    variables = "Proportion", minprop = 0.1, progressbar = FALSE
  )
  km_data$Group <- factor(
    surv_categorize(cutpoint, labels = c("Low", "High"))$Proportion,
    levels = c("High", "Low")
  )
  fit <- survfit(Surv(time, event) ~ Group, data = km_data)

  print(ggsurvplot(
    fit, data = km_data, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    palette = c("#D73027", "#2B83BA"),
    legend.title = cell_type, legend.labs = c("High", "Low"),
    xlim = c(0, 5), break.time.by = 1,
    surv.scale = "percent", xlab = "Time (years)", ylab = "Overall survival (%)"
  ))
}

## 6. Ro/e analysis

for (variable in clinical_vars) {
  roe_data <- caf_long %>%
    filter(!is.na(.data[[variable]]), !is.na(Proportion)) %>%
    mutate(Group = .data[[variable]]) %>%
    group_by(CAF_type) %>%
    mutate(Expected = mean(Proportion)) %>%
    group_by(Group, CAF_type) %>%
    summarise(
      Observed = mean(Proportion), Expected = first(Expected), .groups = "drop"
    ) %>%
    mutate(
      Ro_e = Observed / Expected,
      Group = factor(Group, levels = rev(sort(unique(Group)))),
      Enrichment = cut(
        Ro_e, breaks = c(-Inf, 0.8, 1.3, 1.8, Inf),
        labels = c("+/-", "+", "++", "+++")
      )
    )

  print(ggplot(roe_data, aes(CAF_type, Group, fill = Ro_e)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Enrichment, color = Enrichment), fontface = "bold") +
    scale_fill_gradientn(
      colors = c("#FDF0E6", "#FD8D3C", "#A63603"),
      limits = c(0, 2), oob = scales::squish, name = "Ro/e"
    ) +
    scale_color_manual(
      values = rep("black", 4), drop = FALSE, name = NULL,
      labels = c("+/-: 0-0.8", "+: 0.8-1.3", "++: 1.3-1.8", "+++: >1.8")
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "top") +
    labs(x = NULL, y = variable))
}
