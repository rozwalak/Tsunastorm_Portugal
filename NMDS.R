# Load necessary libraries
library(readxl)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(svglite)

asv_file  <- "Supplementary Table S2.xlsx"
meta_file <- "Supplementary Table S1.xlsx"

low_read_samples <- c("X19", "X58")

# Both supplementary files have a title row above the real header row,
# so the header row itself is skipped
asv_data <- read_excel(asv_file, skip = 1)

metadata <- read_excel(meta_file, skip = 1)
# Supplementary Table S1 stores the tsunami/storm/lagoon/end-member classification
# under the column "TYPE" rather than "GENETIC_TYPE" -- rename for compatibility
# with the rest of the script
metadata <- metadata %>% rename(GENETIC_TYPE = TYPE)

# Keep LAB_ID, GENETIC_TYPE, and SITE (or whatever column encodes m21/15 vs m21/27)
# Adjust column names below to match your actual sheet headers
metadata_full <- metadata[, c("LAB_ID", "GENETIC_TYPE", "FIELD_ID")]

metadata <- metadata_full[, c("LAB_ID", "GENETIC_TYPE")]

# ── Helper: run NMDS on a subset of samples ───────────────────────────────────
run_nmds <- function(asv_data, keep_cols = NULL, exclude_cols = NULL) {
  drop_cols <- c("ASV_ID", "TAXONOMY", "MEAN IDENTITY", "REFEERNCES", "TOTAL READS")
  if (!is.null(exclude_cols)) drop_cols <- c(drop_cols, exclude_cols)
  
  asv_abundance <- asv_data[, !(names(asv_data) %in% drop_cols)]
  
  # If keep_cols specified, retain only those sample columns
  if (!is.null(keep_cols)) {
    present <- intersect(keep_cols, names(asv_abundance))
    asv_abundance <- asv_abundance[, present, drop = FALSE]
  }
  
  asv_abundance <- as.data.frame(lapply(asv_abundance, as.numeric))
  asv_abundance <- asv_abundance[rowSums(asv_abundance, na.rm = TRUE) > 0, ]
  
  relative_abundance <- asv_abundance / rowSums(asv_abundance, na.rm = TRUE)
  relative_abundance[is.nan(as.matrix(relative_abundance))] <- 0
  
  asv_matrix <- t(relative_abundance)
  asv_matrix  <- asv_matrix[rowSums(asv_matrix, na.rm = TRUE) > 0, ]
  asv_matrix[is.na(asv_matrix)] <- 0
  
  set.seed(123)
  metaMDS(asv_matrix, distance = "bray", k = 2, trymax = 100,
          autotransform = FALSE)
}

build_scores <- function(nmds_result, metadata) {
  scores_df <- as.data.frame(scores(nmds_result, display = "sites"))
  scores_df$LAB_ID <- rownames(scores_df)
  merge(scores_df, metadata, by = "LAB_ID")
}

# ── Custom colour palette ─────────────────────────────────────────────────────
pal <- c(
  "tsunami"           = "#d8585c",
  "storm"             = "#32c0d2",
  "lagoon"            = "#e0b165",
  "end_member_beach"  = "#008080",
  "end_member_marine" = "#000080",
  "end_member_dune"   = "#ffeeaa",
  "end_member_stream" = "#808000",
  "end_member_soil"   = "#a05a2c"
)

# ── Identify sample sets ──────────────────────────────────────────────────────
# All sample columns (non-metadata columns in asv_data)
non_sample_cols <- c("ASV_ID", "TAXONOMY", "MEAN IDENTITY", "REFEERNCES", "TOTAL READS")
all_sample_cols <- setdiff(names(asv_data), non_sample_cols)

# Filtered = all samples minus low-read
filtered_cols <- setdiff(all_sample_cols, low_read_samples)

# Tsunami + storm samples (from metadata)
ts_samples <- metadata_full$LAB_ID[metadata_full$GENETIC_TYPE %in% c("tsunami", "storm")]
ts_cols     <- intersect(filtered_cols, ts_samples)

# m21/15 and m21/27 samples — matched by FIELD_ID prefix
field_ids     <- metadata_full$FIELD_ID
m1515_samples <- metadata_full$LAB_ID[grepl("^m21/15", field_ids, ignore.case = TRUE)]
m1527_samples <- metadata_full$LAB_ID[grepl("^m21/27", field_ids, ignore.case = TRUE)]

m1515_cols <- intersect(filtered_cols, m1515_samples)
m1527_cols <- intersect(filtered_cols, m1527_samples)

cat("Panel A (filtered):", length(filtered_cols), "samples\n")
cat("Panel B (tsunami+storm):", length(ts_cols), "samples\n")
cat("Panel C (m21/15):", length(m1515_cols), "samples\n")
cat("Panel D (m21/27):", length(m1527_cols), "samples\n")

# ── Run NMDS for each panel ───────────────────────────────────────────────────
cat("\nRunning NMDS Panel A (filtered)...\n")
nmds_A <- run_nmds(asv_data, keep_cols = filtered_cols)

cat("\nRunning NMDS Panel B (tsunami + storm)...\n")
nmds_B <- run_nmds(asv_data, keep_cols = ts_cols)

cat("\nRunning NMDS Panel C (m21/15)...\n")
nmds_C <- run_nmds(asv_data, keep_cols = m1515_cols)

cat("\nRunning NMDS Panel D (m21/27)...\n")
nmds_D <- run_nmds(asv_data, keep_cols = m1527_cols)

stress_A <- round(nmds_A$stress, 4)
stress_B <- round(nmds_B$stress, 4)
stress_C <- round(nmds_C$stress, 4)
stress_D <- round(nmds_D$stress, 4)
cat(sprintf("Stress A=%.4f  B=%.4f  C=%.4f  D=%.4f\n",
            stress_A, stress_B, stress_C, stress_D))

# ── Build score data frames ───────────────────────────────────────────────────
scores_A <- build_scores(nmds_A, metadata)
scores_B <- build_scores(nmds_B, metadata)
scores_C <- build_scores(nmds_C, metadata)
scores_D <- build_scores(nmds_D, metadata)

# ── Shared plot theme ─────────────────────────────────────────────────────────
base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    legend.position  = "none",
    panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9)
  )

make_panel <- function(scores_df, stress_val, label) {
  ggplot(scores_df, aes(x = NMDS1, y = NMDS2, fill = GENETIC_TYPE)) +
    geom_point(shape = 21, size = 3, alpha = 0.85,
               color = "black", stroke = 0.25) +
    scale_fill_manual(values = pal, name = "Genetic type",
                      guide = guide_legend(override.aes = list(size = 3.5, stroke = 0.4))) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("Stress = ", stress_val),
             hjust = 1.08, vjust = -0.6,
             size = 3.2, fontface = "italic", color = "grey35") +
    annotate("text", x = -Inf, y = Inf,
             label = label,
             hjust = -0.3, vjust = 1.5,
             size = 4, fontface = "bold", color = "black") +
    labs(x = "NMDS1", y = "NMDS2") +
    base_theme
}

p_A <- make_panel(scores_A, stress_A, "A")
p_B <- make_panel(scores_B, stress_B, "B")
p_C <- make_panel(scores_C, stress_C, "C")
p_D <- make_panel(scores_D, stress_D, "D")

# ── Extract shared legend ─────────────────────────────────────────────────────
legend_plot <- p_A +
  theme(legend.position = "right",
        legend.title    = element_text(size = 10, face = "bold"),
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.55, "cm"))

shared_legend <- cowplot::get_legend(legend_plot)

# ── Compose layout: 2 columns × 2 rows + legend strip ────────────────────────
# layout: (A | B) / (C | D) | legend
library(cowplot)

grid_2x2 <- plot_grid(
  p_A, p_B,
  p_C, p_D,
  ncol   = 2,
  nrow   = 2,
  align  = "hv",
  rel_widths  = c(1, 1),
  rel_heights = c(1, 1)
)

combined <- plot_grid(
  grid_2x2, shared_legend,
  ncol        = 2,
  rel_widths  = c(4, 0.75)
)

# ── Save outputs ──────────────────────────────────────────────────────────────
ggsave("nmds_4panel.png",
       combined, width = 13, height = 9, dpi = 300, bg = "white")

ggsave("nmds_4panel.svg",
       combined, width = 13, height = 9, device = svglite::svglite)

cat("\nSaved: nmds_4panel.png and nmds_4panel.svg\n")
