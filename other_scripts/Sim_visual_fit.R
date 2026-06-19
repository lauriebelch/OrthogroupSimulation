#!/usr/bin/env Rscript

library(ape)
library(ggplot2)
library(data.table)
library(ggridges)
library(seqinr)
library(dplyr)

# =============================================================================
# Usage:
#   Rscript Sim_visual_fit.R \
#     --sim-folder      <path>   \
#     --training-folder <path>   \
#     --output          <path>.pdf \
#     [--tool-label     <label>] \
#     [--min-genes      <int>]   \
#     [--max-genes      <int>]
#
# Arguments:
#   --sim-folder      Path to the simulation output folder.
#                     Must contain Simulated_Orthogroup_statistics.txt and
#                     simulation_summaries/.
#   --training-folder Path to the OrthoTrain training output folder.
#                     Script appends TrainingResults/Orthogroup_analysis automatically.
#   --output          Output PDF file path.
#   --tool-label      Label for the training dataset in plots [default: Training].
#   --min-genes       Minimum gene count filter for simulated OGs [default: 3].
#   --max-genes       Maximum gene count filter for simulated OGs [default: 500].
# =============================================================================


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript Sim_visual_fit.R [options]\n\n")
    cat("Required:\n")
    cat("  --sim-folder      <path>   Simulation output folder\n")
    cat("  --training-folder <path>   OrthoTrain output folder\n")
    cat("  --output          <path>   Output PDF file\n\n")
    cat("Optional:\n")
    cat("  --tool-label      <label>  Label for training data [default: Training]\n")
    cat("  --min-genes       <int>    Min gene filter for simulated OGs [default: 3]\n")
    cat("  --max-genes       <int>    Max gene filter for simulated OGs [default: 500]\n")
    quit(status = 0)
  }

  result <- list(
    sim_folder      = NULL,
    training_folder = NULL,
    output          = NULL,
    tool_label      = "Training",
    min_genes       = 3L,
    max_genes       = 500L
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--sim-folder") {
      result$sim_folder <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--training-folder") {
      result$training_folder <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output") {
      result$output <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--tool-label") {
      result$tool_label <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--min-genes") {
      result$min_genes <- as.integer(args[i + 1]); i <- i + 2
    } else if (args[i] == "--max-genes") {
      result$max_genes <- as.integer(args[i + 1]); i <- i + 2
    } else {
      cat("WARNING: unknown argument:", args[i], "\n")
      i <- i + 1
    }
  }

  missing <- c()
  if (is.null(result$sim_folder))      missing <- c(missing, "--sim-folder")
  if (is.null(result$training_folder)) missing <- c(missing, "--training-folder")
  if (is.null(result$output))          missing <- c(missing, "--output")
  if (length(missing) > 0) {
    cat("ERROR: missing required arguments:", paste(missing, collapse = ", "), "\n")
    cat("Run with --help for usage.\n")
    quit(status = 1)
  }

  result
}

params          <- parse_args()
sim_folder      <- params$sim_folder
og_analysis_dir <- file.path(params$training_folder, "TrainingResults", "Orthogroup_analysis")
output_file     <- params$output
tool_label      <- params$tool_label
MIN_GENES       <- params$min_genes
MAX_GENES       <- params$max_genes

cat("Simulation folder:   ", sim_folder, "\n")
cat("Training OG dir:     ", og_analysis_dir, "\n")
cat("Output:              ", output_file, "\n")
cat("Tool label:          ", tool_label, "\n")
cat("Gene filter:         ", MIN_GENES, "to", MAX_GENES, "\n\n")

stopifnot(dir.exists(sim_folder))
stopifnot(dir.exists(og_analysis_dir))


# ---------------------------------------------------------------------------
# X-axis helpers
# ---------------------------------------------------------------------------

# For most variables: 0 to the 99th percentile of the combined data.
# This clips extreme outliers without hardcoding dataset-specific values.
xlim_auto <- function(values, q = 0.99) {
  c(0, quantile(values, q, na.rm = TRUE))
}

# Compute pretty breaks that stay within the axis limits.
breaks_auto <- function(xlim, n = 6) {
  b <- pretty(c(0, xlim[2]), n = n)
  b[b >= 0 & b <= xlim[2] * 1.05]
}


# ---------------------------------------------------------------------------
# Shared plot factory
# ---------------------------------------------------------------------------

make_ridgeplot <- function(df, x_label, xlim, x_breaks,
                           y_expand = c(0.05, 0.9)) {
  ggplot(df, aes(x = value, y = group, fill = group)) +
    ggridges::geom_density_ridges(
      alpha = 0.75,
      scale = 0.95,
      jittered_points = TRUE,
      point_size = 0.1,
      aes(point_color = group, point_fill = group)
    ) +
    labs(x = x_label, y = "") +
    scale_x_continuous(
      limits   = xlim,
      breaks   = x_breaks,
      sec.axis = dup_axis(name = "")
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 16),
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 12),
      panel.grid.minor = element_blank(),
      legend.text      = element_text(size = 14),
      legend.title     = element_blank()
    ) +
    guides(
      fill        = guide_legend(nrow = 2, reverse = TRUE),
      point_color = guide_legend(reverse = TRUE),
      point_fill  = guide_legend(reverse = TRUE)
    ) +
    scale_y_discrete(expand = expansion(add = y_expand))
}

make_df <- function(training_values, sim_values, tool_label) {
  data.frame(
    value = c(training_values, sim_values),
    group = factor(c(
      rep(tool_label,  length(training_values)),
      rep("Simulation", length(sim_values))
    ))
  )
}


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

sim_sum_all <- fread(file.path(sim_folder, "Simulated_Orthogroup_statistics.txt"))
sim_sum     <- sim_sum_all[num_genes >= MIN_GENES & num_genes < MAX_GENES]

sampled_ogs <- list.dirs(og_analysis_dir, full.names = FALSE)
sampled_ogs <- sampled_ogs[grepl("OG", sampled_ogs)]


# ---------------------------------------------------------------------------
# Open output PDF
# ---------------------------------------------------------------------------

pdf(output_file, width = 8, height = 5)
cat("Writing plots to:", output_file, "\n")


# ---------------------------------------------------------------------------
# 1. Duplications
# ---------------------------------------------------------------------------

od <- as.vector(unlist(fread(file.path(og_analysis_dir, "duplications.txt"))))
od <- od - 1
df <- make_df(od, sim_sum$num_duplications, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Number of Duplications", xl, breaks_auto(xl), c(0.01, 1)))
cat("  Duplications x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 2. Gamma shape
# ---------------------------------------------------------------------------

ofg <- fread(file.path(og_analysis_dir, "gamma_shape_values.txt"))
df  <- make_df(ofg$V1, sim_sum$gamma_shape, tool_label)
xl  <- xlim_auto(df$value)
print(make_ridgeplot(df, "Gamma Shape Value", xl, breaks_auto(xl), c(0.05, 0.75)))
cat("  Gamma shape x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 3. Proportion invariant sites   [bounded 0–1: fixed limits]
# ---------------------------------------------------------------------------

ofp <- fread(file.path(og_analysis_dir, "invariable_site_props.txt"))
df  <- make_df(ofp$V1, sim_sum$proportion_invariant, tool_label)
print(make_ridgeplot(df, "Proportion Invariant Sites", c(0, 1), seq(0, 1, 0.1), c(0.05, 1)))


# ---------------------------------------------------------------------------
# 4. Number of species
# ---------------------------------------------------------------------------

num_species <- as.vector(unlist(fread(file.path(og_analysis_dir, "num_species.txt"))))
df <- make_df(num_species, sim_sum$num_species, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Number of Species", xl, breaks_auto(xl), c(0.05, 0.9)))
cat("  Species x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 5. Number of genes
# ---------------------------------------------------------------------------

num_genes_tr <- as.vector(unlist(fread(file.path(og_analysis_dir, "num_genes.txt"))))
df <- make_df(num_genes_tr, sim_sum$num_genes, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Number of Genes", xl, breaks_auto(xl), c(0.05, 1)))
cat("  Genes x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 6. Median root-to-tip branch length
# ---------------------------------------------------------------------------

rtt    <- as.vector(unlist(fread(file.path(og_analysis_dir, "median_rtt.txt"))))
simRTT <- as.vector(unlist(fread(file.path(sim_folder, "simulation_summaries", "median_rtt.txt"))))
df <- make_df(rtt, simRTT, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Median Root-to-tip Branch Length", xl, breaks_auto(xl), c(0.05, 0.9)))
cat("  RTT x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 7. Treeness   [bounded 0–1: fixed limits]
# ---------------------------------------------------------------------------

treeness    <- as.vector(unlist(fread(file.path(og_analysis_dir, "treeness.txt"))))
treeness    <- treeness[!is.na(treeness) & treeness > 0]
simTreeness <- as.vector(unlist(fread(file.path(sim_folder, "simulation_summaries", "treeness.txt"))))
simTreeness <- simTreeness[!is.na(simTreeness)]
df <- make_df(treeness, simTreeness, tool_label)
print(make_ridgeplot(df, "Treeness", c(0, 1), seq(0, 1, 0.1), c(0.05, 0.9)))


# ---------------------------------------------------------------------------
# 8. Wiener index
# ---------------------------------------------------------------------------

wiener    <- as.vector(unlist(fread(file.path(og_analysis_dir, "wiener_index.txt"))))
simWiener <- as.vector(unlist(fread(file.path(sim_folder, "simulation_summaries", "wiener_index.txt"))))
df <- make_df(wiener, simWiener, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Wiener Index", xl, breaks_auto(xl), c(0.05, 0.9)))
cat("  Wiener x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 9. Number of gap regions per sequence
# ---------------------------------------------------------------------------

gap_num     <- as.vector(unlist(fread(file.path(og_analysis_dir, "num_gaps.txt"))))
sim_gap_num <- as.vector(unlist(fread(file.path(sim_folder, "simulation_summaries", "num_gaps.txt"))))
df <- make_df(gap_num, sim_gap_num, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Number of Gap Regions per Sequence", xl, breaks_auto(xl), c(0.05, 0.9)))
cat("  Gap count x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# 10. Median gap region size
# ---------------------------------------------------------------------------

gap_med     <- as.vector(unlist(fread(file.path(og_analysis_dir, "median_gap_size.txt"))))
sim_gap_med <- as.vector(unlist(fread(file.path(sim_folder, "simulation_summaries", "median_gap_size.txt"))))
df <- make_df(gap_med, sim_gap_med, tool_label)
xl <- xlim_auto(df$value)
print(make_ridgeplot(df, "Median Gap Region Size", xl, breaks_auto(xl), c(0.05, 0.9)))
cat("  Gap size x-max:", xl[2], "\n")


# ---------------------------------------------------------------------------
# Close PDF
# ---------------------------------------------------------------------------

dev.off()
cat("\nDone. Output written to:", output_file, "\n")
