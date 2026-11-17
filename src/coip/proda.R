#!/usr/bin/env Rscript
# Differential abundance with proDA, for one contrast.
#
# Called by src/coip/differential.py, not by hand:
#
#   Rscript proda.R <intensities.csv> <design.csv> <treatment> <reference> <out.csv>
#
# intensities.csv  protein_id in the first column, one column per sample,
#                  log2 intensities, empty cells where a protein was not
#                  measured. Unlike limma.R, missing values are the point:
#                  proDA models them with a probabilistic dropout curve
#                  instead of dropping the protein.
#
# Nothing is imputed. A missing value stays NA and contributes through the
# dropout model, which is what distinguishes this from filling in a small
# number and calling it a measurement.
#
# proDA does not normalise between samples. median_normalization() is a
# separate function that the caller applies or does not; moderate_location is
# a prior on each protein's own mean, not a loading correction. So whatever
# scale the input is on is the scale the result is on.
#
# Reference: Ahlmann-Eltze and Anders 2020, bioRxiv 661496; proDA package
# vignette "Introduction to proDA".

suppressPackageStartupMessages(library(proDA))

# proDA's hyper-parameter fit is stochastic. Unseeded, two runs of the same
# data differ by about 0.01 in log2 fold change , far too little to change a
# conclusion, and far too much for the output tables to be comparable between
# runs. The seed makes the result reproducible; it is not a tuning knob.
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: proda.R <intensities.csv> <design.csv> <treatment> <reference> <out.csv>")
}
intensity_file <- args[1]; design_file <- args[2]
treatment <- args[3]; reference <- args[4]; out_file <- args[5]

expr <- as.matrix(read.csv(intensity_file, check.names = FALSE, row.names = 1))
design_table <- read.csv(design_file, check.names = FALSE, colClasses = "character")
if (!identical(colnames(expr), design_table$sample)) {
  stop("column order differs from the design")
}
group <- factor(design_table$group)
for (level in c(treatment, reference)) {
  if (!level %in% levels(group)) stop(sprintf("group '%s' is not in the design", level))
}
if (nrow(expr) == 0L) stop("no proteins to test")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- proDA(expr, design = design, data_is_log_transformed = TRUE, verbose = FALSE)
result <- test_diff(fit, contrast = sprintf("%s - %s", treatment, reference),
                    pval_adjust_method = "BH", sort_by = NULL)

out <- data.frame(
  protein_id = result$name,
  log2FC     = result$diff,
  t          = result$t_statistic,
  se         = result$se,
  df         = result$df,
  p_value    = result$pval,
  p_adjusted = result$adj_pval,
  n_obs      = result$n_obs,
  n_approx   = result$n_approx,
  check.names = FALSE, stringsAsFactors = FALSE)
write.csv(out, out_file, row.names = FALSE)

parameters <- hyper_parameters(fit)
cat(sprintf(paste("proDA %s: %d proteins, contrast %s - %s, %d with a p-value;",
                  "location prior mean %.2f scale %.3f df %.1f;",
                  "dropout curve position %s scale %s\n"),
            as.character(packageVersion("proDA")), nrow(out), treatment, reference,
            sum(!is.na(out$p_value)),
            parameters$location_prior_mean, parameters$location_prior_scale,
            parameters$location_prior_df,
            paste(round(parameters$dropout_curve_position, 2), collapse = "/"),
            paste(round(parameters$dropout_curve_scale, 2), collapse = "/")))
