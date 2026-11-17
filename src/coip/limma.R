#!/usr/bin/env Rscript
# Moderated t-test for one contrast, using limma.
#
# Called by src/coip/differential.py, not by hand:
#
#   Rscript limma.R <intensities.csv> <design.csv> <treatment> <reference> <out.csv>
#
# intensities.csv  protein_id in the first column, one column per sample,
#                  log2 intensities, no missing values
# design.csv       columns 'sample' and 'group'; the sample order must match
#                  the column order of the intensity file
#
# The model is ~ 0 + group, so each group gets its own coefficient and the
# contrast treatment - reference is the log2 fold change. eBayes shares
# variance information across proteins, which is what makes a design with two
# replicates per group estimable at all. P-values are adjusted with
# Benjamini-Hochberg over the proteins actually tested.
#
# Reference: Ritchie et al. 2015, Nucleic Acids Research 43:e47; limma User's
# Guide, section 9.2 (two groups) and 13.1 (contrasts).

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: limma.R <intensities.csv> <design.csv> <treatment> <reference> <out.csv>")
}
intensity_file <- args[1]
design_file    <- args[2]
treatment      <- args[3]
reference      <- args[4]
out_file       <- args[5]

expr <- read.csv(intensity_file, check.names = FALSE, row.names = 1)
design_table <- read.csv(design_file, check.names = FALSE, colClasses = "character")

if (!all(c("sample", "group") %in% colnames(design_table))) {
  stop("the design file needs a 'sample' and a 'group' column")
}
if (!identical(colnames(expr), design_table$sample)) {
  stop(sprintf("column order differs from the design: [%s] vs [%s]",
               paste(colnames(expr), collapse = ", "),
               paste(design_table$sample, collapse = ", ")))
}
if (anyNA(expr)) {
  stop("limma was given missing values; filter to complete cases first")
}
if (nrow(expr) == 0L) stop("no proteins to test")

group <- factor(design_table$group)
for (level in c(treatment, reference)) {
  if (!level %in% levels(group)) {
    stop(sprintf("group '%s' is not in the design", level))
  }
}

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Residual degrees of freedom before any moderation. With two samples per
# group and no further terms this is 2. It is written out so the notebook can
# state what the test actually rests on.
residual_df <- nrow(design) - qr(design)$rank
if (residual_df < 1L) {
  stop("no residual degrees of freedom: the design is saturated")
}

fit <- lmFit(as.matrix(expr), design)
contrast_matrix <- makeContrasts(
  contrasts = sprintf("%s - %s", treatment, reference), levels = design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

result <- topTable(fit, coef = 1, number = Inf, adjust.method = "BH",
                   sort.by = "none")
out <- data.frame(
  protein_id   = rownames(result),
  log2FC       = result$logFC,
  avg_log2     = result$AveExpr,
  t            = result$t,
  p_value      = result$P.Value,
  p_adjusted   = result$adj.P.Val,
  B            = result$B,
  residual_df  = residual_df,
  moderated_df = residual_df + fit$df.prior[1],
  check.names  = FALSE, stringsAsFactors = FALSE)
write.csv(out, out_file, row.names = FALSE)

cat(sprintf("limma %s: %d proteins, contrast %s - %s, residual df %d, prior df %.3f\n",
            as.character(packageVersion("limma")), nrow(out),
            treatment, reference, residual_df, fit$df.prior[1]))
