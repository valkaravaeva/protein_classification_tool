#!/usr/bin/env Rscript
# Stacked barplot of module completeness categories (Present / Partial / Absent) per taxon,
# from the matrices produced by tool_prep_matrices_for_plotting.py.
#
# Usage:
#   Rscript tool_r_stacked_barplot.R <module_matrix.tsv> <output.svg>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript tool_r_stacked_barplot.R <module_matrix.tsv> <output.svg>")
}
input_matrix <- args[1]
output_svg <- args[2]

colors <- c("#55028C", "#BB8FE4", "white")

mat <- read.csv(input_matrix, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
mat_t <- t(mat)

svg(file = output_svg, width = 300, height = 50)

barplot(mat_t, col = colors, border = "white", xlab = "taxon", font.lab = 100)
axis(2, at = 0:5, labels = 0:5)
legend("topright", colnames(mat), fill = colors, bty = "n")

dev.off()
