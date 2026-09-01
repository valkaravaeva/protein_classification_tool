#!/usr/bin/env Rscript
# Heatmap of module/marker presence-partial-absence percentages
# (e.g. the matrix_markersPRESENT.tsv / matrix_markersPARTIAL.tsv files produced by
# tool_prep_matrices_for_plotting.py).
#
# Usage:
#   Rscript tool_r_heatmap.R <matrix.tsv> <output.svg>

suppressMessages(library(pheatmap))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript tool_r_heatmap.R <matrix.tsv> <output.svg>")
}
input_matrix <- args[1]
output_svg <- args[2]

mat <- read.csv(input_matrix, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

svg(file = output_svg, width = 100, height = 160)

breaks1 <- c(0, seq(0.00000001, 101, by = 5))
# FIXED: originally read `length(breaks) - 1`, but the vector was named `breaks1` - `breaks`
# was never defined, so this would raise "object 'breaks' not found" as soon as it ran.
colors <- c("white", colorRampPalette(c("#fff7e0", "#604600"))(length(breaks1) - 1)) # yellow
# colors <- c("white", colorRampPalette(c("#e4fdff", "#015c64"))(length(breaks1) - 1)) # cyan
# colors <- c("white", colorRampPalette(c("#fbfbfb", "#000000"))(length(breaks1) - 1)) # gray
# colors <- c("white", colorRampPalette(c("#fff3f3", "#5d0000"))(length(breaks1) - 1)) # red
# colors <- c("white", colorRampPalette(c("#e5ffe9", "#003e08"))(length(breaks1) - 1)) # green

pheatmap(mat, cluster_row = FALSE, cluster_cols = FALSE, show_colnames = TRUE, show_rownames = TRUE,
         color = colors, fontsize = 70, breaks = breaks1)

dev.off()
