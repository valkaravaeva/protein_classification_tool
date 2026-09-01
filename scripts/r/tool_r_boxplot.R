#!/usr/bin/env Rscript
# Boxplot of "uncharacterized percent" per taxon, from the long-format table produced by
# tool_boxplot_prep.py (columns: Genome, Uncharacterized_percent, Taxon).
#
# Usage:
#   Rscript tool_r_boxplot.R <boxplot_table.tsv> <output.pdf>

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript tool_r_boxplot.R <boxplot_table.tsv> <output.pdf>")
}
input_table <- args[1]
output_pdf <- args[2]

dat <- read.csv(input_table, header = TRUE, sep = "\t", check.names = FALSE)

pdf(file = output_pdf, onefile = TRUE, paper = "a4", width = 300, height = 3000, compress = FALSE)

# my_colors <- rev(c("#cb181d", "#fb6a4a", "#fcae91")) ## red
my_colors <- rev(c("#2171b5", "#6baed6", "#bdd7e7")) ## blue

dat <- dat %>%
  group_by(Taxon) %>%
  mutate(med_unch = median(Uncharacterized_percent))

# cornflowerblue red
ggplot(dat, aes(x = factor(Taxon), y = Uncharacterized_percent, fill = med_unch)) +
  geom_boxplot(outlier.colour = "cornflowerblue", outlier.shape = 6, outlier.size = 2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_gradientn(colors = my_colors) +
  theme_light()

dev.off()
