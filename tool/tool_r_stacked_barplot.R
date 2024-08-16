colors <- c("#55028C","#BB8FE4","white")

file<-paste0("/path/to/your_module_matrix.tsv")
matrix <- read.csv(file, header = TRUE, sep ="\t", row.names = 1, check.names = FALSE)

df_t <- t(matrix)

svgfile <- paste0("/path/to/save/stacked_barplot.svg")
svg(file = svgfile, width = 300, height=50)

barplot(df_t, col=colors, border="white", xlab="taxon",font.lab=100)

axis(2, at = 0:5, labels = 0:5)
legend("topright", colnames(matrix), fill = colors, bty = "n")

dev.off()
