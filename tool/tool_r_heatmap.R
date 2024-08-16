library(pheatmap)

file<-paste0("/path/to/your/matrix_markers.tsv")
matrix <- read.csv(file, header = TRUE, sep ="\t", row.names = 1, check.names = FALSE)

svgfile <-paste0("/path/to/save/heatmap.svg")
svg(file = svgfile, width = 100, height=160)

breaks1 <- c(0,seq(0.00000001, 101, by = 5))
colors <- c("white",colorRampPalette(c("#fff7e0","#604600"))(length(breaks) - 1)) #yellow
# colors <- c("white",colorRampPalette(c("#e4fdff","#015c64"))(length(breaks) - 1)) #cyan
# colors <- c("white",colorRampPalette(c("#fbfbfb","#000000"))(length(breaks) - 1)) #gray
# colors <- c("white",colorRampPalette(c("#fff3f3","#5d0000"))(length(breaks) - 1)) #red
# colors <- c("white",colorRampPalette(c("#e5ffe9","#003e08"))(length(breaks) - 1)) #green

pheatmap(matrix,cluster_row=F,cluster_cols=F,show_colnames = T, show_rownames = T, color=colors, fontsize = 70, breaks = breaks1)

dev.off()
