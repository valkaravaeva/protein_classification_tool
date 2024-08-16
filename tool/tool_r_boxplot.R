library(ggplot2)
library(dplyr)

ftable<-paste0("/path/to/boxplot_table.tsv")
matrix <- read.csv(ftable, header = TRUE, sep ="\t", check.names = FALSE)

pdffile <- paste0("/path/to/save/boxplot.pdf")
pdf(file = pdffile, onefile=TRUE, paper = "a4", width = 300, height=3000, compress=F)

# my_colors <- rev(c("#cb181d","#fb6a4a","#fcae91")) ##red
my_colors <- rev(c("#2171b5","#6baed6","#bdd7e7")) ##blue

mt = matrix %>% mutate(Taxon = Taxon) %>%
  group_by(Taxon) %>% 
  mutate(med_unch = median(Uncharacterized_percent))

#cornflowerblue red
ggplot(mt,aes(x = factor(Taxon),y = Uncharacterized_percent, fill=med_unch)) + geom_boxplot(outlier.colour="cornflowerblue", outlier.shape=6,outlier.size=2) + coord_flip() + scale_y_continuous(limits = c(0,1)) + scale_fill_gradientn(colors = my_colors) + theme_light()

dev.off()

