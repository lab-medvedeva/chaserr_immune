library(ggplot2)
library(ggpubr)
library(png)
library(gridExtra)
library(reshape2)
library(AnnotationHub)
library(biomaRt)
library(magick)

# A. ChIP-seq Human FOXP3 from https://chip-atlas.org/
image <- image_read('chipseq_human.png')
plot <- image_ggplot(image)
p1 <- plot+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.tag = element_text(size = 8, face="plain")) +
    labs(title = "", tag = "A")

# B. ChIP-seq from https://doi.org/10.1038/s41590-018-0291-z and https://doi.org/10.1038/s41590-023-01685-w
image <- image_read('chipseq_mouse.png') 
plot <- image_ggplot(image)
p2 <- plot+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.text.y =  element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.tag = element_text(size = 8, face="plain")) +
    labs(title = "", tag = "B")

g <- arrangeGrob(p1, p2, ncol = 1, nrow = 2, layout_matrix= rbind(c(1), c(2)))
ggsave(paste0("foxp1_foxp3_supplementary.png"), plot = g, units = "mm", width = 170, height = 100)
