#R code
library(ggplot2)
library(reshape)
library(scales)

plotheatmap <- function(filename,graphname,h,w){
  mutfit <- read.table(filename, header=1, check.names=FALSE)
  y <- melt(mutfit)
  y$value[y$value==99] <- NA
  y$value[y$value> 7.5] <- 7.5
  y$value[y$value< -7.5] <- -7.5
  lvls=c('40E','40D','40R','40K','40H','40Q','40N','40S','40T','40P','40G','40C','40A','40V','40I','40L','40M','40F','40Y','40W',
         '39E','39D','39R','39K','39H','39Q','39N','39S','39T','39P','39G','39C','39A','39V','39I','39L','39M','39F','39Y','39W',
         '41E','41D','41R','41K','41H','41Q','41N','41S','41T','41P','41G','41C','41A','41V','41I','41L','41M','41F','41Y','41W',
         '54E','54D','54R','54K','54H','54Q','54N','54S','54T','54P','54G','54C','54A','54V','54I','54L','54M','54F','54Y','54W')
  y$mut <- factor(y$mut, levels=rev(lvls))
  y$variable <- factor(y$variable, levels=lvls)
  p <- ggplot(data=y, aes(x=variable, y=mut, fill=value)) +
              geom_tile() +
              scale_fill_gradientn(colours=c(colors()[114],"white",colors()[367]),
                         values=rescale(c(0, 0.5, 1)),
                         limits=c(-7.5,7.5),
                         guide="colorbar",
                         na.value="grey") +
              theme_classic() +
              theme(panel.border = element_rect(colour = "black", fill=NA, size=4),
                    text = element_text(size=20))
  ggsave(graphname,p, height=h, width=w)
  }

plotheatmap('analysis/Heatmap_G41FV54A','graph/HMG41FV54A.png',6,7)
plotheatmap('analysis/Heatmap_G41LV54H','graph/HMG41LV54H.png',6,7)
plotheatmap('analysis/Heatmap_V39WV54H','graph/HMV39WV54H.png',6,7)
