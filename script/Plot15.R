#R code
library(ggplot2)
library(reshape)
library(scales)

t <- read.table('analysis/AAtransitionmatrix',header=1)
t$value[t$value==-1] <- NA
p <- ggplot(data=t, aes(x=aa1, y=aa2, fill=value)) +
            geom_tile(colour = "black") +
            scale_fill_gradientn(colours=c("white", "yellow", "red"),
                       values=rescale(c(0, 0.5, 1)),
                       guide='none',
                       na.value="grey") +
            theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                  text = element_text(size=6)) +
            xlab('Initial Amino Acid') + 
            ylab('Targeting Amino Acid')
ggsave('ManFig/AAcodontransit.png',p, height=1.8, width=1.8)
