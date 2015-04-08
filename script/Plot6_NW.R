#R code
plotdecompose <- function(t){
  png('graph/FitnessDecomp.png')
  avg <- c(mean(t$pearson_all_1),mean(t$pearson_all_2),mean(t$pearson_all_3),mean(t$pearson_all_4))
  for (i in 1:length(t[,1])){
    if (i%%10000 == 0){print (i)}
    plot(as.numeric(t[i,2:5]),type='l',ylim=c(0,1),xlim=c(1,4),col='grey',lwd=0.3,axes=F)
    par(new=T)
    }
  plot(as.numeric(avg),type='l',ylim=c(0,1),xlim=c(1,4),col='red',lwd=0.8)
  dev.off()
  }

plothist <- function(t){
  png('graph/FitnessDecompHist.png')
  par(mfrow=c(3,1))
  hist(t$pearson_all_1,xlim=c(0,1),breaks=100)
  hist(t$pearson_all_2,xlim=c(0,1),breaks=25)
  hist(t$pearson_all_3,xlim=c(0,1),breaks=10)
  dev.off()
  }

t <- read.table('analysis/FitnessDecompose',header=1)
plothist(t)
plotdecompose(t)
