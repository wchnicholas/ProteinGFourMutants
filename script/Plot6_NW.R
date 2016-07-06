#R code
coloring <- function(v,cutoff){
  if (v > cutoff){return ('grey')}
  else (return ('red'))
  }

plotdecompose <- function(t,graph){
  png(graph,res=50,width=300,height=300)
  cutoff <- sort(t$spectrum_all_2)[length(t$genotype_all)*0.001]
  avg <- c(median(t$spectrum_all_1),median(t$spectrum_all_2),median(t$spectrum_all_3),median(t$spectrum_all_4))
  q1  <- c(quantile(t$spectrum_all_1)[2],quantile(t$spectrum_all_2)[2],quantile(t$spectrum_all_3)[2],quantile(t$spectrum_all_4)[2])
  col <- mapply(coloring,t$spectrum_all_2,rep(cutoff,length(t$spectrum_all_2)))
  for (i in 1:length(t[,1])){
    if (i%%10000 == 0){print (paste('processed',i,'variants',sep=''))}
    plot(as.numeric(t[i,2:5]),type='l',ylim=c(0,1),xlim=c(1,4),col=col[i],lwd=0.3,axes=F)
    par(new=T)
    #if (t$spectrum_all_2[i] < cutoff){print (t$genotype_all[i])}
    }
  plot(as.numeric(avg),type='l',ylim=c(0,1),xlim=c(1,4),col='blue',lwd=3)
  dev.off()
  }

plothist <- function(t,graph){
  png(graph)
  par(mfrow=c(3,1))
  hist(t$spectrum_all_1,xlim=c(0,1),breaks=100)
  hist(t$spectrum_all_2,xlim=c(0,1),breaks=25)
  hist(t$spectrum_all_3,xlim=c(0,1),breaks=10)
  dev.off()
  }

t  <- read.table('analysis/FitnessDecomposeFit',header=1)
t  <- t[order(t$spectrum_all_2,decreasing=T),]
plothist(t,'graph/FitnessDecompHist.png')
plotdecompose(t,'graph/FitnessDecomp.png')
#t2 <- t[which(t$FitRange > 2),]
#plotdecompose(t2,'graph/FitnessDecomp2.png')
