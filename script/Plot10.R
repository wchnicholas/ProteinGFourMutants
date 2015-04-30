#R code
plotbasindist <- function(filename,graphname){
  t <- read.table(filename,header=1)
  png(graphname)
  plot(t$I20fit,log10(t$basin),pch=20, col='blue',ylim=c(0,5))
  abline(v=1, lty=2, lwd=1.5, col='red')
  dev.off()
  totalvarinbasingeq1 <- sum(t[which(t$I20fit>1),"basin"])/160000
  print (filename)
  print (paste('total variants included in those basin of local max > 1 is', totalvarinbasingeq1, sep=' '))
  print (paste('total local max is', length(t$mut), sep=' '))
  }

plotpathdist <- function(filename,graphname,breaks){
  t <- read.table(filename,header=1)
  png(graphname)
  hist(t$steps,xlim=c(0,35),breaks=breaks,col='orange')
  dev.off()
  print (filename)
  print (paste('Longest path has step =', max(t$steps), sep=' '))
  }

plotbasindist('analysis/LocalMaxCompile_greedy','graph/LocalMaxBasinvsFit_greedy.png')
plotbasindist('analysis/LocalMaxCompile_random','graph/LocalMaxBasinvsFit_random.png')
plotbasindist('analysis/LocalMaxCompile_weight','graph/LocalMaxBasinvsFit_weight.png')
plotpathdist('analysis/LocalMaxClimb_greedy','graph/LocalMaxPathtoMax_greedy.png',50)
plotpathdist('analysis/LocalMaxClimb_random','graph/LocalMaxPathtoMax_random.png',200)
plotpathdist('analysis/LocalMaxClimb_weight','graph/LocalMaxPathtoMax_weight.png',100)
