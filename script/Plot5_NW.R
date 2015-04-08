#R code

plotepidiffbg <- function(filename,graph){
  t <- read.table(filename,header=1)
  t <- t[which(t$Count >= 2),]
  print ('Total number of double variants occurs in at least two backgrounds')
  print (length(t$Mut))
  print ('Mean:')
  print (mean(t$Count))
  print ('Median:')
  print (median(t$Count))
  t <- t[which(t$MaxDFit != 'NA' & t$MinDFit != 'NA' & t$WTEpi != 'NA'),]
  print ('Total number of double variants display reciprocal sign epistasis with different signs under different backgrounds')
  print (length(t$Mut))
  EpiDiff <- t$MaxEpi/t$MinEpi
  t <- cbind(t,EpiDiff)
  t <- t[order(t$MaxEpi),]
  ylim <- c(min(log10(t$MinEpi)),max(log10(t$MaxEpi)))
  print ('range of MaxEpi-MinEpi:')
  print (ylim)
  png(graph)
  plot(log10(t$MaxEpi),col=colors()[367],ylim=ylim,pch=20,cex=0.5)
  #par(new=T)
  #plot(log10(t$WTEpi),col='green',ylim=ylim,pch=20,cex=0.5)
  par(new=T)
  plot(log10(t$MinEpi),col=colors()[114],ylim=ylim,pch=20,cex=0.5)
  dev.off()
  }

#plotepidiffbg('analysis/EpiDiffBGI10fit','graph/ContextEpiIGG10.png')
plotepidiffbg('analysis/EpiDiffBGI20fit','graph/ContextEpiIGG20.png')
#plotepidiffbg('analysis/EpiDiffBGI90fit','graph/ContextEpiIGG90.png')
