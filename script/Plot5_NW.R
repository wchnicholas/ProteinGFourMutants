#R code

plotepidiffbg <- function(filename,graph){
  t <- read.table(filename,header=1)
  t <- t[which(t$PosCount > 0 & t$NegCount > 0),]
  t <- t[order(t$EpiRange),]
  print ('Total number of double variants have both positive and negative epistasis')
  print (length(t$Mut))
  print ('Total number of double variants have both log(MaxEpi) > 0.5 and log(MinEpi) < -0.5')
  print (length(t[which(log(t$MaxEpi) > 1 & log(t$MinEpi) < -1),1]))
  print ('Mean:')
  print (mean(t$Count))
  print ('Median:')
  print (median(t$Count))
  t <- t[which(t$WTEpi != 'NA'),]
  png(graph,height=400,width=400,res=60)
  par(fig=c(0,1,0,0.7),mar=c(4.2,4.2,0.3,0.3))
  ylim <- c(-10,10)
  plot(log(t$WTEpi),col=colors()[200],ylim=ylim,pch=20,cex=0.1,las=2)
  par(new=T)
  plot(log(t$MaxEpi),col=colors()[367],ylim=ylim,pch=20,cex=0.1,axes=F)
  par(new=T)
  plot(log(t$MinEpi),col=colors()[114],ylim=ylim,pch=20,cex=0.1,axes=F)
  #abline(h=1,col=colors()[258])
  #abline(h=-1,col=colors()[258])
  par(fig=c(0,1,0.7,0.85),mar=c(0.3,4.2,0.3,0.3),new=TRUE)
  plot(t$EpiSD,col=colors()[258],pch=20,cex=0.1,xaxt='n',las=2,ylim=c(0,3))
  par(fig=c(0,1,0.85,1),mar=c(0.3,4.2,0.3,0.3),new=TRUE)
  plot(t$EpiRange,col=colors()[624],pch=20,cex=0.05,xaxt='n',las=2,ylim=c(0,15))
  dev.off()
  }

plotepidiffbg('analysis/EpiDiffBGI20fit','graph/ContextEpiIGG20.png')
