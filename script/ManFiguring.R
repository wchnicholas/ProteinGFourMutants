#R code
#Path length to the max

###############   SUBROUTINE   ###############
EntropyCal <- function(l){
  l  <- l/1000
  en <- 0
  for (i in l){en <- en + i**2}
  return (en)
  }

ColoringbyHD <- function(HD){
  if (HD == 1){return ('purple')}
  else if (HD == 2){return ('green')}
  else{return ('black')}
  }

###########   END OF SUBROUTINE   ############

plotrepeatabilityoutcome <- function(){
  w <- read.table('analysis/LocalMaxDes_weight',header=1)
  r <- read.table('analysis/LocalMaxDes_random',header=1)
  en <- apply(w[,3:199],1,EntropyCal)
  w <- cbind(w,en)
  en <- apply(r[,3:199],1,EntropyCal)
  r <- cbind(r,en)
  xlim <- c(0,1)
  ylim <- c(0,12000)
  col1 <- rgb(0,0,1,1/4)
  col2 <- rgb(1,0,0,1/4)
  #png('graph/LocalMaxDesReprodAll_random.png')
  hist(r$en,breaks=100,ylim=ylim,xlim=xlim,col=col1)
  par(new=T)
  hist(w$en,breaks=100,ylim=ylim,xlim=xlim,col=col2)
  #dev.off()
  }

plotpathlengthdistribution <- function(){
  graphname <- 'ManFig/PathLenDistribution.png'
  w <- read.table('analysis/LocalMaxDist_weight',header=1)
  r <- read.table('analysis/LocalMaxDist_random',header=1)
  g <- read.table('analysis/LocalMaxClimb_greedy',header=1)
  xlim <- c(0,round(max(c(w$avgdist2max,r$avgdist2max,g$steps))))
  ylim <- c(0,12000)
  col1 <- rgb(0,0,1,1/4)
  col2 <- rgb(1,0,0,1/4)
  col3 <- rgb(0,1,0,1/4)
  png(graphname,res=50,width=300,height=300)
  hist(w$avgdist2max,xlim=xlim,breaks=50,ylim=ylim,col=col1)
  par(new=T)
  hist(r$avgdist2max,xlim=xlim,breaks=100,ylim=ylim,col=col2,axes=F)
  par(new=T)
  hist(g$steps,xlim=xlim,breaks=50,col=col3,axes=F,ylim=c(0,41000))
  axis(4)
  #legend(0,12000,c('',''),col=c(col1,col2),pch=15,bty='n')
  dev.off()
  }

plotcorrelation <- function(){
  graphname <- 'ManFig/FitCorrelations.png'
  t <- read.table('result/Mutfit',header=1)
  t <- t[which(t$AOfit != 'NA' & t$I10fit != 'NA'),]
  png(graphname)#,res=50,width=300,height=300)
  plot(t$AOfit,t$I20fitRaw,pch=20,cex=0.5,col=mapply(ColoringbyHD,t$HD))
  dev.off()
  }

plotDFE <- function(){
  library(stringr)
  t <- read.table('result/Mutfit',header=1)
  t <- t[which(as.character(str_sub(t$mut,1,1)) != '_'
             & as.character(str_sub(t$mut,2,2)) != '_'
             & as.character(str_sub(t$mut,3,3)) != '_'
             & as.character(str_sub(t$mut,4,4)) != '_'),]
  t <- t[which(as.character(t$I10fit) != 'NA'),]
  F1 <- t[which(t$HD==1),]
  F2 <- t[which(t$HD==2),]
  F3 <- t[which(t$HD==3),]
  F4 <- t[which(t$HD==4),]
  png('ManFig/DFEIGG20.png')
  par(mfrow=c(1,4))
  ylim <- c(0,9)
  boxplot(F1$I20fit, ylim=ylim, main='HD = 1',pch=20,cex=0.5)
  boxplot(F2$I20fit, ylim=ylim, main='HD = 2',pch=20,cex=0.5)
  boxplot(F3$I20fit, ylim=ylim, main='HD = 3',pch=20,cex=0.5)
  boxplot(F4$I20fit, ylim=ylim, main='HD = 4',pch=20,cex=0.5)
  dev.off()
  }

plotEvolutionPot <- function(){
  #graphname <- 'ManFig/ShortPathLen2Ben.png'
  #t <- read.table('analysis/LocalMaxEvolvePotWT',header=1)
  graphname <- 'ManFig/ShortPathLen2Ben_pair.png'
  t <- read.table('analysis/LocalMaxEvolvePotWT_pair',header=1)
  F1 <- t[which(t$HD==1),]
  F2 <- t[which(t$HD==2),]
  F3 <- t[which(t$HD==3),]
  F4 <- t[which(t$HD==4),]
  png(graphname,res=50,width=300,height=300)
  par(mfrow=c(3,1))
  hist(F2[which(F2$PathLength!=-1),4],xlim=c(0,max(t$PathLength)),breaks=100,col='black')
  hist(F3[which(F3$PathLength!=-1),4],xlim=c(0,max(t$PathLength)),breaks=100,col='black')
  hist(F4[which(F4$PathLength!=-1),4],xlim=c(0,max(t$PathLength)),breaks=100,col='black')
  dev.off()
  }

plotImputeCor <- function(){
  library(stringr)
  graphname <- 'ManFig/ImputeCorrelation.png'
  t <- read.table('result/Mutfit',header=1)
  t <- t[which(as.character(str_sub(t$mut,1,1)) != '_'
             & as.character(str_sub(t$mut,2,2)) != '_'
             & as.character(str_sub(t$mut,3,3)) != '_'
             & as.character(str_sub(t$mut,4,4)) != '_'),]
  t <- t[which(as.character(t$I10fit) != 'NA'),]
  #t <- t[which(t$I20fit!=0),] #FOR CALCULATING CORRELATION
  r <- read.table('result/regression_all_WT',header=1)
  r <- r[which(r$genotype_comb %in% t$mut),]
  t <- t[order(t$mut),]
  r <- r[order(r$genotype_comb),]
  r1 <- r[which(t$HD==1),]
  r2 <- r[which(t$HD==2),]
  r3 <- r[which(t$HD==3),]
  r4 <- r[which(t$HD==4),]
  t1 <- t[which(t$HD==1),]
  t2 <- t[which(t$HD==2),]
  t3 <- t[which(t$HD==3),]
  t4 <- t[which(t$HD==4),]
  png(graphname)
  par(mfrow=c(2,2))
  xlim <- c(-10,2)
  ylim <- c(-10,2)
  plot(r1$fitness_comb,log(t1$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  plot(r2$fitness_comb,log(t2$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  plot(r3$fitness_comb,log(t3$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  plot(r4$fitness_comb,log(t4$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  dev.off()
  }

plotEpistasis <- function(){
  graphname <- 'ManFig/EpisDistribution.png'
  library(stringr)
  t <- read.table('result/EpistasisCol',header=1)
  t <- t[order(t$EpiE),]
  other <- t[which(t$color=='black'),]
  roi   <- t[which(t$color=='red'),]
  ylim  <- c(-5,5.5)
  png(graphname)
  plot(other$EpiE,ylim=ylim,type='l',xaxt='n',lwd=2)
  par(new=T)
  plot(roi$EpiE,ylim=ylim,type='l',col='red',axes=F,lwd=2)
  dev.off()
  }

plotCorevsfit <- function(){
  graphname <- 'ManFig/Corevsfit.png'
  t <- read.table('analysis/BulkvsFit',header=1)  
  xlim <- c(180,580)
  ylim <- c(0,9)
  png(graphname)
  plot(jitter(t$coremass,300),t$fit,cex=0.6,pch=20,xlim=xlim,ylim=ylim)
  abline(v=255.3,lty=2,col='green')
  dev.off()
  }

plotlocalmaxpredictfit <- function(){
  graphname <- 'ManFig/PreditFitLocalMax.png'
  t_lmax <- read.table('analysis/LocalMaxCompile_greedy',header=1)
  t_lmax <- t_lmax[which(t_lmax$I20fit > 1),]
  t_lmax <- t_lmax[order(t_lmax$mut),]
  r <- read.table('result/regression_all_WT',header=1)
  r <- r[which(r$genotype_comb %in% t_lmax$mut),]
  r <- r[order(r$genotype_comb),]
  t <- cbind(t_lmax,r)
  png(graphname,res=50,width=300,height=300)
  plot(t$I20fit,exp(t$fitness_comb),pch=20,xlim=c(0,9),ylim=c(0,9),cex=1)
  abline(h=exp(-0.630609982250683),col='red',lty=2)
  abline(v=1,col='blue',lty=2)
  abline(0,1,col='grey')
  dev.off()
  }

plotsubgraphpathprob <- function(){
  #graphname <- 'ManFig/SubgraphPathBias_prop.png'
  #t <- read.table('analysis/PathwayParamResult.prop',header=1)
  #png(graphname,res=50,width=300,height=300)
  graphname <- 'ManFig/SubgraphPathBias_equal.png'
  t <- read.table('analysis/PathwayParamResult.equal',header=1)
  png(graphname)
  xlim <- c(1,24)
  ylim <- c(0,1)
  for (i in 1:length(t[,1])){
    plot(as.numeric(t[i,8:31]),type='l',axes=F,xlim=xlim,ylim=ylim)
    par(new=T)
    }
  plot(seq(1/24,1,1/24),type='l',xlim=xlim,ylim=ylim,col='grey')
  dev.off()
  }

plotheatdesmap <- function(){
  library(gplots)
  library(fastcluster)
  library(marray)
  t <- read.table('analysis/LocalMaxDes_weight',header=1)
  k <- t[,3:17]
  pal<-maPalette(low="white", high="red", mid="yellow")
  heatmap.2(t(k[1:100,]),trace="none",dendrogram="none",Rowv=F,Colv=T,col=pal)
  }
