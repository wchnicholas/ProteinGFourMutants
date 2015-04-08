#R code

plotting <-function(t,col,ylim){
  plot(log10(t$MaxFit),pch=20,col=col,ylim=ylim,type='l')
  par(new=T)
  plot(log10(t$MinFit),pch=20,col=col,ylim=ylim,type='l')
  }

realminmax <- function(filename,graph){
  t  <- read.table(filename,header=1)
  t1 <- t[which(t$HD == 1),]
  t2 <- t[which(t$HD == 2),]
  t3 <- t[which(t$HD == 3),]
  t4 <- t[which(t$HD == 4),]
  t1 <- t1[order(t1$Mut),]
  t2 <- t2[order(t2$Mut),]
  t3 <- t3[order(t3$Mut),]
  t4 <- t4[order(t4$Mut),]
  t4 <- t4[order(t1$MaxFit),]
  t3 <- t3[order(t1$MaxFit),]
  t2 <- t2[order(t1$MaxFit),]
  t1 <- t1[order(t1$MaxFit),] #Do not change the order of t1 first, change t1 last
  t4 <- t4[which(t1$MaxFit != 'NA'),]
  t3 <- t3[which(t1$MaxFit != 'NA'),]
  t2 <- t2[which(t1$MaxFit != 'NA'),]
  t1 <- t1[which(t1$MaxFit != 'NA'),]
  png(graph)
  ylim <- c(-3,3)
  plotting(t1,'green',ylim); par(new=T)
  plotting(t2,'blue',ylim); par(new=T)
  plotting(t3,'orange',ylim); par(new=T)
  plotting(t4,'red',ylim)
  dev.off()
  #t.test(t2$MaxFit,t3$MaxFit,paired=T) #T test
  }

realminmax('analysis/MutDiffBGI20fit','graph/DiffBGIGG20.png')

#FITNESS DISTRIBUTION: SINGLE MUTATION ON DIFFERENT HD
t  <- read.table('analysis/FitDiffBG',header=1)
t  <- t[which(t$condition == 'I20fit'),]
t1 <- t[which(t$HD == 1),]
t2 <- t[which(t$HD == 2),]
t3 <- t[which(t$HD == 3),]
t4 <- t[which(t$HD == 4),]
xlim <- log10(c(min(t$fit),max(t$fit)))
png('graph/DFEDiffBG.png')
par(mfrow=c(4,1))
hist(log10(t1$fit),xlim=xlim,breaks=20)
hist(log10(t2$fit),xlim=xlim,breaks=40)
hist(log10(t3$fit),xlim=xlim,breaks=40)
hist(log10(t4$fit),xlim=xlim,breaks=40)
dev.off()
