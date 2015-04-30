#R code
floor <- function(v){
  if (v < 0.01){return (0.01)}
  else{return (v)}
  }

t  <- read.table('analysis/LocalMaxEvolveAll',header=1)
t2 <- read.table('analysis/LocalMaxEvolveAll_pair',header=1)
t  <- t[which(t$AllEnds != 0),]
t2 <- t2[which(t2$AllEnds!= 0),]
t  <- t[order(t$fit,decreasing=T),]
t2 <- t2[order(t2$fit,decreasing=T),]
#t  <- t[1:20000,]
#t2 <- t2[1:20000,]
#t <- t[which(str_sub(t$mut,3,3)=='G'),]
#t2 <- t2[which(str_sub(t2$mut,3,3)=='G'),]
ylim=c(0,1)
png('graph/LocalMaxEvolvePot')
par(fig=c(0,1,0,0.7),mar=c(4.2,4.2,0.3,0.3))
plot(t$ReachEnds/t$AllEnds,cex=0.01,pch=20,col=rgb(0,0,1,1/20),ylim=ylim)
par(new=T)
plot(runmed(t$ReachEnds/t$AllEnds, 10001, endrule = "median"),col='blue',type='l',ylim=ylim)
par(new=T)
plot(t2$ReachEnds/t2$AllEnds,cex=0.01,pch=20,col=rgb(1,0,0,1/20),ylim=ylim)
par(new=T)
plot(runmed(t2$ReachEnds/t2$AllEnds, 10001, endrule = "median"),col='red',type='l',ylim=ylim)
par(fig=c(0,1,0.7,0.85),mar=c(0.3,4.2,0.3,0.3),new=TRUE)
plot(log10(mapply(floor,t2$fit)),col='red',type='l')
par(fig=c(0,1,0.85,1),mar=c(0.3,4.2,0.3,0.3),new=TRUE)
plot(log10(mapply(floor,t$fit)),col='blue',type='l')
dev.off()
