#R code

parsetable <- function(t){
  Fstuck <- t$StuckPath/t$Paths
  FmonoI <- t$MonoIPath/t$Paths
  t <- cbind(t, Fstuck)
  t <- cbind(t, FmonoI)
  return (t)
  }


t10 <- parsetable(read.table('analysis/ShortPaths4mutsCutoff10fold',header=1))
t1  <- parsetable(read.table('analysis/ShortPaths4mutsCutoff1fold',header=1))

png('graph/PathwayDistribute.png')
par(mfrow=c(2,1))
xlim <- c(0,25)
hist(t10$Paths,breaks=100,col='black',xlim=xlim)
hist(t1$Paths,breaks=60,col='black',xlim=xlim)
dev.off()

png('graph/PathwayDistribute_stuck.png') #Fraction of paths that got stuck
par(mfrow=c(2,1))
xlim <- c(0,1)
hist(t10$Fstuck,breaks=100,col='black',xlim=xlim)
hist(t1$Fstuck,breaks=100,col='black',xlim=xlim)
dev.off()

png('graph/PathwayDistribute_mono.png') #Fraction of paths being monotonic
par(mfrow=c(2,1))
xlim <- c(0,1)
hist(t10$FmonoI,breaks=50,col='black',xlim=xlim)
hist(t1$FmonoI,breaks=100,col='black',xlim=xlim)
dev.off()
