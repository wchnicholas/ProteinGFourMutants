#R code
t_prop  <- read.table('analysis/PathwayParamResult.prop',header=1)
t_equal <- read.table('analysis/PathwayParamResult.equal',header=1)
xlim    <- c(1,24)
ylim    <- c(0,1)
plot(t_prop$monopaths,t_prop$giniindex,xlim=xlim,ylim=ylim,col='blue',pch=20)
par(new=T)
plot(t_equal$monopaths,t_equal$giniindex,xlim=xlim,ylim=ylim,col='green',pch=20)
