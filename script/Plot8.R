#R code
t <- read.table('analysis/ThreeWayEpi',header=1)
MaxPairEpi <- mapply(max,t$Epi1vs2,t$Epi1vs3)
t <- cbind(t, MaxPairEpi)

