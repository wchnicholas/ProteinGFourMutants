#R code
filename <- 'result/HD4EpiIGG20'
t <- read.table(filename,header=1)
t <- t[order(t$mutfit),]
eoi  <- t[which(t$mutfit > 1),]
fold <- eoi$mutfit/eoi$Maxexpfit
eoi  <- cbind(eoi,fold)
eoi  <- eoi[order(eoi$fold),]
plot(log10(t$mutfit),log10(t$Maxexpfit),pch=20,cex=0.2)
