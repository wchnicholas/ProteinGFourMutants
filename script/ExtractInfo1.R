#R code
#Extract better than expected quadruple mutants
t <- read.table('result/IGG10HD4Epi',header=1)
fold <- t$mutfit/t$Maxexpfit
t <- cbind(t,fold)
t <- t[which(t$mutfit > 0.9),]
tail(t[order(t$fold),],)
