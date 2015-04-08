#R code
#CONVERT THIS SCRIPT TO PYTHON SCRIPT!!
flooring <- function(fit,C){
  if (fit < 0.01){return (0.01)}
  else if (fit > C){return (C*0.95)}
  else {return (fit)}
  }

fit2relKd <- function(fit,C){
  return ((C-1)/((C/flooring(fit,C))-1))
  }

#Trial 1: 25, 10, 2.5 (Not scientifically justified)
#Trial 2: 10, 5.45, 1.75 (Not the best in terms of SD, but makes sense in terms of fold diff
#Trial 2: 18.4, 7, 1.86 (The best in terms of SD, very light weighting on the fold diff

t <- read.table('result/Mutfit',header=1)
t <- t[which(t$I10fit != 'NA'),]
#t <- t[which(t$AOfit != 'NA' & t$I10fit != 'NA'),]
t <- t[which(t$I10fit != 'NA' & t$Input >= 500 & t$I10fit > 1),]
AO  <- log(mapply(fit2relKd,t$AOfit,rep(8,length(t$mut))))
I10 <- log(mapply(fit2relKd,t$I10fit,rep(18.4,length(t$mut))))
I20 <- log(mapply(fit2relKd,t$I20fit,rep(7,length(t$mut))))
I90 <- log(mapply(fit2relKd,t$I90fit,rep(1.86,length(t$mut))))
Kd_std <- mapply(function(x,y,z) sd(c(x,y,z)), I10, I20, I90)
Kd_avg <- mapply(function(x,y,z) mean(c(x,y,z)), I10, I20, I90)


#OPTIMIZATION PROCEDURE
for (i in 7){
  I20 <- log(mapply(fit2relKd,t$I20fit,rep(i,length(t$mut))))
  for (j in seq(10,20,0.05)){ 
    #I90 <- log(mapply(fit2relKd,t$I90fit,rep(j,length(t$mut))))
    #l1 <- lm(I20~I90)
    #l2 <- lm(I90~I20)
    #Kd_std <- mapply(function(x,y) sd(c(x,y)), I20, I90)
    I10 <- log(mapply(fit2relKd,t$I10fit,rep(j,length(t$mut))))
    Kd_std <- mapply(function(x,y) sd(c(x,y)), I10, I20)
    l1 <- lm(I20~I10)
    l2 <- lm(I10~I20)
    print (paste(j, median(Kd_std),as.numeric(l1$coefficients)[2],as.numeric(l2$coefficients)[2],sep=' '))
    }
  }


