#R code
library(stringr)
flooring <- function(fit){
  if (fit < 0.01){return (0.01)}
  else {return (fit)}
  }

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

png('graph/DFEboxIGG10.png')
par(mfrow=c(1,4))
ylim <- c(0,30)
boxplot(F1$I10fit, ylim=ylim, main='HD = 1')
boxplot(F2$I10fit, ylim=ylim, main='HD = 2')
boxplot(F3$I10fit, ylim=ylim, main='HD = 3')
boxplot(F4$I10fit, ylim=ylim, main='HD = 4')
dev.off()

png('graph/DFEboxIGG20.png')
par(mfrow=c(1,4))
ylim <- c(0,30)
ylim <- c(-5,4)
boxplot(log(mapply(flooring,F1$I20fit)), ylim=ylim, main='HD = 1')
boxplot(log(mapply(flooring,F2$I20fit)), ylim=ylim, main='HD = 2')
boxplot(log(mapply(flooring,F3$I20fit)), ylim=ylim, main='HD = 3')
boxplot(log(mapply(flooring,F4$I20fit)), ylim=ylim, main='HD = 4')
dev.off()

png('graph/DFEboxIGG90.png')
par(mfrow=c(1,4))
ylim <- c(0,30)
boxplot(F1$I90fit, ylim=ylim, main='HD = 1')
boxplot(F2$I90fit, ylim=ylim, main='HD = 2')
boxplot(F3$I90fit, ylim=ylim, main='HD = 3')
boxplot(F4$I90fit, ylim=ylim, main='HD = 4')
dev.off()

unused <- function(){
#Plot DFE below 1
F1_I10 <- F1[which(F1$I10fit <= 1),]
F2_I10 <- F2[which(F2$I10fit <= 1),]
F3_I10 <- F3[which(F3$I10fit <= 1),]
F4_I10 <- F4[which(F4$I10fit <= 1),]
png('graph/DFEhistIGG10.png')
par(mfrow=c(4,1))
xlim=c(0,1)
hist(F1_I10$I10fit,xlim=xlim)
hist(F2_I10$I10fit,xlim=xlim)
hist(F3_I10$I10fit,xlim=xlim)
hist(F4_I10$I10fit,xlim=xlim)
dev.off()

#Plot DFE below 1
F1_I20 <- F1[which(F1$I20fit <= 1),]
F2_I20 <- F2[which(F2$I20fit <= 1),]
F3_I20 <- F3[which(F3$I20fit <= 1),]
F4_I20 <- F4[which(F4$I20fit <= 1),]
png('graph/DFEhistIGG20.png')
par(mfrow=c(4,1))
hist(F1_I20$I20fit,xlim=xlim)
hist(F2_I20$I20fit,xlim=xlim)
hist(F3_I20$I20fit,xlim=xlim)
hist(F4_I20$I20fit,xlim=xlim)
dev.off()

#Plot DFE below 1
F1_I90 <- F1[which(F1$I90fit <= 1),]
F2_I90 <- F2[which(F2$I90fit <= 1),]
F3_I90 <- F3[which(F3$I90fit <= 1),]
F4_I90 <- F4[which(F4$I90fit <= 1),]
png('graph/DFEhistIGG90.png')
par(mfrow=c(4,1))
hist(F1_I90$I90fit,xlim=xlim)
hist(F2_I90$I90fit,xlim=xlim)
hist(F3_I90$I90fit,xlim=xlim)
hist(F4_I90$I90fit,xlim=xlim)
dev.off()
}
