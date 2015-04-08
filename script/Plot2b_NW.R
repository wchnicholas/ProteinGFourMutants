#R code

coloring <- function(HD){
  if (HD == 4){return (rgb(1,0,0,0.5))}
  if (HD == 3){return (rgb(1,0.64,0,0.5))}
  if (HD == 2){return (rgb(0,0,1,0.5))}
  if (HD == 1){return (rgb(0,1,0,0.5))}
  if (HD == 0){return (rgb(0,0,0,0.5))}
  }

t <- read.table('result/AllEpi',header=1)
t <- t[which(t$cap != 'NA'),]
IGG10 <- t[which(t$condition=='I10fit'),]
IGG20 <- t[which(t$condition=='I20fit'),]
IGG90 <- t[which(t$condition=='I90fit'),]
IGG10 <- IGG10[order(IGG10$mut),]
IGG20 <- IGG20[order(IGG20$mut),]
IGG90 <- IGG90[order(IGG90$mut),]
col   <- mapply(coloring,as.numeric(IGG10$HD))

#Plotting Comparison of Epistasis in different conditions
ylim <- c(-6,8)
xlim <- c(-6,8)
png('graph/DEEIGG90vs10.png')
plot(log(IGG10$epis),log(IGG90$epis),col=col,pch=20,xlim=xlim,ylim=ylim,cex=0.2)
abline(0,1,col='black',lty=2)
dev.off()

png('graph/DEEIGG90vs20.png')
plot(log(IGG20$epis),log(IGG90$epis),col=col,pch=20,xlim=xlim,ylim=ylim,cex=0.2)
abline(0,1,col='black',lty=2)
dev.off()

#Plotting DEE (distribution of episatic effect)
png('graph/DEEhistIGG10.png')
par(mfrow=c(1,3))
ylim <- log(c(min(t$epis),max(t$epis)))
boxplot(log(IGG10[which(IGG10$HD==2),]$epis),ylim=ylim)
boxplot(log(IGG10[which(IGG10$HD==3),]$epis),ylim=ylim)
boxplot(log(IGG10[which(IGG10$HD==4),]$epis),ylim=ylim)
dev.off()

png('graph/DEEhistIGG20.png')
par(mfrow=c(1,3))
ylim <- log(c(min(t$epis),max(t$epis)))
boxplot(log(IGG20[which(IGG20$HD==2),]$epis),ylim=ylim)
boxplot(log(IGG20[which(IGG20$HD==3),]$epis),ylim=ylim)
boxplot(log(IGG20[which(IGG20$HD==4),]$epis),ylim=ylim)
dev.off()

png('graph/DEEhistIGG90.png')
par(mfrow=c(1,3))
ylim <- log(c(min(t$epis),max(t$epis)))
boxplot(log(IGG90[which(IGG90$HD==2),]$epis),ylim=ylim)
boxplot(log(IGG90[which(IGG90$HD==3),]$epis),ylim=ylim)
boxplot(log(IGG90[which(IGG90$HD==4),]$epis),ylim=ylim)
dev.off()

png('graph/ExpfitboxIGG10.png')
par(mfrow=c(1,4))
ylim <- c(min(t$expfit),max(t$expfit))
boxplot(IGG10[which(IGG10$HD==1),]$expfit, ylim=ylim, main='HD = 1')
boxplot(IGG10[which(IGG10$HD==2),]$expfit, ylim=ylim, main='HD = 2')
boxplot(IGG10[which(IGG10$HD==3),]$expfit, ylim=ylim, main='HD = 3')
boxplot(IGG10[which(IGG10$HD==4),]$expfit, ylim=ylim, main='HD = 4')
dev.off()

png('graph/ExpfitboxIGG20.png')
par(mfrow=c(1,4))
ylim <- c(min(t$expfit),max(t$expfit))
boxplot(IGG20[which(IGG20$HD==1),]$expfit, ylim=ylim, main='HD = 1')
boxplot(IGG20[which(IGG20$HD==2),]$expfit, ylim=ylim, main='HD = 2')
boxplot(IGG20[which(IGG20$HD==3),]$expfit, ylim=ylim, main='HD = 3')
boxplot(IGG20[which(IGG20$HD==4),]$expfit, ylim=ylim, main='HD = 4')
dev.off()

png('graph/ExpfitboxIGG90.png')
par(mfrow=c(1,4))
ylim <- c(min(t$expfit),max(t$expfit))
boxplot(IGG90[which(IGG90$HD==1),]$expfit, ylim=ylim, main='HD = 1')
boxplot(IGG90[which(IGG90$HD==2),]$expfit, ylim=ylim, main='HD = 2')
boxplot(IGG90[which(IGG90$HD==3),]$expfit, ylim=ylim, main='HD = 3')
boxplot(IGG90[which(IGG90$HD==4),]$expfit, ylim=ylim, main='HD = 4')
dev.off()
