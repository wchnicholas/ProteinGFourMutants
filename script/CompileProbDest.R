#R code
#Analyze the output from CompileProbDest.py
plotbasindist <- function(t,graphname){
  png(graphname,width=300,height=300,res=50)
  plot(as.numeric(t[,3]),log10(as.numeric(t[,2])),pch=20, col='blue',ylim=c(0,5))
  abline(v=1, lty=2, lwd=1.5, col='red')
  dev.off()
  totalvarinbasingeq1 <- sum(as.numeric(t[which(as.numeric(t[,3])>1),2]))/160000
  print (paste('total variants included in those basin of local max > 1 is', totalvarinbasingeq1, sep=' '))
  print (paste('total local max is', length(t[,1]), sep=' '))
  }

plotpathdist <- function(t,graphname,breaks){
  png(graphname)
  hist(t$avgdist2max,xlim=c(0,max(t$avgdist2max)),breaks=breaks,col='orange')
  dev.off()
  print (paste('Longest path has step =', max(t$maxdist2max), sep=' '))
  }

t <- read.table('analysis/LocalMaxDist_random',header=1)
t <- t[which(t$maxdist2max!=0),]
png('graph/LocalMaxDisttoMax_random1000sim.png')
plot(t$fit,t$avgdist2max,pch=20,cex=0.3)
dev.off()

plotpathdist(t,'graph/LocalMaxPathtoMax_random1000sim.png',100)

t <- read.table('analysis/LocalMaxDes_random',header=1)
basin <- apply(t[,3:199],2,sum) #ADJUSTABLE
#basin <- apply(t[,3:9],2,sum) #ADJUSTABLE
basint <- cbind(attributes(basin)$names,as.numeric(basin)/1000)
print (paste('Total base is =',sum(as.numeric(basint[,2])),sep=' '))
maxfit <- read.table('analysis/LocalMaxCompile_greedy',header=1)
maxfit <- maxfit[order(maxfit$mut),]
basint <- basint[order(as.character(basint[,1])),]
outt   <- cbind(basint,maxfit$I20fit)
outt   <- outt[order(as.numeric(outt[,3]),decreasing = TRUE),]
outt   <- as.table(outt)
colnames(outt) <- c('mut','basin','I20fit')
write.table(outt,file='analysis/LocalMaxCompile_random1000sim',sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
plotbasindist(outt,'graph/LocalMaxBasinvsFit_random1000sim.png')

#Plot Accumulative probability for basin
p  <- as.numeric(as.character(outt[,2]))/160000
p <- sort(p,decreasing=T)
ap <- c()
for (i in 1:length(p)){ap <- c(ap, (sum(p[1:i])))}
png('graph/LocalMaxBasinAcProb_random1000sim.png',width=300,height=300,res=50)
plot(ap,xlim=c(1,15),type='o',ylim=c(0,1),pch=20)
dev.off()
