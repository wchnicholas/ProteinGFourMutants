#R code

CountAccessibleNodesEachPeaks <- function(l){
  l <- l[which(l==-1)]
  return (length(l))
  }

CountAccessiblePeaksEachNodes <- function(l){
  return (length(l[which(l!=-1)]))
  }

checklistitembyitem <- function(a,b){
  if (as.character(a) != as.character(b)){return ('BAD')}
  else {return ('GOOD')}
  }

EntropyCal <- function(l){
  l  <- l/1000
  en <- 0
  for (i in l){en <- en + i**2} 
  return (en)
  }

t <- read.table('analysis/LocalMaxPathLen',header=1)
t <- t[order(t$mut),]
fittable <- read.table('result/Mutfit',header=1)
fittable <- fittable[which(fittable$mut %in% t$mut & fittable$I20fit != 'NA'),]
mistable <- read.table('result/regression_missing',header=1)
allmuts  <- c(as.character(mistable$genotype_missing),as.character(fittable$mut))
allfits  <- c(exp(mistable$Var2),fittable$I20fit)
fittable <- cbind(allmuts,allfits)
fittable <- fittable[order(fittable[,1]),]

print (paste('Binding of fittable and t:', unique(mapply(checklistitembyitem,t$mut,fittable[,1]))),sep=' ')
t <- cbind(t,fittable[,2])
n <- apply(t[,2:16],1,CountAccessiblePeaksEachNodes)
t <- cbind(t,n)
t[,18] <- as.numeric(as.character(t[,18]))
colnames(t)[18] <- 'fit'
m <- read.table('analysis/LocalMaxMuts')
m <- m[order(m[,2],decreasing=1),]
moi <- t[which(t[,18]>-1),]
moi <- moi[which(!(moi$mut %in% m[,1])),]
max15 <- head(m,15)

#ANALYZING VARIANTS TO WT
tbelow1 <- t[which(t$fit < 1),]
WTno <- tbelow1[which(tbelow1$VDGV == -1),]
WTyes <- tbelow1[which(tbelow1$VDGV != -1),]
#png('analysis/LocalMaxToWTpie.png')
#pie(c(length(WTno$mut),length(WTyes$mut)),col=c('red','blue'),init.angle=90,labels=c('Inaccessible','Accessible'))
#dev.off()

#PARITIONS ALL VARIANTS BASED ON POTENTIAL EVOLUTIONARY SPACE & PLOT THE FITNESS DISTRIBUTION IN DIFFERENT PARTITIONED
png('graph/LocalMax15accessbox.png')
boxplot(moi[,18]~moi$n,pch=20,cex=0.4,cex.axis=0.9,ylim=c(0,9),xlim=c(1,16))
par(new=T)
plot(max15[,2],pch=20,ylim=c(0,9),xlim=c(0,15),col='red',type='b',axes=F)
dev.off()

#PLOT THE DISTRIBUTION OF POTENTIAL EVOLUTIONARY SPACE
png('graph/LocalMax15accesshist.png')
#hist(moi$n,breaks=50)
plot(log10(table(moi$n)),lwd=10)
dev.off()
print (table(moi$n))

#PLOT THE RELATIONSHIP BETWEEN ENTROPY AND FITNESS
des <- read.table('analysis/LocalMaxDes_random',header=1)
en <- apply(des[,3:199],1,EntropyCal)
des <- cbind(des,en)
breaks <- 100
des_highfit <- des[which(des$fit > 1),]
png('graph/LocalMaxDesReprodAll_random.png')
hist(des$en,xlim=c(0,1),breaks=breaks)
dev.off()
png('graph/LocalMaxPathReprodAll_random.png')
hist(des$PathReprod,xlim=c(0,1),breaks=breaks)
dev.off()
png('graph/LocalMaxDesReprodBenMut_random.png')
hist(des_highfit$en,xlim=c(0,1),breaks=breaks)
dev.off()
png('graph/LocalMaxPathReprodBenMut_random.png')
hist(des_highfit$PathReprod,xlim=c(0,1),breaks=breaks)
dev.off()
print (cor.test(des$fit,des$en,method='spearman'))
print (cor.test(des$fit,des$PathReprod,method='spearman'))

#ANALYZE THOSE VARIANTS THAT GOT STUCK AT LOW PEAK FITNESS
v <- moi[which(moi$n==0),]
des_v <- des[which(des$mut %in% v$mut),]
des_v <- apply(des_v[3:199],2,sum)
des_v <- cbind(attributes(des_v)$names,as.numeric(as.character(des_v)))
des_v <- des_v[which(des_v[,2]!=0),]
des_v <- des_v[order(des_v[,1]),]
fit_v <- fittable[which(fittable[,1] %in% des_v[,1]),]
fit_v <- fit_v[order(fit_v[,1]),]
print (paste('Binding of fit_v and des_v:', unique(mapply(checklistitembyitem,fit_v[,1],des_v[,1]))),sep=' ')
des_v <- cbind(des_v,fit_v[,2])
des_v[,2] <- as.numeric(as.character(des_v[,2]))
des_v[,3] <- as.numeric(as.character(des_v[,3]))
des_v <- des_v[order(as.numeric(as.character(des_v[,2])),decreasing=T),]
basin_v <- as.numeric(as.character(des_v[,2]))
png('graph/LocalMaxDes_lowpeak.png')
barplot(basin_v/sum(basin_v),names.arg=des_v[,1],col='orange',las=2,cex.names=0.5)
par(new=T)
plot(des_v[,3],axes=F,xlim=c(0.6,62.6),pch=20,col='red')
axis(side=4)
dev.off()
