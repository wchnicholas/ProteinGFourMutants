#R code
plotting <- function(t,n,xlim,ylim,r,g,b){
  s <- t[which(t$steps==n),]
  s <- s[sample(nrow(s)),]
  med <- apply(head(s[,3:(3+n)]),2,median)
  print (paste('Total Number of', n, 'Steps Trajectroies =', length(s[,1]), sep=' '))
  for (i in 1:length(s[,1])){
    plot(0:n,s[i,3:(3+n)],type='l',ylim=ylim,xlim=xlim,col=rgb(r,g,b,0.05),lwd=1,xlab='',ylab='',axes=F)
    par(new=T)
    if (i == 100){break}
    }
  plot(0:n,med,type='l',ylim=ylim,xlim=xlim,col='red',xlab='',ylab='Entropy',lwd=1.5)
  }

t <- read.table('analysis/RepTraj_weight')
colnames(t) <- c('mut','steps','step1','step2','step3','step4',
                 'step5','step6','step7','step8','step9','step10',
                 'step11','step12','step13','step14','step15','step16',
                 'step17','step18','step19','step20')

png('graph/TrackingEn_weight.png',res=120,height=800,width=500)
par(mfrow=c(5,1),mar=c(2.2,4.2,0.3,0.3))
xlim <- c(0,9)
ylim <- c(0,2)
plotting(t,5,xlim,ylim,0.5,0.5,0.5)
plotting(t,6,xlim,ylim,0.5,0.5,0.5)
plotting(t,7,xlim,ylim,0.5,0.5,0.5)
plotting(t,8,xlim,ylim,0.5,0.5,0.5)
plotting(t,9,xlim,ylim,0.5,0.5,0.5)
dev.off()
