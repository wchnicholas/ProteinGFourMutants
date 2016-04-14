#R code
#Path length to the max

###############   SUBROUTINE   ###############
EntropyCal <- function(l){
  l  <- l/1000
  en <- 0
  for (i in l){en <- en + i**2}
  return (en)
  }

ColoringbyHD <- function(HD){
  if (HD == 1){return ('purple')}
  else if (HD == 2){return ('green')}
  else{return ('black')}
  }

###########   END OF SUBROUTINE   ############

plotrepeatabilityoutcome <- function(){
  w <- read.table('analysis/LocalMaxDes_weight',header=1)
  r <- read.table('analysis/LocalMaxDes_random',header=1)
  en <- apply(w[,3:199],1,EntropyCal)
  w <- cbind(w,en)
  en <- apply(r[,3:199],1,EntropyCal)
  r <- cbind(r,en)
  xlim <- c(0,1)
  ylim <- c(0,12000)
  col1 <- rgb(0,0,1,1/4)
  col2 <- rgb(1,0,0,1/4)
  #png('graph/LocalMaxDesReprodAll_random.png')
  hist(r$en,breaks=100,ylim=ylim,xlim=xlim,col=col1)
  par(new=T)
  hist(w$en,breaks=100,ylim=ylim,xlim=xlim,col=col2)
  #dev.off()
  }

plotpathlengthdistribution <- function(){
  graphname <- 'ManFig/PathLenDistribution.png'
  w <- read.table('analysis/LocalMaxDist_weight',header=1)
  r <- read.table('analysis/LocalMaxDist_random',header=1)
  g <- read.table('analysis/LocalMaxClimb_greedy',header=1)
  xlim <- c(0,round(max(c(w$avgdist2max,r$avgdist2max,g$steps))))
  ylim1 <- c(0,60000)
  ylim2 <- c(0,13000)
  col1 <- rgb(0,0,1,1/4)
  col2 <- rgb(1,0,0,1/4)
  col3 <- rgb(0,1,0,1/4)
  png(graphname,res=50,width=300,height=300)
  hist(g$steps,xlim=xlim,breaks=50,col=col3,ylim=ylim1)
  par(new=T)
  hist(r$avgdist2max,xlim=xlim,breaks=100,ylim=ylim2,col=col2,axes=F)
  par(new=T)
  hist(w$avgdist2max,xlim=xlim,breaks=50,ylim=ylim2,col=col1,axes=F)
  axis(4)
  #legend(0,12000,c('',''),col=c(col1,col2),pch=15,bty='n')
  dev.off()
  }

plotcorrelation <- function(){
  graphname <- 'ManFig/FitCorrelations.png'
  t <- read.table('result/Mutfit',header=1)
  t <- t[which(t$AOfit != 'NA' & t$I10fit != 'NA'),]
  png(graphname)#,res=50,width=300,height=300)
  plot(t$AOfit,t$I20fitRaw,pch=20,cex=0.5,col=mapply(ColoringbyHD,t$HD))
  dev.off()
  }

plotDFE <- function(){
  flooring <- function(fit){
    if (fit == 0){return (0.0001)}
    else {return (fit)}
    }
  library(stringr)
  library(sm)
  library(vioplot)
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
  png('ManFig/DFEIGG20.png',res=65,width=408,height=360)
  #ylim <- c(0,9)
  ylim <- c(-10,4)
  vioplot(log(mapply(flooring,F1$I20fit)), log(mapply(flooring,F2$I20fit)), 
          log(mapply(flooring,F3$I20fit)), log(mapply(flooring,F4$I20fit)),ylim=ylim,col='gold',range=0)
  dev.off()
  #boxplot(log(mapply(flooring,c(F1$I20fit, F2$I20fit, F3$I20fit, F4$I20fit)))~
  #        c(rep('HD1',length(F1$I20fit)),rep('HD2',length(F2$I20fit)),rep('HD3',length(F3$I20fit)),rep('HD4',length(F4$I20fit))),
  #        pch=20,cex=0.5,ylim=ylim)
  }

plotEvolutionPot <- function(){
  graphname <- 'ManFig/ShortPathLen2Ben.png'
  t <- read.table('analysis/LocalMaxEvolvePotWT',header=1)
  #graphname <- 'ManFig/ShortPathLen2Ben_pair.png'
  #t <- read.table('analysis/LocalMaxEvolvePotWT_pair',header=1)
  #png(graphname,res=50,width=300,height=300)
  png(graphname)
  F1 <- t[which(t$HD==1),]
  F2 <- t[which(t$HD==2),]
  F3 <- t[which(t$HD==3),]
  F4 <- t[which(t$HD==4),]
  par(mfrow=c(3,1))
  hist(F2[which(F2$PathLength!=-1),4],xlim=c(0,max(t$PathLength)),breaks=50,col='black')
  hist(F3[which(F3$PathLength!=-1),4],xlim=c(0,max(t$PathLength)),breaks=100,col='black')
  hist(F4[which(F4$PathLength!=-1),4],xlim=c(0,max(t$PathLength)),breaks=100,col='black')
  dev.off()
  }

plotImputeCor <- function(){
  library(stringr)
  graphname <- 'ManFig/ImputeCorrelation.png'
  t <- read.table('result/Mutfit',header=1)
  t <- t[which(as.character(str_sub(t$mut,1,1)) != '_'
             & as.character(str_sub(t$mut,2,2)) != '_'
             & as.character(str_sub(t$mut,3,3)) != '_'
             & as.character(str_sub(t$mut,4,4)) != '_'),]
  t <- t[which(as.character(t$I10fit) != 'NA'),]
  #t <- t[which(t$I20fit!=0),] #FOR CALCULATING CORRELATION
  r <- read.table('result/regression_all_WT',header=1)
  r <- r[which(r$genotype_comb %in% t$mut),]
  t <- t[order(t$mut),]
  r <- r[order(r$genotype_comb),]
  r1 <- r[which(t$HD==1),]
  r2 <- r[which(t$HD==2),]
  r3 <- r[which(t$HD==3),]
  r4 <- r[which(t$HD==4),]
  t1 <- t[which(t$HD==1),]
  t2 <- t[which(t$HD==2),]
  t3 <- t[which(t$HD==3),]
  t4 <- t[which(t$HD==4),]
  png(graphname)
  par(mfrow=c(2,2))
  xlim <- c(-10,2)
  ylim <- c(-10,2)
  plot(r1$fitness_comb,log(t1$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  plot(r2$fitness_comb,log(t2$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  plot(r3$fitness_comb,log(t3$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  plot(r4$fitness_comb,log(t4$I20fit),pch=20,cex=0.1,xlim=xlim,ylim=ylim)
  dev.off()
  }

plotEpistasis <- function(){
  graphname <- 'ManFig/EpisDistribution.png'
  library(stringr)
  t <- read.table('result/EpistasisCol',header=1)
  t <- t[order(t$EpiE),]
  other <- t[which(t$color=='black'),]
  roi   <- t[which(t$color=='red'),]
  ylim  <- c(-5,5.5)
  png(graphname)
  plot(other$EpiE,ylim=ylim,type='l',xaxt='n',lwd=2)
  par(new=T)
  plot(roi$EpiE,ylim=ylim,type='l',col='red',axes=F,lwd=2)
  dev.off()
  }

plotCorevsfit <- function(){
  flooring <- function(fit){
    if (fit < 0.01){return (0.01)}
    else {return (fit)}
    }

  computefracben <- function(t,tmass,masspoint,span){
    fracben    <- c()
    for (i in masspoint){
      p       <- t[which(tmass >= i-span & tmass <= i+span),]
      fracben <- c(fracben,length(p$fit[which(p$fit > 1)])/length(p$fit))
      }
    return (fracben)
    }

  plotting <- function(mass,masspoint,fit,fracben,WTmass,xlim,ylim,graphname){
    png(graphname)
    plot(mass,fit,cex=0.6,pch=20,xlim=xlim,ylim=ylim)
    abline(v=WTmass,lty=2,col='blue',lwd=2)
    par(new=T)
    plot(masspoint,fracben,type='l',xlim=xlim,col='red',yaxt='n',lwd=2)
    axis(4)
    dev.off()
    }

  #graphname <- 'ManFig/Corevsfit.png'
  graphname <- 'ManFig/CorevsfitHD4.png'
  t <- read.table('analysis/BulkvsFit',header=1)  
  #t <- t[which(t$HD==4),]
  masspoint_core  <- seq(min(t$coremass),max(t$coremass),1)
  masspoint_res39 <- sort(unique(t$mass39))
  masspoint_res40 <- sort(unique(t$mass40))
  masspoint_res41 <- sort(unique(t$mass41))
  masspoint_res54 <- sort(unique(t$mass54))
  fracben_core   <- computefracben(t,t$coremass,masspoint_core,20)
  fracben_res39  <- computefracben(t,t$mass39,masspoint_res39,0.5)
  fracben_res40  <- computefracben(t,t$mass40,masspoint_res40,0.5)
  fracben_res41  <- computefracben(t,t$mass41,masspoint_res41,0.5)
  fracben_res54  <- computefracben(t,t$mass54,masspoint_res54,0.5)
  plotting(t$coremass,masspoint_core,log(mapply(flooring,t$fit)),fracben_core,258,c(150,500),c(-5,3),graphname)
  }

plotddgvsfit <- function(){
  flooring <- function(fit){
    if (fit < 0.01){return (0.01)}
    else {return (fit)}
    }
  plotddgvscore <- function(t,graphname2){
    png(graphname2)
    regline <- lm(t$ddg~t$coremass)
    plot(t$coremass,t$ddg,cex=0.1,pch=20,ylim=c(0,50),xlim=c(150,500))
    abline(regline,col='purple',lwd=1.5)
    dev.off()
    print (summary(regline))
    }
  graphname1 <- 'ManFig/DDGvsfit.png'
  graphname2 <- 'ManFig/DDGvsMass.png'
  graphname3 <- 'ManFig/DDGvsMassHD4.png'
  t <- read.table('analysis/BulkvsFit',header=1)
  d <- read.table('result/1PGAddg')
  t <- t[which(t$mut %in% d[,1]),]
  d <- d[which(d[,1] %in% t$mut),]
  t <- t[order(t$mut),]
  d <- d[order(d[,1]),]
  t <- cbind(t,d[,2])
  colnames(t)[9] <- 'ddg'
  t$fit <- mapply(flooring,t$fit)
  HD1 <- t[which(t$HD == 1),]
  HD2 <- t[which(t$HD == 2),]
  HD3 <- t[which(t$HD == 3),]
  HD4 <- t[which(t$HD == 4),]
  #t <- t[which(t$fit > 0.01),]
  #t <- t[which(t$fit > 0.01 & t$HD == 4),]

  png(graphname1)
  regline <- lm(t$ddg~log(t$fit))
  plot(log(t$fit),t$ddg,cex=0.3,pch=20)
  abline(regline,col='purple',lwd=1.5)
  dev.off()
  
  plotddgvscore(t,graphname2)
  plotddgvscore(HD4,graphname3)
  }

plotlocalmaxpredictfit <- function(){
  graphname <- 'ManFig/PreditFitLocalMax.png'
  t_lmax <- read.table('analysis/LocalMaxCompile_greedy',header=1)
  t_lmax <- t_lmax[which(t_lmax$I20fit > 1),]
  t_lmax <- t_lmax[order(t_lmax$mut),]
  r <- read.table('result/regression_all_WT',header=1)
  r <- r[which(r$genotype_comb %in% t_lmax$mut),]
  r <- r[order(r$genotype_comb),]
  t <- cbind(t_lmax,r)
  png(graphname,res=50,width=300,height=300)
  plot(t$I20fit,exp(t$fitness_comb),pch=20,xlim=c(0,9),ylim=c(0,9),cex=1)
  abline(h=exp(-0.630609982250683),col='red',lty=2)
  abline(v=1,col='blue',lty=2)
  abline(0,1,col='grey')
  dev.off()
  }

plotsubgraphpathprob <- function(){
  plotting <- function(t){
    xlim <- c(1,24)
    ylim <- c(0,1)
    for (i in 1:length(t[,1])){
      plot(as.numeric(t[i,8:31]),type='l',axes=F,xlim=xlim,ylim=ylim,lwd=0.4)
      par(new=T)
      }
    plot(seq(1/24,1,1/24),type='l',xlim=xlim,ylim=ylim,col='purple',lwd=2)
    }

  graphname  <- 'ManFig/SubgraphPathBias_prop.png'
  t <- read.table('analysis/PathwayParamResult.prop',header=1)
  png(graphname)#,res=50,width=300,height=300)
  plotting(t)
  dev.off()

  graphname  <- 'ManFig/SubgraphPathBias_equal.png'
  t <- read.table('analysis/PathwayParamResult.equal',header=1)
  png(graphname)
  plotting(t)
  dev.off()

  graphname  <- 'ManFig/SubgraphPathBiasToWT_prop.png'
  t <- read.table('analysis/PathwayParamResultToWT.prop',header=1)
  png(graphname)
  plotting(t)
  dev.off()

  graphname  <- 'ManFig/SubgraphPathBiasToWT_equal.png'
  t <- read.table('analysis/PathwayParamResultToWT.equal',header=1)
  png(graphname)
  plotting(t)
  dev.off()
  }

plotginihist <- function(){
  barplotting <- function(graphname,t,res,width,height){
    png(graphname,res=res,width=width,height=height)
    barplot(sort(t$giniindex,decreasing=T),col='black',ylab='gini index',space=1.1)
    dev.off()
    }

  histplotting <- function(graphname,t,res,width,height,breaks){
    png(graphname,res=res,width=width,height=height)
    hist(sort(t$giniindex,decreasing=T),col='grey',xlab='gini index',ylab='Frequency',breaks=breaks,xlim=c(0,1))
    dev.off()
    }

  t <- read.table('analysis/PathwayParamResult.prop',header=1)
    graphname  <- 'ManFig/GiniHist_prop.png'
    #barplotting(graphname,t,50,230,230)
    histplotting(graphname,t,50,230,230,5)
  t <- read.table('analysis/PathwayParamResult.equal',header=1)
    graphname  <- 'ManFig/GiniHist_equal.png'
    #barplotting(graphname,t,50,230,230)
    histplotting(graphname,t,50,230,230,10)
    graphname2 <- 'ManFig/CountPathHist.png'
    png(graphname2)
    barplot(sort(t$monopaths,decreasing=T),col='black',ylab='gini index',space=1.1)
    dev.off()

  t <- read.table('analysis/PathwayParamResultToWT.prop',header=1)
    graphname  <- 'ManFig/GiniHistToWT_prop.png'
    histplotting(graphname,t,50,230,230,10)
  t <- read.table('analysis/PathwayParamResultToWT.equal',header=1)
    graphname  <- 'ManFig/GiniHistToWT_equal.png'
    histplotting(graphname,t,50,230,230,20)
    graphname2 <- 'ManFig/CountPathHistToWT.png'
    png(graphname2)
    barplot(table(c(t$monopaths,seq(1,24,1)))-1,col='grey')
    dev.off()

  graphname3 <- 'ManFig/CountPathBox.png'
    t1 <- read.table('analysis/PathwayParamResult.equal',header=1)
    t2 <- read.table('analysis/PathwayParamResultToWT.equal',header=1)
    png(graphname3,width=350,height=400)
    boxplot(c(t1$monopaths,t2$monopaths)~c(rep(1,length(t1$monopaths)),rep(2,length(t2$monopaths))),pch=20,xaxt='n')
    dev.off()
  }

plotheatdesmap <- function(){
  library(gplots)
  library(fastcluster)
  library(marray)
  t <- read.table('analysis/LocalMaxDes_weight',header=1)
  k <- t[,3:17]
  pal<-maPalette(low="white", high="red", mid="yellow")
  heatmap.2(t(k[1:100,]),trace="none",dendrogram="none",Rowv=F,Colv=T,col=pal)
  }

plotmonopathdist <- function(){
  graphname <- 'ManFig/DistMonoSubgraphBen4.png'
  t <- read.table('analysis/PathwayParamResult.prop',header=1)
  png(graphname,res=50,width=300,height=300)
  plot(sort(t$monopaths),type='h',lwd=10)
  dev.off()
  }

plotlocalmaxdist( <- function(){
  graphname <- 'ManFig/DistLocalMax.png'
  t <- read.table('SciMimic/SciMimic_WTbenNoFillin_WorkSpace.tsv',header=1)
  m <- read.table('analysis/LocalMaxCompile_greedy',header=1)
  m <- m[which(m$I20fit>1),]
  t <- t[which(t$Id %in% m$mut),]
  p <- table(c(t$Modularity,0:6))-1
  png(graphname,res=60,width=350,height=250)
  barplot(p,col=c('red','orange','yellow','green','cyan','blue','purple'))
  dev.off()
  }

plotSciMimicPath <- function(){
  vioplotting <- function(t,col,graphname,switch){
    library(sm)
    library(vioplot)
    Class1 <- t[which(t$mut %in% read.table('SciMimic/Class1Vars')[,1]),]
    Class2 <- t[which(t$mut %in% read.table('SciMimic/Class2Vars')[,1]),]
    Class3 <- t[which(t$mut %in% read.table('SciMimic/Class3Vars')[,1]),]
    Class4 <- t[which(t$mut %in% read.table('SciMimic/Class4Vars')[,1]),]
    Class5 <- t[which(t$mut %in% read.table('SciMimic/Class5Vars')[,1]),]
    Class6 <- t[which(t$mut %in% read.table('SciMimic/Class6Vars')[,1]),]
    Class7 <- t[which(t$mut %in% read.table('SciMimic/Class7Vars')[,1]),]
    yes <- rbind(Class1,Class2,Class4,Class6)
    no  <- rbind(Class3,Class5,Class7)
    if (switch == 1){yes <- yes$avgdist2max; no <- no$avgdist2max}
    if (switch == 2){yes <- yes$steps; no <- no$steps}
    png(graphname,res=50,width=120,height=240)
    ylim <- c(min(yes,no),max(yes,no))
    vioplot(yes,no,ylim=ylim,col=col,range=0)
    dev.off()
    print (mean(no))
    print (mean(yes))
    ttest <-  t.test(yes,no)
    print (ttest)
    }
  w <- read.table('analysis/LocalMaxDist_weight',header=1)
  r <- read.table('analysis/LocalMaxDist_random',header=1)
  g <- read.table('analysis/LocalMaxClimb_greedy',header=1) 
  col1 <- rgb(0,0,1,1/4) #Weight
  col2 <- rgb(1,0,0,1/4) #Random
  col3 <- rgb(0,1,0,1/4) #Greedy
  vioplotting(w,col1,'ManFig/SciMimicPathLen_Weight.png',1)
  vioplotting(r,col2,'ManFig/SciMimicPathLen_Random.png',1)
  vioplotting(g,col3,'ManFig/SciMimicPathLen_Greedy.png',2)
  }

plotbootstrap15 <- function(){
  graphname <- 'ManFig/Random15HD.png'
  l <- scan('analysis/RandomSampled15_All')
  png(graphname,res=60,width=400,height=400)
  hist(l,breaks=100,xlim=c(3.3,max(l)))
  abline(v=3.4,col='purple',lwd=2)
  dev.off()
  }
