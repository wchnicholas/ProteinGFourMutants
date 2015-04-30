#R code
flooring <- function(fit,C){
  if (fit < 0.01){return (0.01)}
  else if (fit > C){return (C*0.95)}
  else {return (fit)}
  }

fit2relKd <- function(fit,C){
  return ((C-1)/((C/flooring(fit,C))-1))
  }

relKd2fit <- function(Kd,C){
  return ((Kd*C)/(Kd+C-1))
  }

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
  }

#PLOTTING INTER-CORRELATION AGAINST ANDERS DATA FOR ALL SINGLE AND DOUBLE MUTANTS
t <- read.table('result/Mutfit',header=1)
AO <- t[which(t$AOfit != 'NA' & t$I20fitRaw != 'NA'),]
AO_sub <- cbind(AO$AOfit, AO$I10fit, AO$I20fitRaw, AO$I90fit)
colnames(AO_sub) <- c('AOfit', 'I10fit', 'I20fitRaw', 'I90fit')
png('graph/G4Cor_AO.png')
#pairs(AO_sub, pch=20, lower.panel=panel.smooth, upper.panel=panel.cor,
pairs(~AOfit+I10fit+I20fitRaw+I90fit,data=AO_sub, pch=20, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()

#PLOTTING INTER-CORRELATION BETWEEN DIFFERENT CONDITIONS FOR ALL MUTANTS
t <- t[which(t$I20fitRaw != 'NA'),]
fittable <- cbind(t$I10fit, t$I20fitRaw, t$I90fit)
colnames(fittable) <- c('I10fit', 'I20fitRaw', 'I90fit')
png('graph/G4Cor_NW.png')
pairs(fittable, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()

#Using Kd to convert to fitness under other conditions
#plot(AO$I90fit, mapply(relKd2fit,mapply(fit2relKd,AO$I20,10),2))
