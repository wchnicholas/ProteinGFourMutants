#R code
#Put different parameters in to a vector, order = Input, IGG10, IGG20, IGG90
library(stringr)

countnumvar <- function(countlist,n){return (length(countlist[which(countlist>=n)]))}

basic_stats <- function(t,main){
  totalread   <- c(sum(t$Input),sum(t$IGG10),sum(t$IGG20),sum(t$IGG90))
  totalnumvar <- c(countnumvar(t$Input,1),countnumvar(t$IGG10,1),countnumvar(t$IGG20,1),countnumvar(t$IGG90,1))
  geq10numvar <- c(countnumvar(t$Input,10),countnumvar(t$IGG10,10),countnumvar(t$IGG20,10),countnumvar(t$IGG90,10))
  maxcountvar <- c(max(t$Input),max(t$IGG10),max(t$IGG20),max(t$IGG90))
  mincountvar <- c(min(t$Input),min(t$IGG10),min(t$IGG20),min(t$IGG90))
  avgcountvar <- c(mean(t$Input),mean(t$IGG10),mean(t$IGG20),mean(t$IGG90))
  medcountvar <- c(median(t$Input),median(t$IGG10),median(t$IGG20),median(t$IGG90))
  outtable    <- cbind(totalread,totalnumvar,geq10numvar,maxcountvar,mincountvar,avgcountvar,medcountvar)
  row.names(outtable) <- c('Input','IGG10','IGG20','IGG90')
  print (paste('Input table:',main),sep=' ',quote=F)
  print (outtable)
  }

fitness_stats <- function(t,main){
  F1  <- t[which(t[,1]==1),2]
  F2  <- t[which(t[,1]==2),2]
  F3  <- t[which(t[,1]==3),2]
  F4  <- t[which(t[,1]==4),2]
  HD1 <- c(mean(F1),median(F1),min(F1),max(F1))
  HD2 <- c(mean(F2),median(F2),min(F2),max(F2))
  HD3 <- c(mean(F3),median(F3),min(F3),max(F3))
  HD4 <- c(mean(F4),median(F4),min(F4),max(F4))
  outtable <- cbind(HD1,HD2,HD3,HD4)
  row.names(outtable) <- c('mean','median','min','max')
  print (paste('Condition:',main,sep=' '),quote=F)
  print (outtable)
  }

coverage_stats <- function(t){
  HD1 <- filter_t[which(filter_t$HD == 1),]
  HD2 <- filter_t[which(filter_t$HD == 2),]
  HD3 <- filter_t[which(filter_t$HD == 3),]
  HD4 <- filter_t[which(filter_t$HD == 4),]
  print (paste('HD1',length(HD1$Input),length(HD1$Input[which(HD1$Input>=10)]),sep=' '),quote=F)
  print (paste('HD2',length(HD2$Input),length(HD2$Input[which(HD2$Input>=10)]),sep=' '),quote=F)
  print (paste('HD3',length(HD3$Input),length(HD3$Input[which(HD3$Input>=10)]),sep=' '),quote=F)
  print (paste('HD4',length(HD4$Input),length(HD4$Input[which(HD4$Input>=10)]),sep=' '),quote=F)
  }

t <- read.table('result/Mutfit',header=1)
filter_t <- t[which(as.character(str_sub(t$mut,1,1)) != '_'
                  & as.character(str_sub(t$mut,2,2)) != '_'
                  & as.character(str_sub(t$mut,3,3)) != '_'
                  & as.character(str_sub(t$mut,4,4)) != '_'),]

basic_stats(t,'Unfilter') #Analyze the raw data
basic_stats(filter_t,'Filtered all stops') #Analyze those data filtered with stop
coverage_stats(filter_t) #Analyze input coverage for different HD

#Continue to fitness analysis for mutants with different HD
geq10t <- filter_t[which(as.character(filter_t$I10fit) != 'NA'),]
IGG10  <- cbind(geq10t$HD,geq10t$I10fit)
IGG20  <- cbind(geq10t$HD,geq10t$I20fit)
IGG90  <- cbind(geq10t$HD,geq10t$I90fit)
fitness_stats(IGG10,'IGG10')
fitness_stats(IGG20,'IGG20')
fitness_stats(IGG90,'IGG90')
