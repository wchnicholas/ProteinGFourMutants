#R code

coranalysis <- function(filename){
  print (filename)
  t <- read.table(filename,header=1)
  cort <- cor(t[,2:18])
  print (cort)
  print ("\n")
  }

coranalysis('result/HD4EpiIGG10')
coranalysis('result/HD4EpiIGG20')
coranalysis('result/HD4EpiIGG90')
