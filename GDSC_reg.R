#==================================================================
# This script is to do GDSC data analysis.
# Please run script "GDSC_pre.R" to get the preprocessed data first.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 07-Feb-2019
#==================================================================

library(doParallel)
library(glmnet)
library(tgp)
library(penalizedSVM)
# Please load package IPFStructPenalty after penalizedSVM
library(IPFStructPenalty)

load("GDSC_complete.RData")
# delete the cell line with extreme log(IC50)=-36.49 for drug "AP-24534"
GDSC$y <- GDSC$y[-116,]
GDSC$x <- GDSC$x[-116,]
n <- dim(GDSC$y)[1]
m <- dim(GDSC$y)[2]
p <-  GDSC$p
num.nonpen <- GDSC$num.nonpen
# run 10 repetitions
MultSeed <- c(2067, 6253, 8009, 10893, 32940, 43024, 59327, 72, 82702, 148392)
for(s in 1:10){
  cat("Seed ", s, ": ", MultSeed[s], "\n")
  set.seed(MultSeed[s])
  
  # select train dataset of 80% cell lines of each tissue
  trainID <- NULL
  for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  # make sure at least one mutated cell lines of each cancer tissue for every mutation features
  repeat{
    if(min(colSums(GDSC$x[trainID, num.nonpen+(1+p[1]+p[2]):sum(p)]))>=1) break
    trainID <- NULL
    for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  }
  
  x_test <- GDSC$x[-trainID,]
  x <- GDSC$x[trainID,]
  y_test <- GDSC$y[-trainID,]
  y <- GDSC$y[trainID,]
  
  x_test[,num.nonpen+1:p[1]] <- log(GDSC$x[-trainID, num.nonpen+1:p[1]])
  x[,num.nonpen+1:p[1]] <- log(GDSC$x[trainID, num.nonpen+1:p[1]])
  
  foldid <- sample(rep(seq(5),length=dim(x)[1]))
  methods <- c("lasso", "elastic-net", "IPF-lasso", "sIPF-elastic-net", "IPF-elastic-net", "tree-lasso", "IPF-tree-lasso")
  
  i<-1 # specify the method
  cat("Method ", i, ": ", methods[i], "\n")
  fit<-IPFStructPenaltyReg(x, y, x_test, y_test, p, foldid, N=10, min.ite=2,method=methods[i], num.nonpen=num.nonpen)
  save(fit, file=paste("GDSC",i,"_",MultSeed[s],".RData",sep=""))
} 
