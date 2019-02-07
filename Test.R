#================================================================================================================
# This script illustrates the use of the "IPFStructPenalty” R package with an example analysis of simulated data.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 07-Feb-2019
# depends on: sim.data.R
#================================================================================================================

# load libraries
library(doParallel)
library(glmnet)
library(tgp)
library(penalizedSVM)
# Please load package IPFStructPenalty after penalizedSVM
library(IPFStructPenalty)

source("sim.dat.R")
set.seed(1234)
n <- 100
m <- 24
p <- c(150, 150)
sim <- sim.dat(p=c(150, 150), n=n, m=m)
foldid <- sample(rep(seq(5),length=dim(sim$X)[1])) 
p<-sim$p
y<-scale(sim$Y)[,]
x<-scale(sim$X)[,]

# generate external testing datasets
sim_test<-sim.dat(p=p, n=n, m=m)
y_test<-scale(sim_test$Y)[,]
x_test<-scale(sim_test$X)[,]

methods <- c("lasso", "elastic-net", "IPF-lasso", "sIPF-elastic-net", "IPF-elastic-net", "tree-lasso", "IPF-tree-lasso")
fit<-IPFStructPenaltyReg(x, y, x_test, y_test, p, foldid, method=methods[1])
