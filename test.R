#================================================================================================================
# This script illustrates the use of the "IPFStructPenalty” R package with an example analysis of simulated data.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 07-Feb-2019
# depends on: sim.data.R
#================================================================================================================

# Build the package
# devtools::document("IPFStructPenalty")
# devtools::build("IPFStructPenalty")

## Install the package
#install.packages("IPFStructPenalty_0.1.1.tar.gz",repos = NULL,type = "source")
library(devtools)
install_github("zhizuio/IPFStructPenalty/IPFStructPenalty")


library(IPFStructPenalty)

source("sim.dat.R")
set.seed(1234)
n <- 100
m <- 24
p <- c(150, 150)
# simulate learning dataset
sim <- sim1(p=c(150, 150), n=n, m=m)
y <- sim$Y
x <- sim$X

# simulate validation dataset
sim_test<-sim1(p=p, n=n, m=m)
y_test <- sim_test$Y
x_test <- sim_test$X

# generate cross-validation ID
foldid <- sample(rep(seq(5),length=dim(sim$X)[1])) 

methods <- c("lasso", "elastic-net", "IPF-lasso", "sIPF-elastic-net", "tree-lasso", "IPF-tree-lasso")
# run the penalized regression
fit <- IPFStructPenaltyReg(x, y, x_test, y_test, p, foldid, method=methods[1])
