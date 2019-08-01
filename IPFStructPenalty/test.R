remove.packages("IPFStructPenalty")
devtools::document("/Users/zhiz/Downloads/IPFStructPenalty")
devtools::build("/Users/zhiz/Downloads/IPFStructPenalty")

## Install the package
install.packages("/Users/zhiz/Downloads/IPFStructPenalty_0.1.1.tar.gz",repos = NULL,type = "source")

# load libraries
#library(doParallel)
#library(penalizedSVM)
#library(penalized)
library(IPFStructPenalty)

# generate data
y <- rep(c(1,0), 100)
X <- cbind(matrix (rnorm(4000, 0, 1), ncol = 20), matrix (rnorm(16000, 2, 0.6), ncol = 80)) # pure noise
pairs<- sort(rep(1:100, 2))

# re-build a two columns' matrix which creates a survival object where cases have events and controls are censored, all at the same time
y <- cbind(rep(1,200), rep(c(1,0),100)) 

# give the lambda sequence for the first data source
lambda <- seq(0.8, 1.5, length=2)
# give the searching range of penalty factors 
bounds <- t(data.frame(ipf=c(0.1,10)))
colnames(bounds) <- c("lower", "upper")

# use EPSGO algorithm to search the optimal penalty factors based on the given lambda corresponding to the first data source
fit <- epsgo(Q.func = "tune.clogit.interval", x = X, y = y, strata = pairs, nfold = 5, bounds = bounds, lambda=lambda, p=c(20,80), 
             N=10, min.iter = 5, seed = 123, fminlower = -1000, verbose=TRUE, parallel=FALSE)
sumint <- summary(fit, verbose=F)
opt.lambda.base <- sumint$opt.lambda
opt.ipf <- c(1,sumint$opt.ipf)
opt.fit <- sumint$opt.models$model$cvreg$fullfit

