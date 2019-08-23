#================================================================================================================
# This is the full script for the analysis of simulated data in the pre-print Zhao & Zucknick (2019).
# The analysis was originally run on cluster. It will be extremely slow to run this in the PC.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 23-August-2019
# depends on: sim.data.R
#================================================================================================================

library(IPFStructPenalty)
# load the simulation functions
source("sim.dat.R")
# create a folder to save all results
dir.create("results/")

# set two kinds of data sources
para<-matrix(c(150, 150, 500, 150),nrow=2,byrow=T)

# this is the function to do simulation analysis
SimAnalysis <- function(num=1, scenario=1, hete=1, method="IPF-lasso", lambda=NULL, bounds=NULL,
                        threshold=0, N=21, min.iter=10, seed=1234, parallel=FALSE, verbose=TRUE){
  set.seed(num)
  # generate learning dataset
  repeat{
    sim <- do.call(paste("sim",scenario,sep=""), list(p=para[hete,]))
    foldid <- sample(rep(seq(5),length=dim(sim$X)[1]))
    min.height <- mheight(sim$Y, foldid)
    if(min.height<0.7) break                
  }
  Beta<-sim$Beta
  p<-sim$p
  y<-scale(sim$Y)[,]
  x<-scale(sim$X)[,]
  
  # generate validation dataset
  repeat{
    sim_test <- do.call(paste("sim",scenario,sep=""), list(p=para[hete,]))
    min.height_test <- mheight(sim_test$Y)
    if(min.height_test<0.7) break                
  }
  y_test <- scale(sim_test$Y)[,]
  x_test <- scale(sim_test$X)[,]
  
  cvm <- array(0,dim=c(6000,24))
  methods <- c("lasso", "elastic-net", "IPF-lasso", "sIPF-elastic-net", "tree-lasso", "IPF-tree-lasso")
  fit <- IPFStructPenaltyReg(x, y, x_test, y_test, p, foldid, method=methods[method.idx], lambda=lambda, bounds=bounds, 
                               threshold=threshold, N=N, min.iter=min.iter, seed=seed, parallel=parallel, verbose=verbose, lib.loc=lib.loc)
  if(method.idx>1){
    if(file.exists(paste("results/SimScenario", scenario, "_", hete, "_", num, ".RData", sep=""))){
      load(paste("results/SimScenario", scenario, "_", hete, "_", num, ".RData", sep=""))
    }else{
      stop("Warning: Please run the lasso model (method.idx=1) first!")
    }
  }
  cvm[method.idx,1:5] <- c(fit$cvm, fit$cvm_cv, fit$lambda, fit$ipf,
                   ifelse(length(fit$alpha)==1,fit$alpha,round(fit$alpha[1],dig=2)*100+fit$alpha[2]))
  cvm[(11+100*(method.idx-1)):(10+100*method.idx),1:24] <- array(as.vector(fit$pred),dim=c(100,24))
  cvm[(1001+(sum(p)+1)*(method.idx-1)):(1000+(sum(p)+1)*method.idx),1:24] <- fit$Beta
  
  save(cvm,file=paste("results/SimScenario", scenario, "_", hete, "_", num, ".RData", sep=""))
}

# Doing 50 simulations for the 3 scenarios, 2 feature cases for all 6 methods
# NOTE: It's better to parallelize the following loops on cluster as many as possible to speed the running
# be careful to give proper arguments "lambda", "bounds", "N" and "min.iter" in some situations
for(num in 1:50)
  for(scenario in 1:3)
    for(hete in 1:2)
      for(method.idx in 1:6)
        start.time=Sys.time();SimAnalysis(num, scenario, hete, method=method.idx, lambda=NULL, bounds=NULL,parallel=T);(consum.time=Sys.time()-start.time)
      

#=========================
# load all results of the simulation analysis and show their prediction performance
#=========================
B <- 50
cvm1 <- cvm2 <- cvm3 <- cvm4 <- cvm5 <- cvm6 <- array(0,dim=c(B,6000,24))
# two different structures of B1 & B2, b=c(0.6,0.2), p=c(150,150)
filelist<-list.files(path="results/",pattern=paste("SimScenario1_1_*",sep=""))[1:B]
datalist<-lapply(paste("results/",filelist,sep=""), function(x){load(x);cvm})
for(i in 1:B) cvm1[i,,]<-datalist[[i]][,,2]
# two different structures of B1 & B2, b=c(0.6,0.2), p=c(500,150)
filelist<-list.files(path="results/",pattern=paste("SimScenario1_2_*",sep=""))[1:B]
datalist<-lapply(paste("results/",filelist,sep=""), function(x){load(x);cvm})
for(i in 1:B) cvm2[i,,]<-datalist[[i]][,,2]
# two different un-structures of B1 & B2, b=c(0.6,0.2), p=c(150,150)
filelist<-list.files(path="results/",pattern=paste("SimScenario2_1_*",sep=""))[1:B]
datalist<-lapply(paste("results/",filelist,sep=""), function(x){load(x);cvm})
for(i in 1:B) cvm3[i,,]<-datalist[[i]][,,2]
# two different un-structures of B1 & B2, b=c(0.6,0.2), p=c(150,150)
filelist<-list.files(path="results/",pattern=paste("SimScenario2_2_*",sep=""))[1:B]
datalist<-lapply(paste("results/",filelist,sep=""), function(x){load(x);cvm})
for(i in 1:B) cvm4[i,,]<-datalist[[i]][,,2]
# two similar structures of B1 & B2, b=c(0.6,0.2), p=c(150,150)
filelist<-list.files(path="results/",pattern=paste("SimScenario3_1_*",sep=""))[1:B]
datalist<-lapply(paste("results/",filelist,sep=""), function(x){load(x);cvm})
for(i in 1:B) cvm5[i,,]<-datalist[[i]][,,2]
# two similar structures of B1 & B2, b=c(0.6,0.2), p=c(500,150)
filelist<-list.files(path="results/",pattern=paste("SimScenario3_2_*",sep=""))[1:B]
datalist<-lapply(paste("results/",filelist,sep=""), function(x){load(x);cvm})
for(i in 1:B) cvm6[i,,]<-datalist[[i]][,,2]


lab<-c("lasso","elastic","IPF-lasso","sIPFEN","tree-lasso","IPF-tree-lasso")
x0<-rep(lab,each=B)
x1 <-factor(x0, levels=c("lasso","IPF-lasso","elastic","sIPFEN","tree-lasso","IPF-tree-lasso"))
# plot MSE from validation data
layout(matrix(1:8,nrow=2,byrow=F))
par(mar=c(1,4,4,1))
cvm0<-cvm1
y0<-c(cvm0[,1,1],cvm0[,2,1],cvm0[,3,1],cvm0[,4,1],cvm0[,5,1],cvm0[,6,1])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.09,0.22))
points(colMeans(cvm0[,c(1,3,2,4,5,6),1]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["val"]), line=2.2)
title("(a) Scenario 1")
mtext("p=(150,150)",cex=3/4,line=.3)
cvm0<-cvm2
y0<-c(cvm0[,1,1],cvm0[,2,1],cvm0[,3,1],cvm0[,4,1],cvm0[,5,1],cvm0[,6,1])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.09,0.22))
points(colMeans(cvm0[,c(1,3,2,4,5,6),1]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["val"]), line=2.2)
title("(b) Scenario 1")
mtext("p=(500,150)",cex=3/4,line=.3)
cvm0<-cvm3
y0<-c(cvm0[,1,1],cvm0[,2,1],cvm0[,3,1],cvm0[,4,1],cvm0[,5,1],cvm0[,6,1])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.05,0.18))
points(colMeans(cvm0[,c(1,3,2,4,5,6),1]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["val"]), line=2.2)
title("(c) Scenario 2")
mtext("p=(150,150)",cex=3/4,line=.3)
cvm0<-cvm4
y0<-c(cvm0[,1,1],cvm0[,2,1],cvm0[,3,1],cvm0[,4,1],cvm0[,5,1],cvm0[,6,1])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.05,0.18))
points(colMeans(cvm0[,c(1,3,2,4,5,6),1]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["val"]), line=2.2)
title("(d) Scenario 2")
mtext("p=(500,150)",cex=3/4,line=.3)
cvm0<-cvm5
y0<-c(cvm0[,1,1],cvm0[,2,1],cvm0[,3,1],cvm0[,4,1],cvm0[,5,1],cvm0[,6,1])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.3,0.43))
points(colMeans(cvm0[,c(1,3,2,4,5,6),1]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["val"]), line=2.2)
title("(e) Scenario 3")
mtext("p=(150,150)",cex=3/4,line=.3)
cvm0<-cvm6
y0<-c(cvm0[,1,1],cvm0[,2,1],cvm0[,3,1],cvm0[,4,1],cvm0[,5,1],cvm0[,6,1])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.3,0.43))
points(colMeans(cvm0[,c(1,3,2,4,5,6),1]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["val"]), line=2.2)
title("(f) Scenario 3")
mtext("p=(500,150)",cex=3/4,line=.3)

par(mar=c(5,1,4,0))
plot(1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend(0,10,
       legend = c("lasso","IPF-lasso","elastic-net","sIPF-elastic-net","tree-lasso","IPF-tree-lasso"), 
       fill = gray.colors(6), bty="n", y.intersp=2)
points(1, 1 , pch = 24, bg = "darkgrey", cex=1.5)
text(2.4, 1, "mean")

# plot MSE of cross-validation from learning data
layout(matrix(1:8,nrow=2,byrow=F))
par(mar=c(1,4,4,1))
cvm0<-cvm1
y0<-c(cvm0[,1,2],cvm0[,2,2],cvm0[,3,2],cvm0[,4,2],cvm0[,5,2],cvm0[,6,2])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.1,0.26))
points(colMeans(cvm0[,c(1,3,2,4,5,6),10]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["CV"]), line=2.2)
title("(a) Scenario 1")
mtext("p=(150,150)",cex=3/4,line=.3)
cvm0<-cvm2
y0<-c(cvm0[,1,2],cvm0[,2,2],cvm0[,3,2],cvm0[,4,2],cvm0[,5,2],cvm0[,6,2])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.1,0.26))
points(colMeans(cvm0[,c(1,3,2,4,5,6),10]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["CV"]), line=2.2)
title("(b) Scenario 1")
mtext("p=(500,150)",cex=3/4,line=.3)
cvm0<-cvm3
y0<-c(cvm0[,1,2],cvm0[,2,2],cvm0[,3,2],cvm0[,4,2],cvm0[,5,2],cvm0[,6,2])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.06,0.22))
points(colMeans(cvm0[,c(1,3,2,4,5,6),10]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["CV"]), line=2.2)
title("(c) Scenario 2")
mtext("p=(150,150)",cex=3/4,line=.3)
cvm0<-cvm4
y0<-c(cvm0[,1,2],cvm0[,2,2],cvm0[,3,2],cvm0[,4,2],cvm0[,5,2],cvm0[,6,2])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.06,0.22))
points(colMeans(cvm0[,c(1,3,2,4,5,6),10]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["CV"]), line=2.2)
title("(d) Scenario 2")
mtext("p=(500,150)",cex=3/4,line=.3)
cvm0<-cvm5
y0<-c(cvm0[,1,2],cvm0[,2,2],cvm0[,3,2],cvm0[,4,2],cvm0[,5,2],cvm0[,6,2])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.305,0.465))
points(colMeans(cvm0[,c(1,3,2,4,5,6),10]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["CV"]), line=2.2)
title("(e) Scenario 3")
mtext("p=(150,150)",cex=3/4,line=.3)
cvm0<-cvm6
y0<-c(cvm0[,1,2],cvm0[,2,2],cvm0[,3,2],cvm0[,4,2],cvm0[,5,2],cvm0[,6,2])
boxplot(y0~x1,boxwex=0.7,col=gray.colors(6),xaxt="n",ylim=c(0.305,0.465))
points(colMeans(cvm0[,c(1,3,2,4,5,6),10]), pch = 24, bg = "darkgrey")
title(ylab=expression("MSE"["CV"]), line=2.2)
title("(f) Scenario 3")
mtext("p=(500,150)",cex=3/4,line=.3)

par(mar=c(5,1,4,0))
plot(1:10, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend(0,10,
       legend = c("lasso","IPF-lasso","elastic-net","sIPF-elastic-net","tree-lasso","IPF-tree-lasso"), 
       fill = gray.colors(6), bty="n", y.intersp=2)
points(1, 1 , pch = 24, bg = "darkgrey", cex=1.5)
text(2.4, 1, "mean")

  