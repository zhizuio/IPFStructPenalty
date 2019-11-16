#' IPFStructPenalty
#' @title Wrapper function for tree-lasso objects.
#' @description
#' Wrapper function for tree-lasso objects used by epsgo function. This function is mainly used within the function \code{epsgo}.
#' 
#' @importFrom parallel detectCores makeCluster	stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' 
#' @param parms tuning parameter alpha for the tree-lasso object.
#' @param x,y \code{x} is a matrix where each row refers to a sample a each column refers to a gene; \code{y} is a factor which includes the class for each sample
#' @param lambda A user supplied lambda sequence. 
#' @param nfolds number of cross-validation's folds, default 5.
#' @param foldid an optional vector of values between 1 and nfold identifying what fold each observation is in. If supplied, nfold can be missing.
#' @param num.nonpen number of predictors forced to be estimated (i.e., nonpenalization).
#' @param seed random seed
#' @param intercept should  intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE}).
#' @param standardize.response standardization for the response variables. Default: \code{TRUE}.
#' @param p the number of predictors from different data source.
#' @param verbose print the middle search information, default is \code{TRUE}.
#' @param parallel  If \code{TRUE}, use parallel foreach to fit each lambda. If \code{c(TRUE,TRUE)}, use parallel foreach to fit each lambda and each fold. 
#' @param search.path save the visited points, default is \code{FALSE}.
#' @param threshold threshold for estimated coefficients of the tree-lasso models.
#' @export
tune.tree.interval<-function(parms, x, y,
                                lambda = NULL, 
                                nfolds = 5,
                                foldid=NULL,
                                num.nonpen = 0,
                                seed=12345, 
                                intercept=TRUE,
                                standardize.response=FALSE,
                                p=NULL,
                                verbose=FALSE,
                                parallel=FALSE,
                                search.path=FALSE,
                                threshold=threshold,
                                ...){
  
  # 1. decode the parameters ############################################################
  
  parms<-round(parms,digits=2)
  if (verbose) print(paste("IPF=",paste(as.character(parms),collapse=":")))

  if(standardize.response) y<-scale(y)[,]
  set.seed(seed)
  if(is.null(foldid)) foldid<-sample(rep(seq(nfolds),length=dim(x)[1]))
  
  #=========
  # using augmented data
  #=========
  if(is.null(lambda)) stop("No given lambda sequence!")
  lambda <- matrix(lambda,ncol=1)
  x2 <- x
  for(i in 1:(length(p)-1)) x2[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- x2[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]]/parms[i]
 
  tree.parm0 <- tree.parms(y=y)
  cv5 <- function(xx,la){
    sum((y[foldid==xx,] - cbind(rep(ifelse(intercept,1,NULL),sum(foldid==xx)),x2[foldid==xx,]) %*% tree.lasso(x=x2[!foldid==xx,], y=y[!foldid==xx,],lambda=la,tree.parm=tree.parm0,num.nonpen=num.nonpen,intercept=intercept,threshold=threshold))^2)/(prod(dim(y[foldid==xx,])))
  }
  la.seq <- function(la) {mean(sapply(1:max(foldid), cv5, la))}
  
  if(sum(parallel)){
    if(sum(parallel)>1){
      la.xx <- cbind(rep(lambda,times=max(foldid)), rep(1:max(foldid),each=length(lambda)))
      cores <- length(lambda)*max(foldid)
      cvm0 <- numeric(cores)
      cl <- makeCluster(cores)
      #clusterEvalQ(cl, library(IPFStructPenalty))
      registerDoParallel(cl)
      cvm0[1:cores] <- foreach(i = 1:cores, .combine=c, .packages= c('base','Matrix','MASS')) %dopar%{
        cv5(la.xx[i,2], la.xx[i,1])
      }
      stopCluster(cl)
      cvm0 <- rowMeans(matrix(cvm0, ncol=max(foldid)))
    }else{
      cvm0 <- numeric(length(lambda))
      nCores <- min(length(lambda),16,detectCores()-1)
      cl <- makeCluster(nCores)
      # clusterEvalQ(cl, library(mlegp,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(tgp,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(statmod,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(corpcor,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(e1071,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(lhs,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(iterators,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(foreach,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(doParallel,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(penalizedSVM,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(penalized,lib.loc="/home/zhiz/RLibrary"))
      # clusterEvalQ(cl, library(IPFStructPenalty,lib.loc="/home/zhiz/RLibrary"))
      #clusterEvalQ(cl, library(IPFStructPenalty,lib.loc=NULL))
      registerDoParallel(cl)
      cvm0[1:min(length(lambda),nCores)] <- foreach(i = 1:min(length(lambda),nCores), .combine=c, .packages= c('base','Matrix','MASS')) %dopar%{
        la.seq(lambda[i])
      }
      stopCluster(cl)
      
      if(length(lambda)>nCores){
        cvm0[(nCores+1):length(lambda)]<-apply(matrix(lambda[(nCores+1):length(lambda)],ncol=1), 1, la.seq)
      }
    }
    
  }
  if(sum(parallel)<1){
    cvm0<-apply(lambda, 1, la.seq)
  }
  #=========
  
  opt.lambda<-lambda[which.min(cvm0)]
  # q.val= mean cross-validated error over the folds
  q.val<-min(cvm0)
  
  if(!search.path){
    ret<-list(q.val=q.val, model=list(lambda=opt.lambda, ipf=parms, nfolds=nfolds))
  }else{
    ret<-list(q.val=q.val, model=list(lambda=opt.lambda, ipf=parms, nfolds=nfolds, search.cvm=c(lambda,cvm0)))
  }
  return(ret)
}



