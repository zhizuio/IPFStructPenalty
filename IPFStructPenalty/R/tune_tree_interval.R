tune.tree.interval<-function(parms, x, y,
                                lambda = NULL, 
                                nfolds = 5,
                                num.nonpen = 0,
                                seed=12345, 
                                foldid=NULL,
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
  if(is.null(lambda)) cat("No given lambda sequence!")
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
      clusterEvalQ(cl, library(iterators))
      clusterEvalQ(cl, library(foreach))
      clusterEvalQ(cl, library(doParallel))
      clusterEvalQ(cl, library(IPFStructPenalty))
      registerDoParallel(cl)
      cvm0[1:cores] <- foreach(i = 1:cores, .combine=c, .packages= c('base','Matrix','MASS')) %dopar%{
        cv5(la.xx[i,2], la.xx[i,1])
      }
      stopCluster(cl)
      cvm0 <- rowMeans(matrix(cvm0, ncol=max(foldid)))
    }else{
      cvm0 <- numeric(length(lambda))
      cl <- makeCluster(min(length(lambda),15))
      clusterEvalQ(cl, library(iterators))
      clusterEvalQ(cl, library(foreach))
      clusterEvalQ(cl, library(doParallel))
      clusterEvalQ(cl, library(IPFStructPenalty))
      registerDoParallel(cl)
      cvm0[1:min(length(lambda),15)] <- foreach(i = 1:min(length(lambda),15), .combine=c, .packages= c('base','Matrix','MASS')) %dopar%{
        la.seq(lambda[i])
      }
      stopCluster(cl)
      
      if(length(lambda)>15){
        cvm0[(15+1):length(lambda)]<-apply(matrix(lambda[(15+1):length(lambda)],ncol=1), 1, la.seq)
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



