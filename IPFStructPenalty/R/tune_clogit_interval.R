#' IPFStructPenalty
#' @title Wrapper function for conditional logistic lasso objects.
#' @description
#' Wrapper function for conditional logistic lasso objects used by epsgo function. This function is mainly used within the function \code{epsgo}. See the \code{R} package \pkg{c060} for details.
#' @export
tune.clogit.interval<-function(parms, x=x, y=y,
                                         x_test=NULL,
                                         y_test=NULL,
                                         strata=NULL,
                                         lambda = NULL,
                                         nfolds = 3,
                                         foldid =NULL,
                                         seed=12345, 
                                         standardize.response=FALSE,
                                         p=NULL,
                                         parallel=FALSE,
                                         verbose=FALSE,
                                         search.path=FALSE,
                                         ...){
  
  # 1. decode the parameters ############################################################
  set.seed(seed)
  X <- data.matrix(x)
  y <- data.matrix(y)
  if(is.null(p) | length(p)==1) stop("The argument p must be a vector!")
  if(is.null(lambda)) stop("No given lambda sequence!")
  lambda <- matrix(lambda,ncol=1)
  if(!is.null(foldid)){
    nfolds <- max(foldid)
  }else{
    foldid <- rep(sample(rep(seq(nfolds),length(table(strata))/nfolds)), each=2)
  }
  
  if(!identical(grep("alpha", names(parms)), integer(0))){
    alpha<-round(parms[grep("alpha", names(parms))],dig=2)
    if(verbose) print(paste("alpha=",paste(as.character(alpha),collapse=",")))
  }else{
    alpha<-0
  }
  if(!identical(grep("ipf", names(parms)), integer(0))){
    ipf<-round(parms[grep("ipf", names(parms))],dig=2)
    if(verbose) print(paste("IPF=",paste(as.character(ipf),collapse=":")))
    if(length(p) != (length(ipf)+1)) stop("The arguments p and number of penalty factor searching ranges have to be the same length!")
  }else{
    stop("Please provide a searching range for the penalty factor!")
  }
  
  
  #  2. find optimal lambda for given alpha (lambda2) and penalty factors ###########################
  # find optimal lambda1 for given alpha (lambda2)
  
  CV_k <- function(la){
    
    lambda1 <- rep(la, p[1])
    for(i in 2:length(p)) lambda1 <- c(lambda1, rep(la,p[i])*ipf[i-1])
    resp <- Surv(y[,1], event= y[,2])
    fit.loglik <- cvl(resp~strata, X, lambda1=lambda1, lambda2=alpha, model="cox", fold=foldid)$fullfit@loglik
    return(fit.loglik)
    
  }
  
  if(sum(parallel)){
      cvm <- numeric(length(lambda))
      cl <- makeCluster(min(length(lambda),15))
      clusterEvalQ(cl, library(iterators))
      clusterEvalQ(cl, library(foreach))
      clusterEvalQ(cl, library(doParallel))
      clusterEvalQ(cl, library(penalized))
      registerDoParallel(cl)
      cvm[1:min(length(lambda),15)] <- foreach(i = 1:min(length(lambda),15), .combine=c, .packages= c('base','Matrix','MASS')) %dopar%{
        CV_k(lambda[i])
      }
      stopCluster(cl)
      
      if(length(lambda)>15){
        cvm[(15+1):length(lambda)]<-apply(matrix(lambda[(15+1):length(lambda)],ncol=1), 1, CV_k)
      }
    }else{
    cvm<-apply(lambda, 1, CV_k)
  }


  opt.lambda<-lambda[which.max(cvm)]
  # q.val= mean cross-validated error over the folds
  q.val<-max(cvm)
  
  lambda1 <- rep(lambda, p[1])
  for(i in 2:length(p)) lambda1 <- c(lambda1, rep(opt.lambda,p[i])*ipf[i-1])
  resp <- Surv(y[,1], event= y[,2])
  cv <- cvl(resp~strata, X, lambda1=opt.lambda, lambda2=alpha, model="cox", fold=foldid)
  
  
  if(!search.path){
    ret<-list(q.val=q.val, model=list(lambda=opt.lambda, alpha=alpha, ipf=ipf, p=p, nfolds=nfolds, cvreg=cv))
  }else{
    ret<-list(q.val=q.val, model=list(lambda=opt.lambda, alpha=alpha, ipf=ipf, p=p, nfolds=nfolds, cvreg=cv, search.cvm=c(lambda,cvm)))
  }
  
  return(ret)
}



