#' IPFStructPenalty
#' @title Structured penalized regression
#' @description
#' Function producing results of the structured penalized regression
#' @param x,y \code{x} is the input design matrix; \code{y} is the input response matrix
#' @param x_test,y_test \code{x} is the input validated design matrix; \code{y} is the input validated response matrix
#' @param p the number of predictors from different data source.
#' @param foldid an vector of values for the cross-validation.
#' @param num.nonpen number of predictors forced to be estimated (i.e., nonpenalization).
#' @param method specify the the method to optimize its penalty parameters. The penalty parameters of \code{elastic-net}, \code{IPF-lasso}, \code{sIPF-elastic-net}, \code{IPF-elastic-net}, \code{IPF-tree-lasso} and \code{clogitLasso} are optimzed by the EPSGO algorithm. The penalty parameter of \code{lasso} and \code{tree-lasso} is optimzed by cross-validation. The default method is \code{IPF-lasso} 
#' @param lambda optional user-supplied \code{lambda} sequence; default is NULL, and \code{espsgo} chooses its own sequence except the tree-lasso methods.
#' @param bounds bounds for the interval-searching parameters
#' @param strata.surv stratification variable for the Cox survival model.
#' @param search.path save the visited points, default is \code{FALSE}.
#' @param EI.eps he convergence threshold for the expected improvement between fmin and the updated point 
#' @param fminlower minimal value for the function Q.func, default is 0.
#' @param threshold threshold for estimated coefficients of the tree-lasso methods.
#' @param N define the number of start points depending on the dimensionality of the parameter space.
#' @param min.iter the minimus iterations after the initial \code{N} iterations.
#' @param seed random seed.
#' @param parallel If \code{TRUE}, use parallel foreach to fit each fold except parallelizing each lambda for the tree-lasso methods. If \code{c(TRUE,TRUE)}, use parallel foreach to fit each fold and each lambda. 
#' @param verbose print the middle search information, default is \code{TRUE}.
#' @return An object of list "\code{IPFStructPenaltyReg}" is returned:
#'  \item{cvm}{the mean cross-validated error}  
#'  \item{cvm_cv}{the mean cross-validated error if providing external dataset "\code{x_test}" and "\code{y_test}". } 
#'  \item{alpha}{optimized \code{alpha}}  
#'  \item{lambda}{optimized \code{lambda}}  
#'  \item{pred}{the prediction of the responses}  
#'  \item{ipf}{optimzed penalty factors}  
#'  \item{Beta}{estimate of the coefficients}  
#'  \item{cv}{number of nonzero coefficients}  
#' 
#' @references Zhao, Z. & Zucknick, M. (2019). \emph{Stuctured penalized regression for drug sensitivity prediction.} arXiv: 1905.00095.
#' @export
IPFStructPenaltyReg <- function(x, y, x_test=NULL, y_test=NULL, p, foldid, num.nonpen=0, method="IPF-lasso", 
                        lambda=NULL, bounds=NULL, strata.surv=NULL, search.path=FALSE, EI.eps=0.01, fminlower = 0,
                        threshold=0, N=NULL, min.iter=20, seed=1234,parallel=FALSE, verbose=TRUE,...){
  if((method!="lasso") & (method!="tree-lasso") & is.null(bounds)){
    if(method=="elastic-net"){
      bounds <- t(data.frame(alpha=c(0,1)))
    }else{
      # the default setting for bounds is only for two data sources
      if(method=="IPF-lasso" | method=="clogitLasso"){
        bounds <- t(data.frame(ipf=c(0.1,10)))
      }else{
        if(method=="sIPF-elastic-net"){
          bounds <- t(data.frame(alpha=c(0,1), ipf1 = c(0.1, 10)))
        }else{
          if(method=="IPF-elastic-net"){
            bounds <- t(data.frame(alpha1=c(0,1), alpha2=c(0,1), ipf1 = c(0.5, 8)))
          }else{
            if(method=="IPF-tree-lasso"){
              bounds <- t(data.frame(tree.ipf1 = c(0.1, 5)))
            }else{
              cat("Ooops! Please give searching bounds for EPSGO algorithm!")
            }
          }
        } 
      }
    } 
    colnames(bounds) <- c("lower", "upper")
  }
  
  x <- data.matrix(x)
  y <- data.matrix(y)
  if(is.null(x_test) & is.null(y_test)){
    x_test <- x
    y_test <- y
  }else{
    x_test <- data.matrix(x_test)
    y_test <- data.matrix(y_test) 
  }   
  m <- dim(y)[2]
  if(method!="clogitLasso") foldid2 <- rep(foldid,m)
  y_star <- as.vector(y)
  x_star <- kronecker(Diagonal(m), cbind(rep(1,dim(x)[1]),x))
  Beta <- as.vector(rbind(matrix(0,nrow=1+num.nonpen,ncol=m), matrix(10,ncol=dim(y)[2],nrow=sum(p,na.rm=T))))
  adpen <- rep(1, dim(x_star)[2])
  adpen[Beta==0] <- 0
  
  if(method=="lasso"){
    lambda <- glmnet(x_star, y_star, family="gaussian", nlambda=20, intercept=F, standardize.response=F, penalty.factor=adpen)$lambda
    la <- cv.glmnet(x_star, y_star, family="gaussian", lambda=lambda, intercept=F, foldid=foldid2, standardize.response=F, penalty.factor=adpen, parallel=parallel)
    mse_cv <- min(la$cvm)
    alpha <- 1
    lambda <- la$lambda.min
    
    # prediction
    fit <- glmnet(x_star, y_star, family="gaussian", intercept=F, lambda=seq(lambda*0.8,lambda*1.2,length=11), standardize.response=F, penalty.factor=adpen)
    Beta <- matrix(fit$beta[,ceiling(11/2)], nrow=dim(x)[2]+1, ncol=dim(y)[2])
    ypred <- cbind(rep(1,dim(x_test)[1]),x_test)%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
    ipf <- 1
  }
  
  #=============
  # Elastic-net
  #=============
  if(method=="elastic-net"){
    colnames(bounds) <- c("lower", "upper")
    fit <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds, lambda=lambda, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, num.nonpen=num.nonpen,
                 family = "gaussian", foldid = foldid, type.min = "lambda.min", p=p, intercept=T, standardize.response=F,
                 type.measure = "mse", min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, ...)
    sumint <- summary(fit, verbose=F)
    mse_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    lambda <- sumint$opt.lambda
    
    # prediction
    adpen <- adpen*(dim(x_star)[2])/sum(adpen)
    fit <- glmnet(x_star, y_star, family="gaussian", intercept=F, alpha=alpha, lambda=seq(lambda*0.8,lambda*1.2,length=11), standardize.response=F, penalty.factor=adpen)
    Beta <- matrix(fit$beta[,ceiling(11/2)], nrow=dim(x)[2]+1, ncol=dim(y)[2])
    ypred <- data.matrix(cbind(rep(1,dim(x_test)[1]),x_test))%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
    ipf <- 1
  }
  
  #=============
  # IPF-lasso, sIPFEN
  #=============
  if((method=="IPF-lasso") | (method=="sIPF-elastic-net")){
    fit <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds, lambda=lambda, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, num.nonpen=num.nonpen,
                 family = "gaussian", foldid = foldid, type.min = "lambda.min", p=p, intercept=T, standardize.response=F,
                 type.measure = "mse", min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, ...)
    sumint <- summary(fit, verbose=F)
    mse_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    lambda <- sumint$opt.lambda
    ipf <- sumint$opt.ipf
    
    beta<-matrix(0, nrow=1+num.nonpen+sum(p), ncol=m)
    for(s in 1:length(p)) beta[1+num.nonpen+ifelse(s==1,0,sum(p[1:(s-1)]))+1:p[s],]<-10+10*s
    adpen <- rep(1, dim(x_star)[2])
    adpen[as.vector(beta)==0] <- 0
    adpen[as.vector(beta)==20] <- 1
    for(i in 1:length(ipf)) adpen[as.vector(beta)==(10+10*(i+1))] <- ipf[i]
    adpen <- adpen*(dim(x_star)[2])/sum(adpen)
    
    # prediction
    Beta <- matrix(glmnet(x=x_star,y=y_star, family="gaussian",
                          alpha=alpha,
                          offset = NULL,
                          lambda = seq(lambda*0.8,lambda*1.2,length=11),
                          penalty.factor=adpen,
                          intercept=F,
                          standardize.response=F)$beta[,ceiling(11/2)], nrow=dim(x)[2]+1, ncol=dim(y)[2])
    ypred <- data.matrix(cbind(rep(1,dim(x_test)[1]),x_test))%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
  }
  
  #=============
  # IPFEN with varying \alpha
  #=============
  if(method=="IPF-elastic-net"){
    fit <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds, lambda=lambda,  N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, num.nonpen=num.nonpen,
                 family = "gaussian", foldid = foldid, type.min = "lambda.min", p=p, intercept=T, standardize.response=F,
                 type.measure = "mse", min.iter=min.iter, verbose=verbose, parallel=parallel, search.path=T, EI.eps=EI.eps, ...)
    sumint <- summary(fit, verbose=F)
    mse_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    alpha <- alpha
    lambda <- sumint$opt.lambda# * sum(dim(x))*m*2
    ipf <- sumint$opt.ipf
    
    beta<-matrix(0, nrow=1+num.nonpen+sum(p), ncol=m)
    for(s in 1:length(p)) beta[1+num.nonpen+ifelse(s==1,0,sum(p[1:(s-1)]))+1:p[s],]<-10+10*s
    y2 <- as.vector(rbind(matrix(y_star, ncol=m), matrix(rep(0, m*sum(p)),ncol=m)))
    sd.y2 <- sd(y2)*sqrt((length(y2)-1)/length(y2))
    y2 <- y2/sd.y2
    adpen00 <- adpen
    adpen[as.vector(beta)==20] <- 1
    for(i in 1:length(ipf)) adpen[as.vector(beta)==(10+10*(i+1))] <- ipf[i]
    adpen0 <- adpen
    adpen <- adpen*(dim(x_star)[2])/sum(adpen)
    
    adpen00 <- adpen00*(dim(x_star)[2])/sum(adpen00)
    lambda2 <- lambda*c(1,ipf)*m*dim(x)[1]
    pf.r <- c(1, ipf) # ratios of penalty factors
    x2_test <- x2 <- aug <- NULL
    # generate the transformed X
    for(s in 1:length(p)){
      if(s==1){
        x2 <- pf.r[s]/alpha[s]*x[,num.nonpen+(ifelse(s==1,1,cumsum(p)[s-1]+1)):cumsum(p)[s]]
        x2_test <- pf.r[s]/alpha[s]*x_test[,num.nonpen+(ifelse(s==1,1,cumsum(p)[s-1]+1)):cumsum(p)[s]]
        aug <- pf.r[s]/alpha[s]*sqrt(1/2*lambda2[s]*(1-alpha[s]))*diag(p[s])
      } 
      if(s>1){
        x2 <- cbind(x2, pf.r[s]/alpha[s]*x[,num.nonpen+(ifelse(s==1,1,cumsum(p)[s-1]+1)):cumsum(p)[s]])
        x2_test <- cbind(x2_test, pf.r[s]/alpha[s]*x_test[,num.nonpen+(ifelse(s==1,1,cumsum(p)[s-1]+1)):cumsum(p)[s]])
        aug <- bdiag(aug, pf.r[s]/alpha[s]*sqrt(1/2*lambda2[s]*(1-alpha[s]))*diag(p[s]))
      } 
    }
    org.const <- c(rep(1,dim(x2)[1]),rep(0,sum(p)))
    if(num.nonpen>0){
      x2 <- cbind(x[,c(1:num.nonpen)], x2)
      x2_test <- cbind(x_test[,c(1:num.nonpen)], x2_test)
      aug <- cbind(matrix(0,ncol=num.nonpen,nrow=sum(p)), aug)
    }
    x2_star <- kronecker(Diagonal(m), cbind(org.const,rbind(x2, aug))/sd.y2)
    lam_star <- lambda2[1]/(sd.y2^2*sum(dim(x))*m*2)
    Beta <- matrix(glmnet(x2_star,y2,family="gaussian",lambda=seq(lam_star*0.8,lam_star*1.2,length=11),intercept=F,penalty.factor=adpen00)$beta[,ceiling(11/2)], ncol=dim(y)[2])
    for(s in 1:length(p)) Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,] <- 1/alpha[s]*pf.r[s]*Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,]
    ypred <- data.matrix(cbind(rep(1,dim(x_test)[1]),x_test)) %*% Beta
    
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
  }
  
  #==================
  # Tree-lasso
  #==================
  
  if(method=="tree-lasso"){
    if(is.null(lambda)){
      cat("Warning: Please provide a proper lambda sequence!")
      lambda <- seq(1,2.5,length=10)
      # fun.lambda <- function(x,y){
      #   sum(matrix(x,nrow=1)%*%(y-matrix(rep(apply(y,2,mean),each=dim(y)[1]),ncol=dim(y)[2])))
      # }
      # lambda.max <- sqrt(max(sapply(split(x, rep(1:ncol(x), each=nrow(x))), fun.lambda, y=y)))
      # lambda <- seq(lambda.max*0.1, lambda.max, length=5)
    }
    cvm0 <- numeric(length(lambda))
    tree.parm0 <- tree.parms(y=y)
    cv5<-function(xx,la) {sum((y[foldid==xx,] - cbind(rep(1,sum(foldid==xx)),x[foldid==xx,]) %*% tree.lasso(x=x[!foldid==xx,], y=y[!foldid==xx,],lambda=la,tree.parm=tree.parm0,num.nonpen=num.nonpen, threshold=threshold))^2)/(sum(foldid==xx)*(dim(y)[2]))}
    la.seq<-function(la) {mean(sapply(1:max(foldid), cv5, la))}
    la.xx <- cbind(rep(lambda,times=max(foldid)), rep(1:max(foldid),each=length(lambda)))
    if(parallel){
      cores <- length(lambda) * max(foldid)
      cvm0 <- numeric(cores)
      cl <- makeCluster(cores)
      clusterEvalQ(cl, library(iterators))
      clusterEvalQ(cl, library(foreach))
      clusterEvalQ(cl, library(doParallel))
      clusterEvalQ(cl, library(IPFStructPenalty))
      registerDoParallel(cl)
      cvm0[1:cores] <- foreach(i = 1:cores, .combine=c, .packages= c('base','foreach','Matrix','MASS','IPFStructPenalty')) %dopar%{
        cv5(la.xx[i,2], la.xx[i,1])
      }
      stopCluster(cl)
    }else{
      lambda <- matrix(lambda,ncol=1)
      cvm0<-apply(lambda, 1, la.seq)
    }
    
    cvm0 <- rowMeans(matrix(cvm0,ncol=max(foldid)))
    mse_cv <- min(cvm0)
    alpha <- 1
    lambda <- lambda[which.min(cvm0)]
    
    # prediction
    Beta <- tree.lasso(x=x, y=y, lambda=lambda, tree.parm=tree.parm0, num.nonpen=num.nonpen)
    ypred <- data.matrix(cbind(rep(1,dim(x_test)[1]),x_test))%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
    ipf <- 1
  }
  
  #==================
  # Tree-IPF-lasso
  #==================
  
  if(method=="IPF-tree-lasso"){
    if(is.null(lambda)){
      cat("Warning: Please provide a proper lambda sequence!")
      lambda <- seq(1,3,length=10)
      # fun.lambda <- function(x,y){
      #   sum(matrix(x,nrow=1)%*%(y-matrix(rep(apply(y,2,mean),each=dim(y)[1]),ncol=dim(y)[2])))
      # }
      # lambda.max <- sqrt(max(sapply(split(x, rep(1:ncol(x), each=nrow(x))), fun.lambda, y=y)))
      # lambda <- seq(lambda.max*0.1, lambda.max, length=5)
    }
    tree.parm0 <- tree.parms(y=y)
    fit <- epsgo(Q.func = "tune.tree.interval", bounds = bounds, lambda=lambda, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, 
                 intercept=TRUE, foldid = foldid, p=p, standardize.response=F, num.nonpen=num.nonpen,
                 min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, threshold=threshold, ...)
    sumint <- summary(fit, verbose=F)
    
    ipf <- sumint$opt.ipf
    mse_cv <- sumint$opt.error
    alpha <- 1
    lambda <- sumint$opt.lambda
    
    # prediction
    Xtemp <- x
    for(i in 1:(length(p)-1)) Xtemp[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- Xtemp[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]]/ipf[i]
    Beta <- tree.lasso(x=Xtemp, y=y, lambda=lambda, tree.parm=tree.parm0, num.nonpen=num.nonpen, threshold=threshold)
    pf.r <- c(1, ipf)
    for(s in 1:length(p)) Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,] <- 1/pf.r[s]*Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,]  
    ypred <- data.matrix(cbind(rep(1,dim(x_test)[1]),x_test))%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
  }
  
  #==================
  # conditional logistic-lasso
  #==================
  if(method=="clogitLasso"){
    fit <- epsgo(Q.func = "tune.clogit.interval", strata.surv=strata.surv, bounds = bounds, lambda=lambda, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y,
                 foldid = foldid, p=p, min.iter=min.iter, 
                 verbose=verbose, parallel=parallel, EI.eps=EI.eps, ...)
    sumint <- summary(fit, verbose=F)
    mse_cv <- sumint$opt.error
    lambda <- sumint$opt.lambda
    ipf <- sumint$opt.ipf
    
    # beta<-matrix(0, nrow=1+num.nonpen+sum(p), ncol=m)
    # for(s in 1:length(p)) beta[1+num.nonpen+ifelse(s==1,0,sum(p[1:(s-1)]))+1:p[s],]<-10+10*s
    # adpen <- rep(1, dim(x_star)[2])
    # adpen[as.vector(beta)==0] <- 0
    # adpen[as.vector(beta)==20] <- 1
    # for(i in 1:length(ipf)) adpen[as.vector(beta)==(10+10*(i+1))] <- ipf[i]
    # adpen <- adpen*(dim(x_star)[2])/sum(adpen)
    
    ## prediction
    # Beta <- clogitLasso(X=x_star,y=y_star, strata.surv=strata.surv,
    #                       fraction = lambda,
    #                       adpative = TRUE,
    #                       p.fact = adpen)$beta
    fit.clogit <- penalized(Surv(y[,1], event= y[,2])~strata.surv, X, lambda1=lambda, lambda2=alpha, model="cox")
    Beta <- matrix(fit.clogit@penalized, ncol=1)
    ypred <- fit.clogit@lin.pred
    mse_val <- NA
    vs <- sum(Beta!=0)
    alpha <- 1
  }
  
  if(method %in% c("lasso","tree-lasso")){
    return(list(cvm=mse_val, cvm_cv=mse_cv, alpha=alpha, lambda=lambda, pred=ypred, ipf=ipf, Beta=Beta, cv=vs))
  }else{
    return(list(cvm=mse_val, cvm_cv=mse_cv, alpha=alpha, lambda=lambda, pred=ypred, ipf=ipf, Beta=Beta, cv=vs, search.fit=fit))
  }
}
