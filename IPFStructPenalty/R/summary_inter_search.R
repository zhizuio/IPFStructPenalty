#' IPFStructPenalty
#' @title Summary results
#' @description
#' Summarize the results from the oject of \code{epsgo}. 
#' @return An object with list "\code{summary.intsearch}"
#' \item{info}{the visited points of alpha, lambda, ipf and their resulting deviance and selected number of features}
#' \item{opt.alpha}{the optimized alpha}
#' \item{opt.lambda}{the optimized (first) lambda}
#' \item{opt.ipf}{the optimized penalty factors}
#' \item{opt.model}{the optimized model}
#' \itemize{alpha }{ the optimzed alpha}
#' \itemize{lambda}{ the optimzed (first) lambda}
#' \itemize{ipf   }{ the optimzed penalty factors}
#' \itemize{p     }{ a vector of the numbers of features from multiple data sources}
#' \itemize{nfolds}{ number of folds used for the cross-validation}
#' \itemize{cvreg }{ the cross-validation results}
#' @export
###########################################################################################################
summary.intsearch<-function(object,digits = max(3, getOption("digits") - 3), verbose=TRUE, first.n=5, ...){
  fit <- object
  lambdas <- unlist(sapply(sapply(fit$model, "[", "model"), "[", "lambda"))
  if(!identical(grep("tree",names(fit$model.list[[1]]$model$ipf)[1]), integer(0))){
    #fit <- object
    #lambdas <- unlist(sapply(sapply(fit$model, "[", "model"), "[", "lambda"))
    ipf <- data.matrix(fit$Xtrain)
    deviances <- fit$Ytrain
    opt.lambda <- lambdas[which.min(deviances)]
    opt.ipf <- ipf[which.min(deviances),]
    opt.error <- fit$fmin
    out <- list(info=data.frame(lambda=lambdas,ipf=ipf,deviance=deviances),
                opt.lambda=opt.lambda, opt.ipf=opt.ipf, opt.error=opt.error)
    class(out) <- "sum.intsearch"
    if(verbose){
      cat("Summary interval search \n\n")
      cat(paste("show the first", first.n,"out of",nrow(out$info),"entries\n"))
      print(out$info[1:first.n,])
      cat("\n..............................")
      
      cat("\n\n Optimal parameters found are: \n\n")
      cat(paste("lambda = ",round(out$opt.lambda,digits),
                "\t",
                "ipf = ",paste(as.character(opt.ipf),collapse=":"),
                "deviance = ",round(out$opt.error,digits)))
    }
    invisible(out)
  }else{
    #fit <- object
    #lambdas <- unlist(sapply(sapply(fit$model, "[", "model"), "[", "lambda"))
    if(identical(grep("alpha",colnames(fit$Xtrain)), integer(0))){
      fit$Xtrain <- data.frame(alpha=rep(1,dim(fit$Xtrain)[1]), fit$Xtrain)
      colnames(fit$Xtrain)[1] <- "alpha"
    } 
    alphas <- data.matrix(fit$Xtrain[,grep("alpha",colnames(fit$Xtrain))])
    if(identical(grep("ipf",colnames(fit$Xtrain)), integer(0))){
      ipf <- NA
    }else{
      ipf <- data.matrix(fit$Xtrain[,grep("ipf",colnames(fit$Xtrain))])
    }
    
    deviances <- fit$Ytrain
    
    # optimal models   
    opt.models <- sapply(fit$model.list, "[", "model") [which.min(fit$Ytrain)]
    if(length(grep("alpha",colnames(fit$Xtrain)))==1){
      if(names(fit$model[[1]]$model$cvreg)[1]!="cvl"){
        tmp.models<-sapply(sapply(sapply(fit$model, "[", "model"), "[", "cvreg"), "[", "glmnet.fit")
        n.features<-mapply(function(List, lam) List$df[which(List$lambda %in% lam)], tmp.models, lambdas)
      }else{
        tmp.models<-sapply(sapply(sapply(fit$model, "[", "model"), "[", "cvreg"), "[", "fullfit")
        n.features<-mapply(function(List) sum(List@penalized[-1]!=0), tmp.models)
      }
      opt.alpha <- opt.models[[1]]$alpha
      opt.lambda <- opt.models[[1]]$lambda
      opt.ipf <- opt.models[[1]]$ipf
      opt.error <- fit$fmin
      opt.n.features <- n.features[which.min(deviances)]
      out <- list(info=data.frame(alpha=alphas,lambda=lambdas,ipf=ipf,deviance=deviances,n.features=n.features),
                  opt.alpha=opt.alpha, opt.lambda=opt.lambda, opt.ipf=opt.ipf, opt.error=opt.error,
                  st.resp=fit$model$st.resp, opt.models=opt.models)
    }else{
      opt.alpha <- alphas[which.min(deviances),]
      opt.lambda <- lambdas[which.min(deviances)]
      opt.ipf <- ipf[which.min(deviances),]
      opt.error <- fit$fmin
      out <- list(info=data.frame(lambda=lambdas,alpha=alphas,ipf=ipf,deviance=deviances),
                  opt.alpha=opt.alpha, opt.lambda=opt.lambda, opt.ipf=opt.ipf, opt.error=opt.error)
    }
    
    class(out) <- "sum.intsearch"
    
    if(verbose){
      cat("Summary interval search \n\n")
      cat(paste("show the first", first.n,"out of",nrow(out$info),"entries\n"))
      print(out$info[1:first.n,])
      cat("\n..............................")
      
      cat("\n\n Optimal parameters found are: \n\n")
      cat(paste("alpha = ",paste(as.character(out$opt.alpha),collapse=","),
                "lambda = ",round(out$opt.lambda,digits),
                "\t",
                "ipf = ",paste(as.character(opt.ipf),collapse=":"),
                "deviance = ",round(out$opt.error,digits)))
    }
    invisible(out)
  }
  
}
