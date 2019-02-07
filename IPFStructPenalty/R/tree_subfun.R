accgrad <- function(y, x, lambda, Tree, C, g_idx, TauNorm,  mu, option, num.nonpen=0, intercept=TRUE){
  
  #Y Centered Matrix: N by K
  #X Centered Matrix: N by p
  #lam: lambda
  #Tree: sparse matrix: group info. rows: number of group, cols: number of tasks
  #Tw: n_group by 1: weight for each group
  #C: note \sum_|g| by K
  #g_idx: n_group by 2, group index
  #L1, Lipschitz cond
  #TauNorm: \|\Tau\|_1,2^2 
  #mu: mu in nesterov paper
  #maxiter
  
  Y <- data.matrix(y)
  X <- data.matrix(x)
  C <- data.matrix(C)
  p <- dim(X)[2]
  #K <- dim(Tree)[2]
  
  #obj <- rep(0, option["maxiter"])
  obj <- 0
  
  
  XX <- crossprod(X)
  XY <- crossprod(X, Y)
  
  L <- eigen(XX)$values[1] + lambda^2 * TauNorm/mu
  C <- C * lambda
  if(intercept){
    bw0 <- matrix(0, ncol=dim(Y)[2])
  }else{
    b0 <- bw0 <- matrix(0, nrow=0, ncol=dim(Y)[2])
  }
  
  #We can also use "bw0 <- matrix(0, ncol=K)" since "dim(Y)[2]=K"
  if(num.nonpen==0){
    bw1 <- matrix(0, nrow=p, ncol=dim(Y)[2])
    bx <- rbind(bw0, bw1)
  }else{
    bW0 <- matrix(0, nrow=num.nonpen, ncol=dim(Y)[2])
    bw1 <- matrix(0, nrow=p-num.nonpen, ncol=dim(Y)[2])
    bx <- rbind(bw0, bW0, bw1)
  }
  if(intercept) bw0 <- matrix(0, ncol=dim(Y)[2])
  theta <- 1
  #ptm <- proc.time()
  for(iter in 1:option["maxiter"]){
    #compute grad(f(w_k))
    R <- shrink(C%*%t(bw1)/mu, g_idx)
    
    if(num.nonpen==0){
      if(intercept) grad_bw0 <- matrix(1, ncol=dim(X)[1]) %*% (matrix(1, nrow=dim(X)[1])%*%bw0 + X%*%bw1 - Y)
      grad_bw <- XX %*% bw1 - XY + t(R) %*% C
    }else{
      if(!is.null(bw0)) grad_bw0 <- matrix(1, ncol=dim(X)[1]) %*% (matrix(1, nrow=dim(X)[1])%*%bw0 + X%*%rbind(bW0,bw1) - Y)
      grad_bW0 <- crossprod(X[,1:num.nonpen], matrix(1, nrow=dim(X)[1])%*%bw0 + X%*%rbind(bW0,bw1) - Y)
      grad_bw <- crossprod(X[,-c(1:num.nonpen)], matrix(1, nrow=dim(X)[1])%*%bw0 + X%*%rbind(bW0,bw1) - Y) + t(R) %*% C
    }
    
    if(intercept) b0 <- bw0 - 1/L * grad_bw0
    if(num.nonpen>0) B0 <- bW0 - 1/L * grad_bW0
    bv <- bw1 - 1/L * grad_bw
    #browser()
    b <- abs(bv)-lambda/L
    b[b<0] <- 0
    b <- data.matrix(sign(bv) * b)
    bx_new <- data.matrix(rbind(b0, b))
    if(num.nonpen>0) bx_new <- data.matrix(rbind(b0, B0, b))
    
    if(intercept){
      obj_new <- sum((Y-cbind(rep(1,dim(X)[1]),X)%*%bx_new)^2)/2 + cal2norm(tcrossprod(C,b), g_idx)
    }else{
      obj_new <- sum((Y-X%*%bx_new)^2)/2 + cal2norm(tcrossprod(C,b), g_idx)
    }
    
    theta_new <- 2/(iter+2)
    
    bw <- bx_new + (1-theta)/theta * theta_new * (bx_new-bx)
    if(intercept) bw0 <- bw[1,]
    if(num.nonpen>0){
      bW0 <- bw[nrow(b0)+1:num.nonpen,]
      bw1 <- bw[-(1:(nrow(b0)+num.nonpen)),]
    }else{
      if(intercept){
        bw1 <- bw[-1,]
      }else{
        bw1 <- bw
      }
    }
    
    #time0 <- proc.time() - ptm
    #time[iter] <- time0$user
    
    if((iter>10) && (abs(obj_new-obj)/abs(obj)<option["tol"])) break
    theta <- theta_new
    bx <- bx_new
    obj <- obj_new
  }
  
  bx[abs(bx) < option["threshold"]] <- 0
  Beta <- bx
  #obj <- obj[1:iter]
  #time <- time[1:iter]
  #return(list(Beta=Beta, obj=obj, time=time, iter=iter))
  return(list(Beta=Beta))
}
# function accgrad2 is for updating coefficients of different data sources separately with different penalty factors
accgrad2 <- function(y, x, lambda, Tree, C, g_idx, TauNorm, p, mu, option, num.nonpen=0, intercept=TRUE){
  Tree <- data.matrix(Tree)
  Y <- data.matrix(y)
  X <- data.matrix(x)
  X <- cbind(rep(1, dim(X)[1]), X)
  X1 <- data.matrix(x[,num.nonpen+1:p[1]])
  X2 <- data.matrix(x[,-c(1:(num.nonpen+p[1]))])
  C <- data.matrix(C)
  
  # to be optimized for more data sources 
  XX1 <- crossprod(X1)
  XX2 <- crossprod(X2)
  X12 <- crossprod(X1, X2)
  X21 <- crossprod(X2, X1)
  XY1 <- crossprod(X1, Y)
  XY2 <- crossprod(X2, Y)
  
  L <- eigen(crossprod(X))$values[1] + TauNorm/mu
  
  J1 <- dim(X1)[2]
  J2 <- dim(X2)[2]
  K <- dim(Y)[2]
  
  C1 <- C * lambda[1]
  C2 <- C * lambda[2]
  if(intercept){
    bw0 <- matrix(0, ncol=dim(Y)[2])
  }else{
    b0 <- bw0 <- matrix(0, nrow=0, ncol=dim(Y)[2])
  }
  if(num.nonpen==0){
    bw1 <- bx_new1 <- matrix(0, nrow=J1, ncol=K)
    bw2 <- bx_new2 <- matrix(0, nrow=J2, ncol=K)
    bx <- rbind(bw0, bw1, bw2)
  }else{
    bW0 <- matrix(0, nrow=num.nonpen, ncol=K)
    bw1 <- bx_new1 <- matrix(0, nrow=J1, ncol=K)
    bw2 <- bx_new2 <- matrix(0, nrow=J2, ncol=K)
    bx <- rbind(bw0, bW0, bw1, bw2)
  }
  if(intercept) bw0 <- matrix(0, ncol=dim(Y)[2])
  theta <- 1
  obj <- 0
  #ptm <- proc.time()
  for(iter in 1:option["maxiter"]){
    #compute grad(f(w_k))
    R1 <- shrink(tcrossprod(C1,bw1)/mu, g_idx)
    R2 <- shrink(tcrossprod(C2,bw2)/mu, g_idx)
    
    if(intercept) grad_bw0 <- matrix(1, ncol=dim(X)[1]) %*% (matrix(1, nrow=dim(X)[1])%*%bw0 + X1%*%bw1 + X2%*%bw2 - Y)
    grad_bw1 <- XX1 %*% bw1 + X12 %*% bw2 - XY1 + crossprod(R1,C1)
    grad_bw2 <- XX2 %*% bw2 + X21 %*% bw1 - XY2 + crossprod(R2,C2)     
    
    if(intercept) b0 <- bw0 - 1/L * grad_bw0
    bv1 <- bw1 - 1/L * grad_bw1
    bv2 <- bw2 - 1/L * grad_bw2
    
    b1 <- abs(bv1)-lambda[1]/L
    b1[b1<0] <- 0
    bx_new1 <- data.matrix(sign(bv1) * b1)
    
    b2 <- abs(bv2)-lambda[2]/L
    b2[b2<0] <- 0
    bx_new2 <- data.matrix(sign(bv2) * b2)
    
    bx_new <- rbind(b0, bx_new1, bx_new2)
    obj_new <- sum((Y-X%*%bx_new)^2)/2 + cal2norm(tcrossprod(C1,bx_new1), g_idx) + cal2norm(tcrossprod(C2,bx_new2), g_idx)
    
    theta_new <- 2/(iter+2)
    
    bw <- bx_new + (1-theta)/theta * theta_new * (bx_new-bx)
    if(intercept) bw0 <- bw[1,]
    
    bw1 <- bw[nrow(b0)+1:J1,]
    bw2 <- bw[-c(1:(J1+nrow(b0))),]
    
    if((iter>10) && (abs(obj_new-obj)/abs(obj)<option["tol"])) break
    theta <- theta_new
    bx <- bx_new
    obj <- obj_new
  }
  
  bx[abs(bx) < option["threshold"]] <- 0
  Beta <- bx
  #obj <- obj[1:iter]
  #time <- time[1:iter]
  
  return(list(Beta=Beta))
}
cal2norm <- function(A, g_idx){
  s <- sum(apply(matrix(1:dim(g_idx)[1],ncol=1), 1, function(xx) sum(sqrt(sum(A[g_idx[xx,1]:g_idx[xx,2],]^2)))))
  # V <- dim(g_idx)[1]
  # s <- 0
  # for(v in 1:V){
  #   idx <- g_idx[v,1]:g_idx[v,2]
  #   gnorm <- sqrt(sum(A[idx,]^2))
  #   s <- s + sum(gnorm)
  # }
  return(s)
}
convH2T <- function(H, w_max){
  K <- dim(H)[1] + 1
  Nd <- cbind(rep((K+1):(2*K-1), each = 2), as.vector(t(H[,1:2])))
  W_norm <- H[,3]/max(H[,3])  
  conv0 <- convNd2T(Nd, W_norm, w_max)  
  
  return(conv0)  
}
convNd2T <- function(Nd, w, w_max){
  # Nd : node list
  # w : a vector of weights for internal nodes
  # Tree : VxK matrix
  #	V is the number of leaf nodes and internal nodes
  #	K is the number of tasks
  #	Element (v,k) is set to 1 if task k has a membership to
  #	the cluster represented by node v. Otherwise, it's 0.
  # Tw : V vector
  
  #===========================
  find_leaves <- function(Nd, ch, K, Jt, w, Tw){
    
    for(ii in 1:length(ch)){
      if(Nd[ch[ii], 2] > K){
        leaves0 <- find_leaves(Nd, which(Nd[,1] == Nd[ch[ii], 2]), K, Jt, w, Tw) 
        Jt <- leaves0$Jt
        Tw <- leaves0$Tw
      }else
        Jt <- c(Jt, Nd[ch[ii], 2])
    }
    
    Tw[Nd[ch, 2]] <- Tw[Nd[ch, 2]] * w
    
    return(list(Jt=Jt, Tw=Tw))
  }
  #===========================
  
  # of leaf nodes
  K <- Nd[1,1] - 1
  #V = Nd(size(Nd,1),1);
  #V = Nd(size(Nd,1),1)-1;		# without the root
  if(sum(w < w_max)<1){
    V <- 1 + K
  }else{
    ind0 <- which(w < w_max)    # only the internal nodes with w<w_max
    V <- ind0[length(ind0)] + K 
  }
  
  # for leaf nodes
  I <- 1:K
  J <- 1:K
  
  Tw <- rep(1, V)
  
  # for internal nodes
  for(i in (K+1):V){
    Jt <- NULL
    
    Tw[i] <- Tw[i] * (1 - w[i-K])
    leaves0 <- find_leaves(Nd, which(Nd[,1] == i), K, Jt, w[i-K], Tw)
    Jt <- leaves0$Jt
    Tw <- leaves0$Tw
    
    I <- c(I, rep(1,length(Jt)) * i)
    J <- c(J, Jt)
  }
  
  Tree <- sparseMatrix(i=I, j=J, x=rep(1, length(I)), dims=c(V, K))
  
  return(list(Tree=Tree, Tw=Tw))
  
}
fastCorr <- function(A){
  # n <- dim(A)[1]
  # B <- scale(A)
  C <- crossprod(scale(A))/(dim(A)[1]-1)
  return(C)
}
mheight <- function(y, foldid=NULL){
  if(is.null(foldid)){
    myDist0 <- 1 - abs(fastCorr(y))
    myDist <- myDist0[lower.tri(myDist0)]
    a0 <- dist(t(y))
    a0[1:length(a0)] <- myDist
    # hierarchical clustering for multivariate responses
    myCluster <- hclust(a0, method = "complete")
    min.height <- min(myCluster$height/max(myCluster$height))
  }else{
    min.height <- 0
    for(i in 1:max(foldid)){
      y0 <- y[!foldid==i,]
      myDist0 <- 1 - abs(fastCorr(y0))
      myDist <- myDist0[lower.tri(myDist0)]
      a0 <- dist(t(y0))
      a0[1:length(a0)] <- myDist
      # hierarchical clustering for multivariate responses
      myCluster <- hclust(a0, method = "complete")
      min.height <- max(min(myCluster$height/max(myCluster$height)), min.height)
    }
  }
  return(min.height)
}
pre.grad <- function(Tree, Tw){
  
  V <- dim(Tree)[1]
  K <- dim(Tree)[2]
  
  sum_col_T <- apply(Tree, 1, sum)
  SV <- sum(sum_col_T)
  csum <- cumsum(sum_col_T)
  #g_idx <- cbind(c(1,csum[1:(length(csum)-1)]+1), csum, sum_col_T)
  if(length(csum)!=1){
    g_idx <- cbind(c(1,csum[1:(length(csum)-1)]+1), csum, sum_col_T)
  }else{
    g_idx <- cbind(1, csum, sum_col_T)
  }
  
  J <- rep(0, SV)
  W <- rep(0, SV)
  for(v in 1:V){
    J[g_idx[v,1]:g_idx[v,2]] = which(Tree[v,] == 1)
    W[g_idx[v,1]:g_idx[v,2]] = Tw[v]
  }
  
  C <- sparseMatrix(i=1:SV, j=J, x=W, dims=c(SV, K))
  
  TauNorm0 <- matrix(Tw, ncol=1)
  for(r in 2:K) TauNorm0 <- cbind(TauNorm0, Tw)
  TauNorm <- TauNorm0 * Tree
  TauNorm <- max(apply(TauNorm^2, 2, sum))
  
  return(list(C=C, g_idx=g_idx, TauNorm=TauNorm))
}
# function accgrad2 calculates the norm term in Lipschitz constant of the irectly adapted tree-lasso method for two data sources
pre.grad2 <- function(Tree, Tw, lambda){
  
  if(is.null(dim(Tree))) Tree <- matrix(Tree, nrow=1)
  
  V <- dim(Tree)[1]
  K <- dim(Tree)[2]
  
  sum_col_T <- apply(Tree, 1, sum)
  SV <- sum(sum_col_T)
  csum <- cumsum(sum_col_T)
  if(length(csum)!=1){
    g_idx <- cbind(c(1,csum[1:(length(csum)-1)]+1), csum, sum_col_T)
  }else{
    g_idx <- cbind(1, csum, sum_col_T)
  }
  
  J <- rep(0, SV)
  W <- rep(0, SV)
  for(v in 1:V){
    J[g_idx[v,1]:g_idx[v,2]] = which(Tree[v,] == 1)
    W[g_idx[v,1]:g_idx[v,2]] = Tw[v]
  }
  
  C <- sparseMatrix(i=1:SV, j=J, x=W, dims=c(SV, K))
  
  TauNorm0 <- matrix(Tw, ncol=1)
  for(r in 2:K) TauNorm0 <- cbind(TauNorm0, Tw)
  TauNorm <- TauNorm0 * Tree
  TauNorm <- rbind(lambda[1]*TauNorm, lambda[2]*TauNorm)
  TauNorm <- max(apply(TauNorm^2, 2, sum))
  
  return(list(C=C, g_idx=g_idx, TauNorm=TauNorm))
}
shrink <- function(A, g_idx){
  V <- dim(g_idx)[1]
  R <- matrix(0, nrow=dim(A)[1], ncol=dim(A)[2])
  for(v in 1:V){
    idx <- g_idx[v,1]:g_idx[v,2]
    gnorm <- sqrt(apply(A[idx,]^2, 2, sum))
    gnorm[gnorm<1] <- 1
    gnorm0 <- gnorm
    for(r in 2:g_idx[v,3]) gnorm0 <- rbind(gnorm0, gnorm)
    R[idx,] <- data.matrix(A[idx,]/gnorm0)
  }
  return(R)
}
# function ipf.tree.lasso is the directly adapted tree-lasso method for two data sources
ipf.tree.lasso <- function(x, y, p=NULL, h=0.7, lambda=rep(10,2), tree.parm, num.nonpen=0, intercept=TRUE, mu=0.01, threshold=0){
  X <- data.matrix(x)
  Y <- data.matrix(y)
  m <- dim(Y)[2]
  option <- c(10000, threshold, 1e-6)
  names(option) <- c("maxiter", "threshold", "tol")
  
  pre0 <- pre.grad2(tree.parm$Tree, tree.parm$Tw, lambda)
  
  acc <- accgrad2(Y, X, lambda, tree.parm$Tree, pre0$C, pre0$g_idx, pre0$TauNorm, p, mu, option)
  Beta <- acc$Beta
  obj <- acc$obj
  time <- acc$time
  
  return(Beta)
}
tree.parms <- function(y=y, h=0.7){
  m <- dim(y)[2]
  myDist0 <- 1 - abs(fastCorr(y))
  myDist <- myDist0[lower.tri(myDist0)]
  a0 <- dist(t(y))
  a0[1:length(a0)] <- myDist
  # hierarchical clustering for multivariate responses
  myCluster <- hclust(a0, method = "complete")
  myCluster <- cbind(ifelse(myCluster$merge < 0, - myCluster$merge, myCluster$merge + m), myCluster$height)
  
  conv0 <- convH2T(myCluster, h)
  Tree <- conv0$Tree
  if(is.null(dim(Tree)))
    Tree <- matrix(Tree, nrow=1)
  Tw <- conv0$Tw
  idx <- c(apply(Tree,1,sum) == 1)
  Tree <- Tree[!idx,]
  if(is.null(dim(Tree)))
    Tree <- matrix(Tree, nrow=1)
  Tw <- Tw[!idx]
  
  return(list(Tree=Tree, Tw=Tw))
}
tree.lasso <- function(x=x, y=y, lambda=10, tree.parm, num.nonpen=0, intercept=TRUE, mu=0.01, threshold=0){
    x <- data.matrix(x)
    y <- data.matrix(y)
    option <- c(10000, threshold, 1e-6)
    names(option) <- c("maxiter", "threshold", "tol")
    pre0 <- pre.grad(tree.parm$Tree, tree.parm$Tw)
    Beta <- accgrad(y, x, lambda, tree.parm$Tree, pre0$C, pre0$g_idx, pre0$TauNorm, mu, option, num.nonpen=num.nonpen, intercept=intercept)$Beta
    return(Beta)
}





