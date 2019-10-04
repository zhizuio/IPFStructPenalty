#================================================================================================================
# This script includes three functions to simulate two heterogeneous data sources (covariates) with multivariate response variables.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 08-Feb-2019
#================================================================================================================

sim2 <- function(p=c(500,150),n=100,m=24,rho=.4,B.elem=c(0.2,0.6)){
  b<-10
  if(!is.na(p[2])){
    # generate covariance matrix
    Ap1<-matrix(rep(rho,(p[1]/b)^2),nrow=p[1]/b)
    diag(Ap1)<-rep(1,p[1]/b)
    Ap2<-matrix(rep(rho,(p[2]/b)^2),nrow=p[2]/b)
    diag(Ap2)<-rep(1,p[2]/b)
    Bp12<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[1]/b)
    Bp21<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[2]/b)
    Xsigma1<-Ap1
    Xsigma2<-Ap2
    Xsigma12<-Bp12
    Xsigma21<-Bp21
    for(i in 2:b){
      Xsigma1<-bdiag(Xsigma1,Ap1)
      Xsigma12<-bdiag(Xsigma12,Bp12)
      Xsigma2<-bdiag(Xsigma2,Ap2)
      Xsigma21<-bdiag(Xsigma21,Bp21)
    }
    
    Xsigma<-rbind(cbind(Xsigma1,Xsigma12),cbind(Xsigma21,Xsigma2))
    X<-	mvrnorm(n,mu=rep(0,p[1]+p[2]),Sigma=Xsigma)
    #X1<-X[,1:p[1]]
    #X2<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    X[,(p[1]+1):(p[1]+p[2])]<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    # generate uncorrelated error term
    esd<-diag(m)
    e<-mvrnorm(n,mu=rep(0,m),Sigma=esd)
    
    ## generate beta1 matrix
    Beta1<-matrix(0,nrow=m,ncol=p[1])
    Beta1[,1]<-B.elem[1]
    for(i in 1:2){
      Beta1[((i-1)*m/2+1):(i*m/2),(1+(i-1)*2+1):(1+i*2)]<-B.elem[1]
    }
    for(i in 1:4){
      Beta1[((i-1)*m/4+1):(i*m/4),(1+4*2+(i-1)*4+1):(1+4*2+i*4)]<-B.elem[1]
    }
    for(i in 1:8){
      Beta1[((i-1)*m/8+1):(i*m/8),(1+2*2+4*4+(i-1)*8+1):(1+2*2+4*4+i*8)]<-B.elem[1]
    }
    ## generate beta2 matrix
    Beta2<-matrix(0,nrow=m,ncol=p[2])
    Beta2[,1]<-B.elem[2]
    for(i in 1:3){
      Beta2[((i-1)*m/3+1):(i*m/3),(1+(i-1)*2+1):(1+i*2)]<-B.elem[2]
    }
    for(i in 1:6){
      Beta2[((i-1)*m/6+1):(i*m/6),(1+2*3+(i-1)*4+1):(1+2*3+i*4)]<-B.elem[2]
    }
    for(i in 1:12){
      Beta2[((i-1)*m/12+1):(i*m/12),(1+2*3+4*6+(i-1)*8+1):(1+2*3+4*6+i*8)]<-B.elem[2]
    }
    Beta<-t(cbind(Beta1, Beta2))
  }else{
    cat("Ooops!!! Please specify 2-dim p vector, for example p=c(500,150)\n")
  }
  
  Y<-X%*%Beta+e
  return(list(Y=Y, X=X, Beta=Beta, e=e, p=p))
}

sim1 <- function(p=c(500,150),n=100,m=24,rho=.4,B.elem=c(0.2,0.6)){
  b<-10
  if(!is.na(p[2])){
    # generate covariance matrix
    Ap1<-matrix(rep(rho,(p[1]/b)^2),nrow=p[1]/b)
    diag(Ap1)<-rep(1,p[1]/b)
    Ap2<-matrix(rep(rho,(p[2]/b)^2),nrow=p[2]/b)
    diag(Ap2)<-rep(1,p[2]/b)
    Bp12<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[1]/b)
    Bp21<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[2]/b)
    Xsigma1<-Ap1
    Xsigma2<-Ap2
    Xsigma12<-Bp12
    Xsigma21<-Bp21
    for(i in 2:b){
      Xsigma1<-bdiag(Xsigma1,Ap1)
      Xsigma12<-bdiag(Xsigma12,Bp12)
      Xsigma2<-bdiag(Xsigma2,Ap2)
      Xsigma21<-bdiag(Xsigma21,Bp21)
    }
    
    Xsigma<-rbind(cbind(Xsigma1,Xsigma12),cbind(Xsigma21,Xsigma2))
    X<-	mvrnorm(n,mu=rep(0,p[1]+p[2]),Sigma=Xsigma)
    #X1<-X[,1:p[1]]
    #X2<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    X[,(p[1]+1):(p[1]+p[2])]<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    # generate uncorrelated error term
    esd<-diag(m)
    e<-mvrnorm(n,mu=rep(0,m),Sigma=esd)
    
    ## generate beta1 matrix
    Beta1<-matrix(0,nrow=m,ncol=p[1])
    Beta1[,1]<-B.elem[1]
    for(i in 1:2){
      Beta1[((i-1)*m/2+1):(i*m/2),(1+(i-1)*2+1):(1+i*2)]<-B.elem[1]
    }
    for(i in 1:2){
      Beta1[((i-1)*m/4+1):(i*m/4),(1+4*2+(i-1)*4+1):(1+4*2+i*4)]<-B.elem[1]
    }
    for(i in 1:4){
      Beta1[((i-1)*m/8+1):(i*m/8),(1+2*2+4*4+(i-1)*8+1):(1+2*2+4*4+i*8)]<-B.elem[1]
    }
    ## generate beta2 matrix
    Beta2<-matrix(0,nrow=m,ncol=p[2])
    Beta2[,1]<-B.elem[2]
    for(i in 1:2){
      Beta2[((i-1)*m/2+1):(i*m/2),(1+(i-1)*2+1):(1+i*2)]<-B.elem[2]
    }
    for(i in 3:4){
      Beta2[((i-1)*m/4+1):(i*m/4),(1+2*2+(i-1)*4+1):(1+2*2+i*4)]<-B.elem[2]
    }
    for(i in 5:8){
      Beta2[((i-1)*m/8+1):(i*m/8),(1+2*2+4*4+(i-1)*8+1):(1+2*2+4*4+i*8)]<-B.elem[2]
    }
    Beta<-t(cbind(Beta1, Beta2))
  }else{
    cat("Ooops!!! Please specify 2-dim p vector, for example p=c(500,150)\n")
  }
  
  Y<-X%*%Beta+e
  return(list(Y=Y, X=X, Beta=Beta, e=e, p=p))
}

sim3 <- function(p=c(500,150),n=100,m=24,rho=.4,B.elem=c(0.2,0.6)){
  b<-10
  B.elem[1]<-0.2
  B.elem[2]<-0.6
  if(!is.na(p[2])){
    # generate covariance matrix
    Ap1<-matrix(rep(rho,(p[1]/b)^2),nrow=p[1]/b)
    diag(Ap1)<-rep(1,p[1]/b)
    Ap2<-matrix(rep(rho,(p[2]/b)^2),nrow=p[2]/b)
    diag(Ap2)<-rep(1,p[2]/b)
    Bp12<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[1]/b)
    Bp21<-matrix(rep(rho,p[1]/b*p[2]/b),nrow=p[2]/b)
    Xsigma1<-Ap1
    Xsigma2<-Ap2
    Xsigma12<-Bp12
    Xsigma21<-Bp21
    for(i in 2:b){
      Xsigma1<-bdiag(Xsigma1,Ap1)
      Xsigma12<-bdiag(Xsigma12,Bp12)
      Xsigma2<-bdiag(Xsigma2,Ap2)
      Xsigma21<-bdiag(Xsigma21,Bp21)
    }
    
    Xsigma<-rbind(cbind(Xsigma1,Xsigma12),cbind(Xsigma21,Xsigma2))
    X<-	mvrnorm(n,mu=rep(0,p[1]+p[2]),Sigma=Xsigma)
    #X1<-X[,1:p[1]]
    #X2<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    X[,(p[1]+1):(p[1]+p[2])]<-data.matrix(X[,(p[1]+1):(p[1]+p[2])] > 0) + 0
    # generate uncorrelated error term
    esd<-diag(m)
    e<-mvrnorm(n,mu=rep(0,m),Sigma=esd)
    
    ## generate two coefficient matrices
    Beta1<-matrix(0,nrow=m,ncol=p[1])
    Beta2<-matrix(0,nrow=m,ncol=p[2])
    Beta1[1:5,5:6] <- B.elem[1]
    Beta1[15:20,5:6] <- B.elem[1]
    Beta1[3:6,10:11] <- B.elem[1]
    Beta1[2:8,30:31] <- B.elem[1]
    Beta1[16:21,30:31] <- B.elem[1]
    Beta1[20:m,26:27] <- B.elem[1]
    Beta1[7:13,90:91] <- B.elem[1]
    Beta1[9:17,100:101] <- B.elem[1]
    
    Beta2[1:7,5:6] <- B.elem[2]
    Beta2[5:6,7:8] <- B.elem[2]
    Beta2[6:11,12:13] <- B.elem[2]
    Beta2[8:9,14:15] <- B.elem[2]
    Beta2[10:13,20:21] <- B.elem[2]
    Beta2[14:15,26:27] <- B.elem[2]
    Beta2[16:17,36:37] <- B.elem[2]
    Beta2[18:19,46:47] <- B.elem[2]
    Beta2[20:21,56:57] <- B.elem[2]
    Beta2[22:23,66:67] <- B.elem[2]
    Beta2[11:12,90:91] <- B.elem[2]
    Beta2[3:14,100:101] <- B.elem[2]
    Beta2[3:14,100:101] <- B.elem[2]
    Beta2[13:18,110:111] <- B.elem[2]
    Beta2[17:24,120:121] <- B.elem[2]
    Beta <- t(cbind(Beta1, Beta2))
  }else{
    cat("Ooops!!! Please specify 2-dim p vector, for example p=c(500,150)\n")
  }
  
  Y<-X%*%Beta+e
  return(list(Y=Y, X=X, Beta=Beta, e=e, p=p))
}

