#' IPFStructPenalty
#' @title Efficient Parameter Selection via Global Optimization
#' @description
#' Finds an optimal solution for the \code{Q.func} function.
#' @param Q.func name of the function to be minimized.
#' @param bounds bounds for the interval-searching parameters
#' @param x,y input matrix.
#' @param family response type.
#' @param lambda optional user-supplied \code{lambda} sequence; default is NULL, and \code{espsgo} chooses its own sequence.
#' @param intercept should  intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE}).
#' @param strata.surv stratification variable for the Cox survival model.
#' @param threshold threshold for estimated coefficients of the tree-lasso models.
#' @param foldid an vector of values for the cross-validation.
#' @param standardize.response standardization for the response variables. Default: \code{TRUE}.
#' @param p the numbers of predictors from different data sources.
#' @param num.nonpen number of predictors forced to be estimated (i.e., nonpenalization).
#' @param round.n number of digits after comma, default is \code{5}.
#' @param parms.coding parmeters coding: none or log2, default: \code{none}.
#' @param fminlower minimal value for the function Q.func, default is 0.
#' @param flag.find.one.min do you want to find one min value and stop? Default: \code{FALSE}.
#' @param show show plots of DIRECT algorithm: none, final iteration, all iterations. Default: \code{none}.
#' @param N define the number of start points depending on the dimensionality of the parameter space.
#' @param maxevals the maximum number of DIRECT function evaluations, default: 500.
#' @param espilon the convergence shreshold for the function \code{Q.func}, default is 0.01.
#' @param EI.eps the convergence threshold for the expected improvement between fmin and the updated point 
#' @param min.iter the minimus iterations after the initial \code{N} iterations.
#' @param verbose print the middle search information, default is \code{TRUE}.
#' @param seed random seed
#' @param parallel If \code{TRUE}, use parallel foreach to fit each fold except parallelizing each lambda for the tree-lasso methods. If \code{c(TRUE,TRUE)}, use parallel foreach to fit each fold and each lambda.
#' @param search.path save the visited points, default is \code{FALSE}.
#' @param modelList detailed information of the search process
#' @references Frohlich, H. & Zell, A. (2005). \emph{Efficient Parameter Selection for Support Vector Machines in Clas- sification and Regression via Model-Based Global Optimization.} Proceedings of the International Joint Conference of Neural Networks, pp 1431-1438.
#' @references Sill, M., Hielscher, T., Becker, N. & Zucknick, M. (2014).\emph{c060: Extended Inference with Lasso and elastic net Regularized Cox and Generalized Linear methods.} Journal of Statistical Software, 62(5):1-22.
#' @export
epsgo<- function(
								# function to be minimized
								Q.func, 
								# bounds for parameters
								bounds, 
								x, y,
								family="gaussian",
								lambda=NULL,
								p=NULL,
								intercept=TRUE,
								foldid=NULL,
								# number of first nonpenalised features
								num.nonpen = 0,
								# vector with stratum membership of each observation for conditional logistic lasso
								strata.surv=NULL,
								# threshold for estimated coefficients of the IPF-tree-lass
								threshold=0,
								standardize.response=FALSE,
                # round.n -number of digits after comma
                round.n=5,
								parms.coding="none", # or log2 
								# min value for the function Q.func
								fminlower=0, 
								# do you want to find one min value and stop?
								flag.find.one.min =FALSE,
								# show plots ?   none, final iteration, all iterations 
								#show=c("none", "final", "all"),
								show="none",
								# define the number of start points
								N= NULL, 
								#% maximum # of function evaluations 
								maxevals =   500,
                constantMean=0,
                #% convergent condition
                epsilon = 1e-4,
								Dir.ep = 1e-4,
								Dir.tol = 0.01,
								# threshold different of improvement between fmin and the updated point 
								EI.eps = 0.01,
                #% No changes in the last min.iter iterations, break iterations
                min.iter = 10,
								# plot parameter
                pdf.name=NULL, 
                pdf.width=12, pdf.height=12,
                my.mfrow=c(1,1), 
								parallel=FALSE,
								modelList=NULL,
                # verbose?
                verbose=TRUE,
                seed=123, 
								search.path=FALSE,
                ... ){
	# The EPSGO algorithm (theory from Frohlich and Zell (2005) and original code from Sill, Hielscher, Becker and Zucknick (2014))
	#
	#Input: function Q.func to measure gen. error
	#			  bounds: data frame with parameter bounds
	#								rownames: parameter names
	#								colnames: "lower" and "upper"
	#								Example: bounds=t(data.frame(lambda1=c(0, 10), lambda2=c(0,8)));colnames(bounds)<-c("lower", "upper")
	#								fminlower:  min Q value				
	#Output: Qmin, p*, number of visited points neval 
	
	## Scheme:
	#	1. number of tuning patameters
	#		 D = dim(P)
	#	2. create N = 10D sample points p1, ..., pN
	#		in [bounds] using Latin hypercube sampl.
	#	3. compute X= Q(p_i), i = 1, ...,N
	# 
	# 4. train Online GP
	#	5. number of visited points p_obs 
	#	 	 neval = N
	#	6.	REPEAT
	#					6.1  Find a new p, with max E[I(p)] ( the same as min -E[I(p)] )  
	#					?p = argmaxp E[I(p)] (computed by DIRECT)
	#					# Important! Direct.R calculate global min! --> change Problem function
	#					6.2 compute std. dev. and mean of E[I(p)]
	#					6.3 new p, new Q(p) 
	#						  Qnew = Q(?p)
	#					6.4 Add the new p ?
	#							if Qnew < Qmin
	#								Qmin = Qnew
	#								?p = ?p
	#							end
	#					6.5 update Online GP
	#							neval = neval + 1
	#			UNTIL convergence
	
	if (verbose){
	  print("parms.coding")
	  print(parms.coding)
	}
	 
	###################################################################################################################
		## 1. ## number of tuning patameters D = dim(P)
	###################################################################################################################
	
	D<- nrow(bounds)
	
	###################################################################################################################
		## 2. ##  create N = 10D sample points p1, ..., pN
		#in [l, u] using Latin hypercube sampl.
	###################################################################################################################
	
	set.seed(seed)
	# in ego.m  X = start.points
	# N - number of start points in parameter space
	# wie in ego.m
	
	# define the number of start points 
	if (is.null(N)){
		ns<-c(21, 21, 33, 41, 51, 65, 65)  ##
		# N = 10D, D= number of parameters , but for high dim data with more than 6 dim, restrict the initial set of p to 65
		N<- ifelse ( D <= length(ns),ns[D], 65 )
	}
	
	# start points X (= p_1,..., p_N)
	
	X<- lhs(N, bounds)
    colnames(X)<-rownames(bounds)
	# round.n -number of digits after comma
	X<-round(X,round.n)
	
	if (verbose) print(X)
	
	if(show !="none" & !is.null(pdf.name)) { 
		pdf(pdf.name, pdf.width, pdf.height)
		par(mfrow=my.mfrow)
	}

	if ( (show !="none") & (D<=2 )){
		# 1D plot
		if (D==1) {
			plot(X, xlab="Index", ylab=rownames(bounds)[1], col="orange", pch=20,
						 main=paste( "Latin hypercube sampling, n=",nrow(X)*D) )
			abline(h=seq(bounds[1,1],bounds[1,2],length=(nrow(X)+1)), lty=2, col=3)
		} else {
			# 2D plot
			plot(X, xlab=rownames(bounds)[1], ylab=rownames(bounds)[2], col="orange", pch=20, 
						 main=paste( "Latin hypercube sampling, n=",nrow(X)))
			abline(v=seq(bounds[1,1],bounds[1,2],length=(nrow(X)+1)), lty=2, col=3)
			abline(h=seq(bounds[2,1],bounds[2,2],length=(nrow(X)+1)), lty=2, col=3)
		}
	}


	###################################################################################################################
	## 3. ## compute Q(p_i), i = 1, ...,N
	###################################################################################################################
	
	model.list<-apply(X, 1, eval(Q.func), parms.coding=parms.coding, x=x, y=y, strata.surv=strata.surv, lambda=lambda, maxevals=maxevals, seed=seed, family=family, 
	                    num.nonpen=num.nonpen,foldid=foldid, intercept=intercept, standardize.response=standardize.response, p=p, parallel=parallel,verbose=verbose,search.path=search.path,threshold=threshold, ...)
	#browser()
	# take the Q.values
	Q<- as.numeric(unlist( sapply(model.list, "[", "q.val")))
	Q.min <- min (Q, na.rm=T)
	min.p<- X[which.min(Q), , drop=FALSE]
	if (verbose) print(data.frame(X,Q))
	
	# delete the point(s) with no Q value(s)
	X <- X[!is.na(Q), ,drop=FALSE]
	Q <-Q[!is.na(Q)]
	
	# add start.q.values to the plot
	if (show !="none" & D ==1)  text( c(1: nrow(X)),X, labels=round(Q,6), pos=1, cex=0.5 )
	if (show !="none" & D ==2)  text( X, labels=signif(Q,6), pos=1, cex=0.5 )
	
	#		# 4. train Online GP
	
	# train Gaussian Process 
	# Input:  collection of random variables G(x)(here: x= points in tuning parameter space = start.points X ) 
	#
	###################################################################################################################
	## 4. ## train Gaussian Process
	###################################################################################################################
	
		gp.seed.new<- seed 
			
		# if we have tried 5 times and are still not able to fit --> break
			# the reason is in having a new point in Ytrain very close to one of the old ones.
			# --> matrix is singular! 
			
			if (exists("fit.gp")) rm(fit.gp)
			flag.fit.gp<- FALSE
      tmp.i<-1
      
     	while( ! flag.fit.gp   ) {
				
				if (tmp.i >5) {
					print(print( "At least one of the intial points is very close to the other points. It is not possible to fit the model via Gaussian Process. "))
					break  
				}
				try(fit.gp<-mlegp(X, Q,constantMean=constantMean, seed=gp.seed.new, verbose=0))
				
				flag.fit.gp<- FALSE
      	# if fit.gp exists AND is not null
      	if (exists("fit.gp"))
      		if (!is.null(fit.gp))
      				flag.fit.gp<- TRUE
      
				# if fails to fit change seed
				if(!flag.fit.gp  ) {
					 print("fails to fit gp (fit.gp), change seed !")
					 gp.seed.new<- round(runif(1,min=1, max=10^3))
					 tmp.i<-tmp.i + 1
				}
			} # end of while 
			
			gp.seed<- gp.seed.new
	
	###################################################################################################################
	#	5. number of visited points p_obs 
	###################################################################################################################
	neval <- length(Q)
	
	###################################################################################################################
	##.6.	REPEAT - main block
	###################################################################################################################
	
	# Initialization for observed points in parameter space Xtrain in R^D and 
	# 							 for quality value Ytrain = Q		  
	Xtrain<- X
	Ytrain <- Q
	finished <- FALSE
	EImax <- Inf
	
	loop<-1
	fcalls<-length(Q)
	# Point with max E[I(p)]
	xmax=X[1, ,drop=FALSE]
	#fmin <- Inf
	# change to already calculated!
		
	# info for point with min Q.func
	fminold<-  Inf 
	fmin <- Q.min
	xmin<- min.p
	not_changed = 0
	
	#set.seed(seed)
	
	while (!finished){ 
		if (verbose ) print(paste("loop", loop))
		# calculate Q for the new points 		
		if (loop >1) {
      EIold <- EImax	
		  fminold <- fmin
		
			model.list.new<-apply(X, 1, eval(Q.func), parms.coding=parms.coding, x=x, y=y, strata.surv=strata.surv, lambda=lambda, maxevals=maxevals, seed=seed, family=family, 
		                          num.nonpen=num.nonpen,foldid=foldid, intercept=intercept, standardize.response=standardize.response, spatial=spatial, p=p, parallel=parallel,verbose=verbose,search.path=search.path,threshold=threshold, ...)
		  
      model.list<-c(model.list, model.list.new )
      Q<- as.numeric(unlist( sapply(model.list.new, "[", "q.val")))
			
			fcalls = fcalls + length(Q);
			Xtrain = rbind(Xtrain, X)
			Ytrain = c(Ytrain, Q )	
			
      if (verbose) print(data.frame(Xtrain,Ytrain))
			
			# fmin =  current min of Q.func at point xmin. 
			fmin<- min(Ytrain)
			xmin = Xtrain[which.min(Ytrain),]	
		}
		
		if (verbose) print(paste("fmin=",fmin))
		if (fmin < fminlower){ 
				break
			} 
		
		if (flag.find.one.min){
			### skip it (in case of having more than 1 points with global min )
			# break if reach the min value
			if (fmin <= fminlower){ 
				break
			} 
		}
		
		
		# break if no changes in the last 10 iterations
		if (fmin == fminold){
			not_changed = not_changed + 1
			if (not_changed >= min.iter){
				print(paste("No changes in the last", min.iter, "iterations, break iterations"))
				break
			}
		}else {
			not_changed = 0
		} 
	  
		# train GP
		if (loop >1) {
			# sometimes error:  Error in solve.default(gp$invVarMatrix) :
	  	#					system is computationally singular: reciprocal condition number = 7.502e-17
			# Solution: change seed
						
			gp.seed.new<- c(seed + loop-1) 
			
			
			if (exists("fit.gp")) rm(fit.gp)
			
			# if we have tried 5 times and are still not able to fit --> break
			# the reason is in having a new point in Ytrain very close to one of the old ones.
			# --> matrix is singular! 
			
		  tmp.i<-1
		  flag.fit.gp<- FALSE
      
     	while(!flag.fit.gp   ) {
				
				if (tmp.i >5) {
					print(print( "The new point X is very close to the one of visited points."))
					finished <- TRUE
					break  
				}
				try(fit.gp<-mlegp(Xtrain, Ytrain,constantMean=constantMean, seed=gp.seed.new, verbose=0))
				
				flag.fit.gp<- FALSE
      	# if fit.gp exists AND is not null
      	if (exists("fit.gp"))
      		if (!is.null(fit.gp))
      				flag.fit.gp<- TRUE
      
				# if fails to fit change seed
				if(!flag.fit.gp  ) {
					 print("fails to fit gp (fit.gp), change seed !")
					 gp.seed.new<- round(runif(1,min=1, max=10^3))
					 tmp.i<-tmp.i + 1
				}
			} # end of while  (!flag.fit.gp   ) 
			
			if (finished) break
			
			gp.seed<-c(gp.seed, gp.seed.new)
			#str(fit.gp)
		}
		
		if (show !="none" ) plot(fit.gp, main=paste("Gaussian Process", "iter ", loop))
			
		Problem<- list(f = "ExpImprovement" )
					
		#6.1  Find a new p, with max E[I(p)] ( the same as min -E[I(p)] ) 
		#[EImax xmax history] = Direct(Problem,bounds,options) 
		Dir.list<- Direct(Problem=Problem, 
								bounds=bounds, 
								# options 
								#% maximum of iterations 
								maxits =     50,   
								#% maximum # of function evaluations 
								maxevals = maxevals,       
								# the optimal value is unknown
								testflag= 0,
								# minimum value of function , min (ExpImprovement) = 0
								globalmin =  0, 
								#% print and plot iteration stats
								showits= show, 
								# verbose? 
								verbose=FALSE,
								ep= Dir.ep,
								tol= Dir.tol, 
								pdf.name=NULL, #"test.pdf" # 
								pdf.width=12, pdf.height=12,
								my.mfrow=c(1,1), #c(3,3),
								#
								# additional args for the problem function 
								fmin=fmin,
								fit.gp=fit.gp
								 )
		
		# E[I(p)], xmax - p with max EI
		EImax <- (-1 ) * Dir.list$minval
		xmax<-  Dir.list$final_point.xatmin
		History <-  Dir.list $ History
		
		if (verbose ){
			print(paste("EImax=", EImax ) )
			print( "xmax:")
			print( xmax)
			print( "history:")
			print( History)
		}
		
		# 6.2 compute std. dev. and mean of E[I(p)]
		### other random sample in the same size; Latin hypercube sampling over the 
		## whole parameter space. The idea is to make sure, that the expected improvement over the whole space is almost equally small. 
		set.seed(seed + sum(gp.seed))
		Xsam <- lhs(N, bounds)
		EIall=rep(0,N)	
		for (i in 1:N ){		
			EIall[i] = - ExpImprovement(Xsam[i, , drop=FALSE], fmin, fit.gp, muX=NULL, muY=NULL, EI.eps=EI.eps)
		}
		
		EIstd =  sd( c(EIall, EImax) )
		EImean = mean(c(EIall, EImax))
		
		# if E[I(p)] is max by a random sample --> take it!
		if ( max(EIall) > EImax ){
			xmax = t(Xsam[which.max(EIall),, drop=FALSE])
			rownames(xmax)<- rownames(bounds)
			EImax <- max(EIall)
		}
		
		# 6.4 Add the new p ?
		# if the new p is one from the observed ones --> break
		# OR if the differences in functions is small,
		#  stop iterations
		 
		visited <-  Xtrain - matrix(xmax, nrow=nrow(Xtrain), ncol= ncol(Xtrain), byrow=TRUE) 
		visited<- any (	apply (visited == 0 , 1, all )  )
		
		if ( (visited )   | ((abs(fmin - fminold) < epsilon) & ((EImax - EImean)^2 <= 0.1*EIstd) & (not_changed >= min.iter)) ){
			if (visited) print( "The new point with min E[I(p)] is already in the set of visited points.")
			if ((abs(fmin - fminold) < epsilon) & ((EImax - EImean)^2 <= 0.1*EIstd) ) {
				print( "the differences in functions between 2 last iterations is small, stop iterations") 
			}
			finished = TRUE
		}else{
			X = t(xmax)
      X<-round(X,round.n)
			rownames(X)<- NULL
			neval<- neval + 1
		}
        if (verbose ){
            print(paste("iteration :", loop, "  fmin = ", fmin ))
            print(paste("finished?", finished))
            print("X")
            print(X)
            loop = loop + 1
		
            print(data.frame(Xtrain,Ytrain))
        }
	
	} # end of while (!finished) 
		
	if(show !="none" & !is.null(pdf.name)) dev.off()	
		
	# define the set of points with the same fmin
		tmp.set<-data.frame(Xtrain, f=Ytrain)
		points.fmin<- tmp.set[tmp.set$f == fmin, ,drop=FALSE  ]
	
  out <- list(fmin =fmin,
	     xmin = xmin, 
	     iter = loop,
	     neval =neval,
	     maxevals= maxevals,
	     seed =seed,
	     bounds = bounds,
	     Q.func =Q.func,
	     points.fmin =points.fmin,
	     Xtrain = Xtrain, 
	     Ytrain= Ytrain, 
	     gp.seed=gp.seed,
	     model.list = model.list 
	) 
  
  class(out) <- "intsearch"
	return(out)
	
}
