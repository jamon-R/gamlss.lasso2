# fit glmnet terms using the glmnet function of glmnet package
# which is used in the backfitting 
# TO DO
# (i)  is scalling the X always needed it?
# (ii)  at the moment X is a matrix should we allow formula?
# (iii)
#--------------------------------------------------------------------------------------
require(glmnet)
gnet <- function( 
                 X = NULL, 
            x.vars = NULL,
						lambda = NULL,	
						method = c("IC", "CV"),
						type = c("agg", "sel"),
						ICpen = c("BIC", "HQC", "AIC"),
						CVp = 2,
						k.se = 0,
						adaptive = 1,
						epsilon = 1/sqrt(dim(X)[1]),
            subsets = NULL,
						sparse = FALSE,
          	control = gnet.control(...), ...) {
  #------------------------------------------
  # function starts here
  #------------------------------------------
     scall <- deparse(sys.call(), width.cutoff = 500L)
    method <- match.arg(method)
  # check for standarized matric
     rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1) { 
      position <- i # get the position
      if (rexpr[i]==TRUE) break
    }
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
## get the data
if (sys.call(position)[1]=="predict.gamlss()") { # if predict is used 
      Data <- get("data", envir=gamlss.env)
    } else if (sys.call(position)[1]=="gamlss()") { # if gamlss() is used
      if (is.null(get("gamlsscall", envir=gamlss.env)$data)) { # if no data argument but the formula can be interpreted
        Data <- data.frame(X)	
    } else {# data argument in gamlss 
        Data <- get("gamlsscall", envir=gamlss.env)$data
    }
    } else {
      Data <- get("data", envir=gamlss.env)
    }
    Data <- data.frame(eval(substitute(Data)))

## GET THE DATA NNET
#if (sys.call(position)[1]=="predict.gamlss()")
#{ # if predict is used 
#  Data <- get("data", envir=gamlss.env)
#}
#else if (sys.call(position)[1]=="gamlss()") 
#{ # if gamlss() is used
#  if (is.null(get("gamlsscall", envir=gamlss.env)$data)) 
#  { # if no data argument but the formula can be interpreted
#    Data <- model.frame(formula)  
#  }
#  else
#  {# data argument in gamlss 
#    Data <- get("gamlsscall", envir=gamlss.env)$data
#  }
#}
#else  {Data <- get("data", envir=gamlss.env)}
#Data <- data.frame(eval(substitute(Data)))

  # prepare glmnet
#	if (sys.call(position)[1]!="predict.gamlss()"){ # not predict.gamlss
	if (is.null(X)&&is.null(x.vars)) stop("X or x.vars has to be set in gnet")
	if (is.null(X)&&!is.null(x.vars)) X <- as.matrix(Data[, x.vars])
  if (!is.null(X)&&is.null(x.vars)) warning("For prediction use the x.vars argument")

#	}
	if( sparse ) {
		X <- Matrix::Matrix(X, sparse = TRUE)
	} 
#else {  
#		X <- scale(X, center=) #FZ: is done in glmnet, could be done here to save time but does not work for sparse matrices, as it would destroy the sparsity structure...
#	}
  nobs <- dim(X)[1]
  nvars <- dim(X)[2]
	
#  formula = "glmnet(x = X[subsets[[i.nets]],], y = y[subsets[[i.nets]]], weights = w[subsets[[i.nets]]], lambda=lambda,"

#  for(i in names(control)){
#      if(!is.null(control[[i]])){
#        formula = c(formula, i,"=control$",i,",")
#    }
#  }
#       formula <- c(formula,")")
#       formula <- paste0(formula, collapse = "")
### further control
	method <- match.arg(method)
	type <- match.arg(type)
## TODO think about CV default setting
	if(method=="CV"){	
		if(is.null(subsets)){
			warning("the subsets argument is not set: the results will be different each time, assuming 5 fold-CV")
			nfolds<- 5
			subsets <- lapply(as.data.frame(t(sapply(sample(rep_len(1:nfolds, length.out= dim(X)[1]), replace=FALSE) ,"!=", 1:nfolds))), which) ## apply on data frame not matrix to be sure to get a list output...
		} 
		tmp<- substr(type,0,1)[1] 
		if( tmp == "A" | tmp == "a" ) type = "agg" else type = "sel" ## default IC is selection, IC+agg is also known as AGGHOO
	} else { ## assuming IC :
		method ="IC"
		if(is.null(subsets)){
			subsets<- list(1:dim(X)[1]) #assume standard batch
		}
		if( is.character(ICpen)|is.null(ICpen) ){
			tmp<- substr(ICpen,0,1)[1]
			if( tmp == "A" | tmp == "a"){#AIC
				ICpen=2
			} else if (tmp == "H" | tmp == "h") { #HQC
				ICpen=2*log(log(nobs))
			} else {# default is BIC
				ICpen=log(nobs)
			}
		} else if (ICpen<=0) {
			ICpen<- log(nobs) ## set to BIC
		}
		tmp<- substr(type,0,1)[1] 
		if( tmp == "S" | tmp == "s" ) type = "sel" else type = "agg" ## default IC is aggregation
	}#method=="IC"

  #-------------------------------------------------
     xvar <- rep(0,   nobs) 
  #-------------------------------------------------
#      attr(xvar,"formula") = formula
       attr(xvar,"design") = X
      attr(xvar,"control") = control
        attr(xvar, "data") = as.data.frame(Data)
  attr(xvar, "gamlss.env") = gamlss.env
#        attr(xvar, "data") = as.data.frame(Data)
        attr(xvar, "call") = substitute(gamlss.gnet(data[[scall]], z, w, ...)) 
attr(xvar, "lambda") = lambda
attr(xvar, "method") = method
attr(xvar, "type") = type
attr(xvar, "ICpen") = ICpen
attr(xvar, "CVp") = CVp
attr(xvar, "k.se") = k.se
attr(xvar, "adaptive") = adaptive
attr(xvar, "epsilon") = epsilon
attr(xvar, "subsets") = subsets
attr(xvar, "sparse") = sparse
       attr(xvar, "class") = "smooth"
                 xvar
}
gnet.control = function(family="gaussian",
													offset = NULL, 
                           alpha = 1, 
                         nlambda = 100,
                lambda.min.ratio = 1e-3,  # p<n default setting
#                          lambda = NULL,
                     standardize = TRUE, 
                       intercept = TRUE, 
                          thresh = 1e-07,  
                           dfmax = NULL,
                            pmax = NULL, 
                         exclude = NULL, 
                  penalty.factor = NULL,
                    lower.limits = -Inf, 
                    upper.limits = Inf, 
                           maxit = 100000,
                   type.gaussian = NULL,
                   type.logistic = "Newton"
#										relax = FALSE ##should not be relevant, not relevant, we do it manually for the (selected) solution
)

{
list(family=family, offset=offset, alpha = alpha, nlambda = nlambda,lambda.min.ratio = lambda.min.ratio, standardize=standardize, intercept=intercept, thresh = thresh,  dfmax = dfmax, pmax = pmax, exclude = exclude, penalty.factor = penalty.factor,lower.limits=lower.limits, upper.limits=upper.limits, maxit=maxit, type.gaussian=type.gaussian,type.logistic=type.logistic)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
gamlss.gnet <-function(x, y, w, xeval = NULL, ...) {
  if (is.null(xeval))
  {#fitting
     control <- as.list(attr(x, "control"))
 #    X <- attr(x,"design")
     lambda <- attr(x,"lambda") 
     ICpen <- attr(x,"ICpen") 
     CVp <- attr(x,"CVp") 
     k.se <- attr(x,"k.se") 
     adaptive <- attr(x,"adaptive") 
     epsilon <- attr(x,"epsilon") 
     subsets <- attr(x,"subsets") 
     method <- attr(x,"method") 
     type <- attr(x,"type") 
     sparse <- attr(x,"sparse") 
	
	lmod<- list()
	nfolds<- length(subsets)
  nobs <- dim(attr(x,"design"))[1]
	nvars <- dim(attr(x,"design"))[2]
	if(is.null(control$dfmax) ) control$dfmax<- nvars+1
	if(is.null(control$pmax) ) control$pmax<- min(control$dfmax * 2 + 20,   nvars)
	if(is.null(control$penalty.factor) ) control$penalty.factor<- rep.int(1, nvars)
	if(is.null(control$type.gaussian) ) control$type.gaussian<- ifelse(nvars < 500, "covariance", "naive") 

	#CODE FOR LAMBDA COMPUTATION REFINED AND UNCOMMENTED BY JULIAN AMON ON 23rd July, 2021:
	#Problem was: so far when doing tuning via CV, lambda sequence would be different across
	#all folds, thereby not allowing for a correct choice of lambda via CV (in fact, output
	#was incorrect this way). Now: lambda sequence is chosen jointly for all folds and then
	#handed to glmnet to compute the lasso path in each fold:
  if(is.null(lambda)) {
    if(sparse){
      lambdamax <- max(abs(Matrix::crossprod(attr(x,"design"),y)))/nobs
    } else{
      lambdamax <- max(abs(crossprod(attr(x,"design"),y)))/nobs
    }
  	rangefactor <- 1.25 ## small extension for subset robustness
  	lambda <- exp(seq(log(lambdamax*rangefactor), log(lambdamax*control$lambda.min.ratio/rangefactor),
  	                  length.out=control$nlambda))
  	DF<- matrix(,control$nlambda, nfolds) ## IC for method IC, MSE for CV...
  } else {
    DF<- matrix(,length(lambda), nfolds) ## IC for method IC, MSE for CV...
  }
	
	OPTCRIT<- DF # same size

	pf<- control$penalty.factor
	for(i.adaptive in 1:(as.numeric(!is.null(adaptive)) +1) ){
		for(i.nets in 1:nfolds){
			if(i.adaptive>1) pf<- (abs(lmod[[i.nets]]$beta[, optid])+epsilon)^(-adaptive)
			lmod[[i.nets]] <- 		glmnet(x = attr(x,"design")[subsets[[i.nets]],], y = y[subsets[[i.nets]]], weights = w[subsets[[i.nets]]], lambda=lambda,	family = control$family,
														offset = control$offset, nlambda = control$nlambda, lambda.min.ratio = control$lambda.min.ratio,
		                   standardize = control$standardize,
		                     intercept = control$intercept, 
		                        thresh = control$thresh,
		                         dfmax = control$dfmax,
		                          pmax = control$pmax,
		                       exclude = control$exclude,
		                penalty.factor = pf,
		                  lower.limits = control$lower.limits, 
		                  upper.limits = control$upper.limits, 
		                         maxit = control$maxit,
		                 type.gaussian = control$type.gaussian,
		                 type.logistic = control$type.logistic)
			DF[1:lmod[[i.nets]]$dim[2],i.nets]<- lmod[[i.nets]]$df
			if(method=="IC"){
				RSS<- (1-lmod[[i.nets]]$dev.ratio)* lmod[[i.nets]]$nulldev
				OPTCRIT[1:lmod[[i.nets]]$dim[2],i.nets]<- log(RSS) + ICpen*lmod[[i.nets]]$df/nobs 
			}
			if(method=="CV"){
				res<- y[-subsets[[i.nets]]] - t(t(as.matrix(attr(x,"design")[-subsets[[i.nets]],]) %*% as.matrix(lmod[[i.nets]]$beta)) + lmod[[i.nets]]$a0)
	#			print(res)
				OPTCRIT[1:lmod[[i.nets]]$dim[2],i.nets]<- apply(abs(res)^CVp,2, mean)
			}
		}#i.nets
		## AGGHOO + ICagg 
		DFmean<- apply(DF,1,mean, na.rm=TRUE)
		optcritmean<- apply(OPTCRIT,1,mean)
		if(nfolds>1) optcritsd<- apply(OPTCRIT,1,sd) else optcritsd<- 0
		#NEW DEFAULT INTRODUCED BY JULIAN AMON ON 23rd July, 2021:
		#Confirmed: corresponds exactly to the 1-SE rule in the glmnet package:
		optid <- which(optcritmean <= optcritmean[which.min(optcritmean)] + k.se*optcritsd[which.min(optcritmean)])[1L]
		##OLD DEFAULT FROM CRAN PACKAGE gamlss.lasso AS FALLBACK INCASE THE ABOVE DOES NOT WORK:
		if(is.na(optid)) optid<-which.min(optcritmean+ k.se*optcritsd)
		
		
	}#i.adaptive
	fv<- numeric(nobs)


	BETA<- matrix(,nvars,nfolds)
	A0<- numeric(nfolds)
	for(i.nets in 1:nfolds){
		BETA[,i.nets]<- lmod[[i.nets]]$beta[, optid]
		A0[i.nets]<- lmod[[i.nets]]$a0[optid]
	}
	BETAsel<- numeric(nvars)
	A0sel<- numeric(1)
## another fit on full model

	if(type=="sel"){# CV + IC select
		pf<- control$penalty.factor
		for(i.adaptive in 1:(as.numeric(!is.null(adaptive)) +1) ){
			subsets[[length(subsets)+1]] <- 1:nobs
			i.nets<- nfolds+1
			if(i.adaptive>1) pf<- (abs(lmod[[i.nets]]$beta[, optid])+epsilon)^(-adaptive)
				lmod[[i.nets]] <- 		glmnet(x = attr(x,"design")[subsets[[i.nets]],], y = y[subsets[[i.nets]]], weights = w[subsets[[i.nets]]], lambda=lambda,	family = control$family,
														offset = control$offset, nlambda = control$nlambda, lambda.min.ratio = control$lambda.min.ratio,
		                   standardize = control$standardize,
		                     intercept = control$intercept, 
		                        thresh = control$thresh,
		                         dfmax = control$dfmax,
		                          pmax = control$pmax,
		                       exclude = control$exclude,
		                penalty.factor = pf,
		                  lower.limits = control$lower.limits, 
		                  upper.limits = control$upper.limits, 
		                         maxit = control$maxit,
		                 type.gaussian = control$type.gaussian,
		                 type.logistic = control$type.logistic)
			fv<- lmod[[i.nets]]$a0[optid] + attr(x,"design") %*% lmod[[i.nets]]$beta[,optid]
			BETAsel<- lmod[[i.nets]]$beta[,optid]
			A0sel<- lmod[[i.nets]]$a0[optid]
			dffin<- lmod[[i.nets]]$df[optid]
		}#adaptive
	} else {
	## AGGHOO + IC agg 
#print(",12")
		for(i.nets in 1:nfolds){
			BETAsel<- BETAsel + lmod[[i.nets]]$beta[,optid]/nfolds ## assuming equallty weighted folds, TODO could be generalised thus weighting based on subset size
			A0sel<- A0sel + lmod[[i.nets]]$a0[optid]/nfolds
		}
		fv<- A0sel + attr(x,"design") %*% BETAsel
		if(all(fv==0)) fv= rep.int(1e-15,nobs) ## FZ: numerical issue... don't ask me why... it is hidden add.fit
		dffin<- DFmean[optid]
	}
#	print(DF[optid,])
	## last entry
	if(is.null(lambda)) lambda <- lmod[[length(lmod)]]$lambda # note 
	lmod[[length(lmod)+1]]<- list(lambda=lambda, DF= DF, OPTCRIT=OPTCRIT,  optid=optid,  A0=A0, BETA=BETA, df=dffin, optcrit=optcritmean[optid], a0=A0sel, beta=BETAsel, optlambda= lambda[optid], scaled.beta = apply( attr(x,"design"),2,sd )*BETAsel)

#print(",")
#	Xmean <- apply( attr(x,"design"),1,mean )
#	Xsd <- 

		residuals <- y - fv
    list(fitted.values = fv, 
             residuals = residuals,
                 nl.df = dffin , 
                lambda = optid ,#list(list(optid=optid, lambdaopt=lmod$lambda[optid], ICpen=ICpen)), 	#FZ: small abuse, store everything relevant, esp. return also optindex
               coefSmo = lmod, 
                   var = NA)   # TODO think about computing variance  
  } else { # predict 
    gamlss.env <- as.environment(attr(x, "gamlss.env"))
           obj <- get("object", envir=gamlss.env ) # get the object from predict
            TT <- get("TT", envir=gamlss.env ) # get wich position is now
            SL <- get("smooth.labels", envir=gamlss.env) # all the labels of the smoother
           fit <- eval(parse(text=paste("obj$", get("what", envir=gamlss.env), ".coefSmo[[",as.character(match(TT,SL)), "]]", sep="")))
         OData <- attr(x,"data") 
            ll <- dim(OData)
            MM <- as.matrix(OData[seq(length(y)+1,ll[1]),-(ll[2])])

				if(dim(MM)[2] != length(fit[[length(fit)]]$beta) ) stop("dimension mismatch")
	b<-fit[[length(fit)]]$beta
	a<-fit[[length(fit)]]$a0

 #         pred <- a-sum(b*scale_mu/scale_sd)+as.numeric(crossprod(t(MM)/scale_sd, b)) ##rescaling
         pred <- a-sum(b)+as.numeric(crossprod(t(MM), b)) ##rescaling
  
   pred
  }         
}







