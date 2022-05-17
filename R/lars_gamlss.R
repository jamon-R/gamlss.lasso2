# fit glmnet terms using the glmnet function of glmnet package
# which is used in the backfitting 
# TO DO
# (i)  is scalling the X always needed it?
# (ii)  at the moment X is a matrix should we allow formula?
# (iii)
#--------------------------------------------------------------------------------------
require(lars)
lrs <- function(  ## this is lars
                 X = NULL, 
            x.vars = NULL,
						lambda = NULL,	
						method = c("IC", "CV"),
						type = c("agg", "sel"),
						ICpen = c("BIC", "HQC", "AIC"),
						CVp = 2,
						k.se = 0,
            subsets = NULL,
						lars.type = "lasso",	
						use.gram = TRUE,	
						eps = .Machine$double.eps, 
						max.steps = NULL, ...) {
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
  # prepare glmnet

#	if (sys.call(position)[1]!="predict.gamlss()"){ # not predict.gamlss
	if (is.null(X)&&is.null(x.vars)) stop("X or x.vars has to be set in gnet")
	if (is.null(X)&&!is.null(x.vars)) X <- as.matrix(Data[, x.vars])
  if (!is.null(X)&&is.null(x.vars)) warning("For prediction use the x.vars argument")

#	}
   X <- scale(X) #FZ: could be done in glmnet, but we do it here to save time, No stand. but would only make sense for "mu", BUT this has to be taken into account for pred
  nobs <- dim(X)[1]
  nvars <- dim(X)[2]
	## TODO weights are not in lars, but can be done manually...
#	print("asdasda") 


## TODO check if subsets has correct format ## think about CV as default
	if(method=="CV"){	
		if(is.null(subsets)){
			warning("the subsets argument is not set: the results will be different each time, assuming 5 folds")
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



	if(use.gram){
		gram<- list()
		for(i in 1:length(subsets))	gram[[i]]<- crossprod(X[subsets[[i]],]) ## TODO could be optimized... but should not be the major problem, as for very large problems, lars is not suitable anyway...
		if(type=="sel") gram[[length(subsets)+1]]<- crossprod(X) # for refit
	} else {
		gram<- NULL
	}

	if(is.null(max.steps)) max.steps <- 8 * min(nvars, nobs) # default in lars 1.2

  #-------------------------------------------------
     xvar <- rep(0,   nobs) 
  #-------------------------------------------------
#      attr(xvar,"formula") = formula
       attr(xvar,"design") = X
#      attr(xvar,"control") = control
        attr(xvar, "data") = as.data.frame(Data)
  attr(xvar, "gamlss.env") = gamlss.env
#        attr(xvar, "data") = as.data.frame(Data)
        attr(xvar, "call") = substitute(gamlss.lrs(data[[scall]], z, w, ...)) 
attr(xvar, "lambda") = lambda
attr(xvar, "lars.type") = lars.type
attr(xvar, "gram") = gram
attr(xvar, "method") = method
attr(xvar, "type") = type
attr(xvar, "ICpen") = ICpen
attr(xvar, "CVp") = CVp
attr(xvar, "k.se") = k.se
attr(xvar, "subsets") = subsets
attr(xvar, "max.steps") = max.steps
attr(xvar, "eps") = eps
       attr(xvar, "class") = "smooth"
                 xvar
}
# lars(x, y, type = c("lasso", "lar", "forward.stagewise", "stepwise"), 
##         trace = FALSE, normalize = TRUE, intercept = TRUE, Gram, eps = .Machine$double.eps, max.steps, use.Gram = TRUE)
#lars.control = function(type="lasso",
#													trace = FALSE, 
#                           normalize = TRUE, 
#                         intercept = TRUE,
#                          eps = .Machine$double.eps,
#                          ...)
#{
#list(type=type, trace = trace) ## TODO
#}
##--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
gamlss.lrs <-function(x, y, w, xeval = NULL, ...) {
  if (is.null(xeval))
  {#fitting
#	print("fssfs") 
    control <- as.list(attr(x, "control"))
           X <- as.matrix(attr(x,"design"))

#     formula <- attr(x,"formula") 
#print(attr(x,"lambda") )
     lambda <- attr(x,"lambda") 
     ICpen <- attr(x,"ICpen") 
     CVp <- attr(x,"CVp") 
     k.se <- attr(x,"k.se") 
     subsets <- attr(x,"subsets") 
     method <- attr(x,"method") 
     type <- attr(x,"type") 


	lmod<- list()
	nfolds<- length(subsets)
  nobs <- dim(X)[1]
	nvars <- dim(X)[2]
	if(is.null(lambda)) ncoef<- nvars +1  else ncoef<- length(lambda) ## for lars // lasso if only taking min rss per df

#	lmod<- lars(X,y, intercept = FALSE, normalize=FALSE, Gram=attr(x,"gram")) ## input should be normalized

	OPTCRIT<- matrix(,ncoef, nfolds) ## IC for method IC, MSE for CV...
	DF<- matrix(,ncoef, nfolds) ## IC for method IC, MSE for CV...
	IDSEL<- matrix(,ncoef, nfolds) ## IC for method IC, MSE for CV...

	lmod<- list()
	for(i.subset in 1:nfolds){
		lmod[[i.subset]] <- lars(X[subsets[[i.subset]],],y[subsets[[i.subset]]], intercept = FALSE, normalize=FALSE, type=attr(x,"lars.type")[[i.subset]], Gram=attr(x,"gram")[[i.subset]], eps=attr(x,"eps"), max.steps=attr(x,"max.steps"))#
		if(is.null(lambda)){ ## selection on path
			IDSEL[, i.subset]<- length(lmod[[i.subset]]$df)-match(0:nvars,rev(lmod[[i.subset]]$df))+1
			DF[,i.subset]<- lmod[[i.subset]]$df[IDSEL[, i.subset]] ## should be 0:nvars
			if(method=="IC"){
				OPTCRIT[,i.subset]<- log( lmod[[i.subset]]$RSS[IDSEL[,i.subset]] /nobs) + ICpen*lmod[[i.subset]]$df[IDSEL[, i.subset]]/nobs 
			}
			if(method=="CV"){
				fit<- predict(lmod[[i.subset]], newx=X[-subsets[[i.subset]],], s=IDSEL[, i.subset])$fit
				res<- y[-subsets[[i.subset]]] - fit#+ lmod[[i.subset]]$a0)
				OPTCRIT[,i.subset]<- apply(abs(res)^CVp,2, mean)
			}
		} else { ## selection on lambda grid
			fbeta<- matrix(predict(lmod[[i.subset]],mode="lambda", type="coef", s=lambda)$coef, length(lambda))
			DF[,i.subset]<- 	apply(fbeta!=0,1,sum)
			if(method=="IC"){
				fit<- as.matrix(predict(lmod[[i.subset]], newx=X[subsets[[i.subset]],],mode="lambda", s=lambda)$fit)
				res<- y[subsets[[i.subset]]] - fit
				RSS<- apply(res*res, 2, mean)
				OPTCRIT[,i.subset]<- log( RSS ) + ICpen*DF[,i.subset]/nobs 
			} else {
				fit<- as.matrix(predict(lmod[[i.subset]], newx=X[-subsets[[i.subset]],],mode="lambda", s=lambda)$fit)
				res<- y[-subsets[[i.subset]]] - fit
				OPTCRIT[,i.subset]<- apply(abs(res)^CVp,2, mean) #log( RSS ) + ICpen*fdf/nobs 
			}
#CV				
		}
	}#i.subset
	## AGGHOO + ICagg 
	DFmean<- apply(DF,1,mean, na.rm=TRUE)
	optcritmean<- apply(OPTCRIT,1,mean)
	if(nfolds>1) optcritsd<- apply(OPTCRIT,1,sd) else optcritsd<- 0
	optid<-which.min(optcritmean+ k.se*optcritsd) ## 
	fv<- numeric(nobs)

	BETA<- matrix(,nvars,nfolds)
	for(i.subset in 1:nfolds){
		if(is.null(lambda)){
			BETA[,i.subset]<- lmod[[i.subset]]$beta[optid, ]
		} else {
			BETA[,i.subset] <- matrix(predict(lmod[[i.subset]],mode="lambda", type="coef", s=lambda)$coef,length(lambda))[optid, ]
		}
	}
	BETAsel<- numeric(nvars)
# another fit on full model
	if(type=="sel"){# CV + IC select
		i.subset<- nfolds+1
		lmod[[i.subset]] <- lars(X,y, intercept = FALSE, normalize=FALSE, type=attr(x,"lars.type")[[i.subset]], Gram=attr(x,"gram")[[i.subset]], eps=attr(x,"eps"), max.steps=attr(x,"max.steps")) #
		if(is.null(lambda)){
			idsel <- length(lmod[[i.subset]]$df)-match(optid-1,rev(lmod[[i.subset]]$df))+1 ## optid == index of opt df -1
			fv<-  X %*% lmod[[i.subset]]$beta[idsel,] # lmod[[i.subset]]$a0[optid]
			BETAsel<- lmod[[i.subset]]$beta[idsel,]
			dffin<- lmod[[i.subset]]$df[idsel]
		} else {
			BETAsel<- predict(lmod[[i.subset]],mode="lambda", type="coef", s=lambda[optid])$coef
			fv<- predict(lmod[[i.subset]], newx=X,mode="lambda", s=lambda[optid])$fit
			dffin<- 	sum(BETAsel!=0)
		}
	} else {
	## AGGHOO + IC agg 
		for(i.subset in 1:nfolds){
			if(is.null(lambda)) b<- lmod[[i.subset]]$beta[optid,] else b<- predict(lmod[[i.subset]],mode="lambda", type="coef", s=lambda[optid])$coef
			BETAsel<- BETAsel + b/nfolds ## assuming equallty weighted folds, TODO could be generalised thus weighting based on subset size
		}
		fv<- X %*% BETAsel #
		dffin<- DFmean[optid]
	}
#	print(DF[optid,])
	## last entry
	lmod[[length(lmod)+1]]<- list(lambda=lambda, DF= DF, OPTCRIT=OPTCRIT,  optid=optid, BETA=BETA, df=dffin, optcrit=optcritmean[optid], beta=BETAsel, optlambda= lambda[optid])

		if(all(fv==0)) fv= rep.int(1e-15,nobs) ## FZ: numerical issue... don't ask me why... I think gnet should have the same
		residuals <- y - fv


    list(fitted.values = fv, 
             residuals = residuals,
                 nl.df = dffin , 
#	               lam2 = lmod$lambda[optid], 
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
#	print(str(MM))
#	print(str(fit[[length(fit)]]$beta))
#	glob<<- fit
#	glob2<<- MM
	scale_mu<- attr(attr(x,"design"),"scaled:center")
	scale_sd<- attr(attr(x,"design"),"scaled:scale")

				if(dim(MM)[2] != length(fit[[length(fit)]]$beta) ) stop("dimension mismatch")
	b<-fit[[length(fit)]]$beta

          pred <- -sum(b*scale_mu/scale_sd)+as.numeric(crossprod(t(MM)/scale_sd, b)) ##TODO still glmnet
    pred
  }         
}







