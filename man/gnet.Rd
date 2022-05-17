\name{gnet}
\alias{gnet}
\alias{gnet.control}

\title{
(Adaptive) elastic net in GAMLSS
}
\description{
This function allows estimating the different components of a GAMLSS model (location, shape, scale parameters) using the (adaptive) elastic net (with adaptive lasso as default special case) estimation method via \code{glmnet}. This method is appropriate for models with many variables.
}
\usage{
gnet( X = NULL, x.vars = NULL, lambda = NULL, method = c("IC", "CV"), 
  type = c("agg", "sel"), ICpen = c("BIC", "HQC", "AIC"), CVp = 2, 
  k.se = 0, adaptive = 1, epsilon = 1/sqrt(dim(X)[1]), subsets = NULL, 
  sparse = FALSE, control = gnet.control(...), ...) 
gnet.control(family="gaussian", offset = NULL, alpha = 1, nlambda = 100, 
  lambda.min.ratio = 1e-3, standardize = TRUE, intercept = TRUE, thresh = 1e-07,
  dfmax = NULL, pmax = NULL, exclude = NULL, penalty.factor = NULL, lower.limits = -Inf,
  upper.limits = Inf, maxit = 100000, type.gaussian = NULL, type.logistic = "Newton")
}
\arguments{
  \item{X}{ The data frame containing the explanatory variables.}
  \item{x.vars}{ Indicates the name of the variables that must be included as explanatory variables from data the data object of GAMLSS. The explanatory variables must be included by \code{X} or by \code{x.vars}.}
  \item{lambda}{ The provided lambda grid. By default \code{NULL}.}
  \item{method}{ The method used to calculate the optimal lambda. If \code{method="IC"} information criteria are used, the penalization for the information criterion is selected in \code{ICpen}.If \code{method="CV"} cross validation resp. sampling is used, the penalization for the cross-validation is selected in \code{CVp}.}
  \item{type}{ The way to select the optimal lambda across the subsample fits. If \code{type="sel"} the optimal lambda is computed by selection. If \code{method="agg"} the optimal lambda is computed by aggregation.}
  \item{ICpen}{ The penalization for the information criteria. If \code{ICpen="AIC"} or  \code{ICpen=2} the optimal lambda is computed by Akaike Information Criterion. If \code{ICpen="BIC"} the optimal lambda is computed by Bayesian Information Criterion.If \code{ICpen="HQC"} the optimal lambda is computed by Hannan-Quinn Information Criterion.}
  \item{CVp}{ The penalization for the cross-validation, establishes the power of the error term. By default is equal to 2, i.e. squared error.}
  \item{k.se}{ This parameter establishes how many times the standard deviation is summed to the mean to select the optimal lambda. By default is equal to 0.}
  \item{adaptive}{ This parameter specifies if adaptive lasso shall be used, the default is 1. If \code{NULL} then standard lasso is used, otherwise adaptive lasso with penalty weights \code{(abs(beta)+epsilon)^(-adaptive)} where \code{beta} is chosen from an initial standard lasso estimate and \code{epsilon} is specified by the next parameter. Note, estimating standard lasso requires about half of computation time, but adaptive lasso has smaller bias and satiesfies the oracle property.}
  \item{epsilon}{ This parameter specifies the adaptive lasso penalty weights. The default is \code{1/sqrt(dim(X)[1])}. }
  \item{subsets}{ The subsets for cross-validation, information criteria or bootstraping, by default 5 random folds are selected.}
  \item{sparse}{ If sparse converts input matrix for \code{glmnet} into a sparse Matrix, may reduces computation time for sparse designs.}
  \item{control}{List of further input parameters for glmnet, e.g. alpha for elastic net parameters.}
  \item{\dots}{for extra arguments}
\item{family}{Either a character string representing
one of the built-in families, or else a \code{glm()} family object.}
\item{offset}{A vector of length \code{nobs} that is included in the linear
predictor (a \code{nobs x nc} matrix for the \code{"multinomial"} family).
Useful for the \code{"poisson"} family (e.g. log of exposure time), or for
refining a model by starting at a current fit. Default is \code{NULL}. If
supplied, then values must also be supplied to the \code{predict} function.}
\item{alpha}{The elastic net mixing parameter, with \eqn{0\le\alpha\le 1}.
The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} is the
lasso penalty, and \code{alpha=0} the ridge penalty. Default is lasso.}
\item{nlambda}{Size of the tuning parameter grid, default is 100. It is irrelevant if lambda is explicitly specified.}
	\item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
\code{lambda.max}, the (data derived) entry value (i.e. the smallest value
for which all coefficients are zero). The default is \code{0.001}.  A very small value of
\code{lambda.min.ratio} will lead to a saturated fit in the \code{nobs <
nvars} case. This is undefined for \code{"binomial"} and
\code{"multinomial"} models, and \code{glmnet} will exit gracefully when the
percentage deviance explained is almost 1. It is irrelevant if lambda is explicitly specified.}
	\item{standardize}{Logical flag for \code{X} or \code{x.vars} variable standardization, prior to
fitting the model sequence. The coefficients are always returned on the
original scale. Default is \code{standardize=TRUE} and it is highly recommended.}
	\item{intercept}{Should intercept(s) be fitted (default=TRUE) or set to zero
(FALSE).}
	\item{thresh}{Convergence threshold for coordinate descent. Each inner
coordinate-descent loop continues until the maximum change in the objective
after any coefficient update is less than \code{thresh} times the null
deviance. Defaults value is \code{1E-7}.}
	\item{dfmax}{Limit the maximum number of variables in the model. Useful for
very large \code{nvars}, if a partial path is desired.}
	\item{pmax}{Limit the maximum number of variables ever to be nonzero.}
	\item{exclude}{Indices of variables to be excluded from the model. Default
is none. Equivalent to an infinite penalty factor (next item).}
	\item{penalty.factor}{Separate penalty factors can be applied to each
coefficient. This is a number that multiplies \code{lambda} to allow
differential shrinkage. Can be 0 for some variables, which implies no
shrinkage, and that variable is always included in the model. Default is 1
for all variables (and implicitly infinity for variables listed in
\code{exclude}). Note: the penalty factors are internally rescaled to sum to
nvars, and the lambda sequence will reflect this change.}
	\item{lower.limits}{Vector of lower limits for each coefficient; default
\code{-Inf}. Each of these must be non-positive. Can be presented as a
single value (which will then be replicated), else a vector of length
\code{nvars}.}
	\item{upper.limits}{Vector of upper limits for each coefficient; default
\code{Inf}. See \code{lower.limits}.}
	\item{maxit}{Maximum number of passes over the data for all lambda values;
default is 10^5.}
	\item{type.gaussian}{Two algorithm types are supported for (only)
\code{family="gaussian"}. The default when \code{nvar<500} is
\code{type.gaussian="covariance"}, and saves all inner-products ever
computed. This can be much faster than \code{type.gaussian="naive"}, which
loops through \code{nobs} every time an inner-product is computed. The
latter can be far more efficient for \code{nvar >> nobs} situations, or when
\code{nvar > 500}.}
	\item{type.logistic}{If \code{"Newton"} then the exact hessian is used
(default), while \code{"modified.Newton"} uses an upper-bound on the
hessian, and can be faster.}
}
\details{
The estimation of the lambda is carried out by BIC by default.
If the objective is to predict the model must be defined by \code{x.vars}. 
Different types of subsets must be constructed if bootstrapping and aggregation are applied, as in this case observations might be repeated. 
% TO DO: decide whether we want a comment on alpha value for the elastic net or not.
}
\value{
This function returns a smooth object of the GAMLSS model. It contains the estimated parameters and related characteristics for the \code{glmnet} component in the GAMLSS model we are estimating.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale, and shape,(with discussion), \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC. An older version can be found in https://www.gamlss.com/. 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, \bold{Vol. 23}, Issue 7, Dec 2007, https://www.jstatsoft.org/v23/i07/. 

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) Flexible Regression and Smoothing: Using GAMLSS in R, Chapman and Hall/CRC. 

Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, \emph{Journal of Statistical Software}, \bold{Vol. 39(5)}, 1-13, https://www.jstatsoft.org/v39/i05/.

 Tibshirani, Robert, Bien, J., Friedman, J., Hastie, T.,Simon, N.,Taylor, J. and Tibshirani, Ryan. (2012) Strong Rules for Discarding Predictors in Lasso-type Problems, \emph{JRSSB}, \bold{Vol. 74(2)}, 245-266, https://statweb.stanford.edu/~tibs/ftp/strong.pdf.
 
 Hastie, T., Tibshirani, Robert and Tibshirani, Ryan. Extended Comparisons of Best Subset Selection, Forward Stepwise Selection, and the Lasso (2017), \emph{Stanford Statistics Technical Report}, https://arxiv.org/abs/1707.08692.




}
\author{
Florian Ziel, Peru Muniain and Mikis Stasinopoulos
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# Contructing the data
library(gamlss.lasso2)
set.seed(123)
n<- 500
d<- 50
X<- matrix(rnorm(n*d), n,d)
BETA<- cbind( "mu"=rbinom(d,1,.1), "sigma"= rbinom(d,1,.1)*.3)
ysd<- exp(1 + tcrossprod( BETA[,2],X))
data<- cbind(y=as.numeric(rnorm(n, sd=ysd))+t(tcrossprod( BETA[,1],X)), as.data.frame(X))


# Estimating the model using default setting: adaptive lasso with BIC tuning
mod <- gamlss(y~gnet(x.vars=names(data)[-1]),
              sigma.fo=~gnet(x.vars=names(data)[-1]), data=data,
              family=NO, i.control = glim.control(cyc=1, bf.cyc=1))

# Estimating the model with standard lasso (BIC tuning)
mod.lasso <- gamlss(y~gnet(x.vars=names(data)[-1], adaptive=NULL),
              sigma.fo=~gnet(x.vars=names(data)[-1], adaptive=NULL), data=data, 
              family=NO, i.control = glim.control(cyc=1, bf.cyc=1))

# Estimated paramters are available at
rbind(true=BETA[,1],alasso=tail(getSmo(mod, "mu") ,1)[[1]]$beta,
                    lasso=tail(getSmo(mod.lasso, "mu") ,1)[[1]]$beta) ##beta for mu
rbind(true=BETA[,2],alasso=tail(getSmo(mod, "sigma") ,1)[[1]]$beta,
                    lasso=tail(getSmo(mod.lasso, "sigma") ,1)[[1]]$beta)##beta for sigma

# Estimating with other setting
nfolds<- 6
n<- dim(data)[1]
# folds for cross-validation and bootstrap
CVfolds<- lapply(as.data.frame(t(sapply(sample(rep_len(1:nfolds,length=n),replace=FALSE)
                 ,"!=", 1:nfolds))), which)  
BOOTfolds<- lapply(as.data.frame(matrix(sample(1:n, nfolds*n, replace=TRUE), n)),sort)  

#Bootstrap + Aggrationg = Bagging:
mod1 <- gamlss(y~gnet(x.vars=names(data)[-1], method="CV",type="agg", subsets=BOOTfolds),
               sigma.fo=~gnet(x.vars=names(data)[-1]), data=data, family=NO,
               i.control = glim.control(cyc=1, bf.cyc=1)) 

# Estimated paramters are available at
tail(getSmo(mod1, "mu") ,1)[[1]]$beta ## beta for mu
tail(getSmo(mod1, "sigma") ,1)[[1]]$beta ## beta for sigma

# Cross-validation (with selection):
mod2 <- gamlss(y~gnet(x.vars=names(data)[-1],method="CV",type="sel", subsets=CVfolds),
               sigma.fo=~gnet(x.vars=names(data)[-1],method="CV",type="sel", ICpen=2, 
               subsets=CVfolds), data=data, family=NO,
               i.control = glim.control(cyc=1, bf.cyc=1)) 

# Estimated paramters are available at
tail(getSmo(mod2, "mu") ,1)[[1]]$beta ## beta for mu
tail(getSmo(mod2, "sigma") ,1)[[1]]$beta ## beta for sigma


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
