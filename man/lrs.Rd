\name{lrs}
\alias{lrs}

\title{
Least angle regression and lasso in GAMLSS
}
\description{
This function allows estimating the different components of a GAMLSS model (mean, sd. dev., skewness and kurtosis) using the elastic net (with lasso as default special case) estimation method via glmnet.  This method is appropriate for models with many variables.
}
\usage{
lrs(X = NULL, x.vars = NULL, lambda = NULL,	method = c("IC","CV"), 
  type = c("agg","sel"), ICpen = c("BIC", "HQC", "AIC"), CVp = 2, k.se = 0, 
  subsets = NULL, lars.type= "lasso", use.gram = TRUE, 
  eps = .Machine$double.eps, max.steps = NULL, ...)
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
  \item{subsets}{ The subsets for cross-validation, information criteria or bootstraping, by default 5 random fold are selected.}
  \item{lars.type}{As in \code{lars}, lars type, e.g. "lasso", "lar" (least angle regression), "forward.stagewise" or "stepwise".}
 	\item{use.gram}{States if Gramian should be precomputed, default TRUE - recommended as gamlss will call lars often during the estimation.}
 	\item{eps}{As in \code{lars}, a small constant.}
 	\item{max.steps}{As in \code{lars}, number of updating steps (for "lars" method equal to number of variables, for "lasso" it can be smaller), default NULL.}
  \item{\dots}{for extra arguments}
}
\details{
The estimation of the lambda is carried out by BIC by default.
If the objective is to predict the model must be defined by \code{x.vars}. 
Different types of subsets must be constructed if bootstrapping and aggregation are applied, as in this case observations might be repeated. 
% TO DO: decide whether we want a comment on alpha value for the elastic net or not.
}
\value{
This function returns a smooth object of the GAMLSS model. It contains the estimated parameters and related characteristics for the \code{lars} component in the GAMLSS model we are estimating.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale, and shape,(with discussion), \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F. (2019) Distributions for modeling location, scale, and shape: Using GAMLSS in R, Chapman and Hall/CRC. An older version can be found in https://www.gamlss.com/. 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, \bold{Vol. 23}, Issue 7, Dec 2007, https://www.jstatsoft.org/v23/i07/. 

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) Flexible Regression and Smoothing: Using GAMLSS in R, Chapman and Hall/CRC. 

Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, \emph{Journal of Statistical Software}, \bold{Vol. 39(5)}, 1-13, https://www.jstatsoft.org/v39/i05/.

 Tibshirani, Robert, Bien, J., Friedman, J., Hastie, T.,Simon, N.,Taylor, J. and Tibshirani, Ryan. (2012) Strong Rules for Discarding Predictors in Lasso-type Problems, \emph{JRSSB}, \bold{Vol. 74(2)}, 245-266, https://statweb.stanford.edu/~tibs/ftp/strong.pdf.
 
 Hastie, T., Tibshirani, Robert and Tibshirani, Ryan. Extended Comparisons of Best Subset Selection, Forward Stepwise Selection, and the Lasso (2017), \emph{Stanford Statistics Technical Report}, https://arxiv.org/abs/1707.08692.

 Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}.




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
data<- cbind(y=as.numeric(rnorm(n,sd=ysd)) + t(tcrossprod( BETA[,1],X)),as.data.frame(X))

# Estimating the model with lrs default setting
mod <- gamlss(y~lrs(x.vars=names(data)[-1] ),
              sigma.fo=~lrs(x.vars=names(data)[-1]), data=data, family=NO,
              i.control = glim.control(cyc=1, bf.cyc=1))

# Estimated paramters are available at
rbind(true=BETA[,1],estimate=tail(getSmo(mod, "mu") ,1)[[1]]$beta )## beta for mu
rbind(true=BETA[,2],estimate=tail(getSmo(mod, "sigma") ,1)[[1]]$beta )## beta for sigma

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
