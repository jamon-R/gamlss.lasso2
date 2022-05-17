\name{gamlss.lasso2-package}
\alias{gamlss.lasso2-package}
\alias{gamlss.lasso2}
\docType{package}
\title{
\packageTitle{gamlss.lasso2}
}
\description{
\packageDescription{gamlss.lasso2}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{gamlss.lasso2}
\packageIndices{gamlss.lasso2}

}
\author{
\packageAuthor{gamlss.lasso2}

Maintainer: \packageMaintainer{gamlss.lasso2}
}
\references{
R
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{https://www.gamlss.com/}).

Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). \emph{Least angle regression}. Annals of statistics, 32(2), 407-499.

Friedman, J., Hastie, T., & Tibshirani, R. (2010). \emph{Regularization paths for generalized linear models via coordinate descent}. Journal of statistical software, 33(1), 1.

}

\keyword{package}
\seealso{\code{\link[gamlss]{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}}, \code{\link[gamlss.add]{gamlss.add}}}  
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

# Estimating the model with gnet default setting
mod <- gamlss(y~gnet(x.vars=names(data)[-1] ),
              sigma.fo=~gnet(x.vars=names(data)[-1]), data=data, family=NO,
              i.control = glim.control(cyc=1, bf.cyc=1))

# Estimated paramters are available at
rbind(true=BETA[,1],estimate=tail(getSmo(mod, "mu") ,1)[[1]]$beta )## beta for mu
rbind(true=BETA[,2],estimate=tail(getSmo(mod, "sigma") ,1)[[1]]$beta )## beta for sigma
}
