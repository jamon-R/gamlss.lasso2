\name{gamlss.lrs}
\Rdversion{1.1}
\alias{gamlss.lrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Support for Function lrs() 
}
\description{
This is support for the smoother function lrs() an interface for Brad Efron and Trevor Hastie for \code{lars()} function.
It is not intended to be called directly by users. 
}
\usage{
gamlss.lrs(x, y, w, xeval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the explanatory variables}
  \item{y}{iterative y variable}
  \item{w}{iterative weights}
  \item{xeval}{if xeval=TRUE then predicion is used}
  \item{\dots}{for extra arguments}
}

\value{No return value, called for GAMLSS lrs procedure.} 

\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Ripley, B. D. (1996) \emph{Pattern Recognition and Neural Networks}. Cambridge. 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{https://www.gamlss.com/}).

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer. 


}
\author{
Florian Ziel, Peru Muniain and Mikis Stasinopoulos
}


\seealso{
\code{\link{lrs}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
