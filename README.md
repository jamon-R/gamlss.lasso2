
# gamlss.lasso2

## Overview

**`gamlss.lasso2`** is an R package based on the original `gamlss.lasso`
package available on [CRAN](https://cran.r-project.org/). Compared to
the original package, two major changes have been implemented in the
function `gamlss.gnet`: - Tuning via cross-validation (CV): The problem
was that so far when doing tuning via CV, the lambda sequence would be
different across all folds, thereby not allowing for a correct choice of
lambda via CV (in fact, the output was incorrect this way). Now: the
lambda sequence is chosen jointly for all folds and then handed to
glmnet to compute the lasso path in each fold. - 1-SE rule when tuning
via CV: So far, `gamlss.gnet` would choose the lambda for which the
cross-validated error plus `k.se` time the standard error was minimized.
However, this does not correspond to the standard in the `glmnet`
package, where the lambda is chosen to be the smallest value such that
the cross-validated standard error is still within a `k.se`-standard
error band around the minimum. This procedure was implemented also here.

## Installation

The source code is hosted on
[github](https://github.com/jamon-R/gamlss.lasso2). To install the *R*
package, run the following code:

``` r
devtools::install_github("jamon-R/gamlss.lasso2")
```

The package requires functionality from the packages `gamlss`, `glmnet`,
`lars` and `Matrix`. All of those dependencies are available on
[CRAN](https://cran.r-project.org/).
