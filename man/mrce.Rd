\name{mrce}
\alias{mrce}
\title{Do multivariate regression with covariance estimation (MRCE)}
\description{
Let \eqn{S_{+}^q} be the set of \eqn{q} by \eqn{q} symmetric and positive definite matrices
and let \eqn{y_i\in R^q} be the measurements of the \eqn{q} responses for the \eqn{i}th subject
(\eqn{i=1,\ldots, n}).
The model assumes that \eqn{y_i}
is a realization of the \eqn{q}-variate random vector
\deqn{Y_i = \mu + \beta'x_i + \varepsilon_i, \ \ \ \ i=1,\ldots, n}  
where \eqn{\mu\in R^q} is an unknown intercept vector;
\eqn{\beta\in R^{p\times q}} is an unknown regression coefficient matrix;
\eqn{x_i \in R^p} is the known vector of values for \eqn{i}th subjects's predictors,
and \eqn{\varepsilon_1,\ldots, \varepsilon_n} are \eqn{n} independent copies of a \eqn{q}-variate Normal random 
vector with mean 0 and unknown inverse covariance 
matrix \eqn{\Omega \in S_{+}^q}.

This function computes penalized likelihood estimates of the unknown parameters 
\eqn{\mu}, \eqn{\beta}, and \eqn{\Omega}.
Let \eqn{\bar y=n^{-1} \sum_{i=1}^n y_i} and \eqn{\bar{x} = n^{-1}\sum_{i=1}^n x_i}.
These estimates are
\deqn{
(\hat{\beta}, \hat\Omega) = \arg\min_{(B, Q)\in R^{p\times q}\times S_{+}^q}
\left\{g(B, Q) +\lambda_1 
\left(\sum_{j\neq k} |Q_{jk}| + 1(p\geq n) \sum_{j=1}^q |Q_{jj}| \right)
+ 2\lambda_{2}\sum_{j=1}^p\sum_{k=1}^q  |B_{jk}|\right\}}
and \eqn{\hat\mu=\bar y - \hat\beta'\bar x}, where 
\deqn{
g(B, Q) = {\rm tr}\{n^{-1}(Y-XB)'(Y-XB) Q\}-\log|Q|,}
\eqn{Y\in R^{n\times q}} has \eqn{i}th row \eqn{(y_{i}-\bar y)'},
and \eqn{X\in R^{n\times p}} has \eqn{i}th row \eqn{(x_{i}-\bar{x})'}.
}
\usage{
mrce(X,Y, lam1=NULL, lam2=NULL, lam1.vec=NULL, lam2.vec=NULL,
     method=c("single", "cv", "fixed.omega"),
     cov.tol=1e-4, cov.maxit=1e3, omega=NULL, 
     maxit.out=1e3, maxit.in=1e3, tol.out=1e-8, 
     tol.in=1e-8, kfold=5, silent=TRUE, eps=1e-5, 
     standardize=FALSE, permute=FALSE)
}
\arguments{
  \item{X}{An \eqn{n} by \eqn{p} matrix of the values for the prediction variables. 
  The \eqn{i}th row of \code{X} is \eqn{x_i} defined above (\eqn{i=1,\ldots, n}). 
  Do not include a column of ones.}
  \item{Y}{An \eqn{n} by \eqn{q} matrix of the observed responses. 
  The \eqn{i}th row of \code{Y} is \eqn{y_i} defined above (\eqn{i=1,\ldots, n}).}
  \item{lam1}{A single value for \eqn{\lambda_1} defined above.  This
  argument is only used if \code{method="single"}}
  \item{lam2}{A single value for \eqn{\lambda_2} defined above
  (or a \eqn{p} by \eqn{q} matrix with \eqn{(j,k)}th entry \eqn{\lambda_{2jk}}
  in which case the penalty \eqn{2\lambda_{2}\sum_{j=1}^p\sum_{k=1}^q  |B_{jk}|} becomes
\eqn{2\sum_{j=1}^p\sum_{k=1}^q \lambda_{2jk}|B_{jk}|}). This
  argument is not used if \code{method="cv"}.}
  \item{lam1.vec}{A vector of candidate values for \eqn{\lambda_1} from which the cross validation procedure
  searches: only used when \code{method="cv"}
  and must be specified by the user when \code{method="cv"}.  Please arrange in decreasing order.}
  \item{lam2.vec}{A vector of candidate values for \eqn{\lambda_2} from which the cross validation procedure
  searches: only used when \code{method="cv"}
  and must be specified by the user when \code{method="cv"}. Please arrange in decreasing order.}
  \item{method}{
  There are three options:
  \itemize{
\item \code{method="single"} computes the MRCE estimate of the regression coefficient matrix
  with penalty tuning parameters \code{lam1} and \code{lam2};
\item \code{method="cv"} performs \code{kfold} cross
  validation using candidate tuning parameters in \code{lam1.vec} and \code{lam2.vec};
\item \code{method="fixed.omega"} computes the regression coefficient matrix estimate for which \eqn{Q} (defined above) 
  is fixed at \code{omega}.
}
}  
\item{cov.tol}{Convergence tolerance for the glasso algorithm that minimizes the objective function (defined above)
  with \eqn{B} fixed.}
  \item{cov.maxit}{The maximum number of iterations allowed for the glasso algorithm that minimizes the objective function 
  (defined above)
  with \eqn{B} fixed.}
  \item{omega}{A user-supplied fixed value of \eqn{Q}.  Only used when 
  \code{method="fixed.omega"} in which case the minimizer of the objective function (defined above) with \eqn{Q}
  fixed at \code{omega} is returned.  
  }
  \item{maxit.out}{The maximum number of iterations allowed for the outer loop of the exact MRCE algorithm.}
  \item{maxit.in}{The maximum number of iterations allowed for the algorithm that minimizes the objective function, defined above,
  with \eqn{\Omega} fixed.}
  \item{tol.out}{Convergence tolerance for outer loop of the exact MRCE algorithm.}
  \item{tol.in}{Convergence tolerance for the algorithm that minimizes the objective function, defined above,
  with \eqn{\Omega} fixed.}
  \item{kfold}{The number of folds to use when \code{method="cv"}.}
  \item{silent}{Logical: when \code{silent=FALSE} this function displays progress updates to the screen.}
  \item{eps}{The algorithm will terminate if the minimum diagonal entry of the current iterate's residual
  sample covariance is less than \code{eps}.  This may need adjustment depending on the scales of the variables.}
  \item{standardize}{Logical: should the columns of \code{X} be standardized so each has unit length and zero average.  The parameter estimates are returned on the original unstandarized scale.
  The default is \code{FALSE}.}
  \item{permute}{Logical: when \code{method="cv"}, should the subject indices be permutted?  The default is \code{FALSE}.}
}
\details{
Please see Rothman, Levina, and Zhu (2010)
for more information on the algorithm and model.
This version of the software uses the glasso algorithm (Friedman et al., 2008) through the R package \code{glasso}.
If the algorithm is running slowly, track its progress with \code{silent=FALSE}.  
In some cases, choosing \code{cov.tol=0.1} and \code{tol.out=1e-10} allows the algorithm to make
faster progress.  If one uses a matrix for \code{lam2}, consider setting \code{tol.in=1e-12}.

When \eqn{p \geq n},
the diagonal of the optimization variable corresponding 
to the inverse covariance matrix of the error is penalized.
Without diagonal penalization, if there exists a \eqn{\bar B} such
that the \eqn{q}th column of \eqn{Y} is equal to the \eqn{q}th
column of \eqn{X\bar B},
then a global minimizer of the objective function 
(defined above) does not exist.

The algorithm that minimizes the objective function, defined above,
with \eqn{Q} fixed uses a similar update strategy and termination 
criterion to those used by Friedman et al. (2010) in the corresponding R package \code{glmnet}.
}
\value{
A list containing
\item{Bhat}{This is \eqn{\hat\beta \in R^{p\times q}} defined above.  If \code{method="cv"}, 
then \code{best.lam1} and \code{best.lam2} defined below are used for \eqn{\lambda_1} and \eqn{\lambda_2}.}
\item{muhat}{This is the intercept estimate \eqn{\hat\mu \in R^q} defined above.
If \code{method="cv"}, 
then \code{best.lam1} and \code{best.lam2} defined below are used for \eqn{\lambda_1} and \eqn{\lambda_2}.}
\item{omega}{This is \eqn{\hat\Omega \in S_{+}^q} defined above.  If \code{method="cv"}, 
then \code{best.lam1} and \code{best.lam2} defined below are used for \eqn{\lambda_1} and \eqn{\lambda_2}.}
\item{mx}{This is \eqn{\bar x \in R^p} defined above.}
\item{my}{This is \eqn{\bar y \in R^q} defined above.}
\item{best.lam1}{The selected value for \eqn{\lambda_1} by cross validation.   Will be \code{NULL} unless \code{method="cv"}.}
\item{best.lam2}{The selected value for \eqn{\lambda_2} by cross validation.   Will be \code{NULL} unless \code{method="cv"}.}
\item{cv.err}{Cross validation error matrix with \code{length(lam1.vec)} rows and
 \code{length(lam2.vec)} columns.  Will be \code{NULL} unless \code{method="cv"}.}
}
\references{
Rothman, A. J., Levina, E., and Zhu, J. (2010) Sparse multivariate regression with covariance estimation.
  Journal of Computational and Graphical Statistics. 19: 947--962.
 
Jerome Friedman, Trevor Hastie, Robert Tibshirani (2008). 
Sparse inverse covariance estimation with the graphical lasso. 
Biostatistics, 9(3), 432-441.

Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
Regularization Paths for Generalized Linear Models via Coordinate Descent. 
Journal of Statistical Software, 33(1), 1-22.
}
\author{
Adam J. Rothman
}
\note{
The algorithm is fastest when \eqn{\lambda_1} and \eqn{\lambda_2} are large.  
Use \code{silent=FALSE} to 
check if the algorithm is converging before the total iterations exceeds \code{maxit.out}.
}
\examples{
set.seed(48105)
n=50
p=10
q=5

Omega.inv=diag(q)
for(i in 1:q) for(j in 1:q)
  Omega.inv[i,j]=0.7^abs(i-j)
out=eigen(Omega.inv, symmetric=TRUE)
Omega.inv.sqrt=tcrossprod(out$vec*rep(out$val^(0.5), each=q),out$vec)
Omega=tcrossprod(out$vec*rep(out$val^(-1), each=q),out$vec)

X=matrix(rnorm(n*p), nrow=n, ncol=p)
E=matrix(rnorm(n*q), nrow=n, ncol=q)\%*\%Omega.inv.sqrt
Beta=matrix(rbinom(p*q, size=1, prob=0.1)*runif(p*q, min=1, max=2), nrow=p, ncol=q)
mu=1:q

Y=rep(1,n)\%*\%t(mu) + X\%*\%Beta + E

lam1.vec=rev(10^seq(from=-2, to=0, by=0.5))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.5))
cvfit=mrce(Y=Y, X=X, lam1.vec=lam1.vec, lam2.vec=lam2.vec, method="cv")
cvfit

fit=mrce(Y=Y, X=X, lam1=10^(-1.5), lam2=10^(-0.5), method="single")
fit

lam2.mat=1000*(fit$Bhat==0)
refit=mrce(Y=Y, X=X, lam2=lam2.mat, method="fixed.omega", omega=fit$omega, tol.in=1e-12) 
refit
}
 