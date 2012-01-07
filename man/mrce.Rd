\name{mrce}
\alias{mrce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute multivariate regression with covariance estimation (MRCE) estimators}
\description{Computes MRCE estimators and optionally does cross-validation}
\usage{
mrce(X, Y, lam1, lam2, lam0 = NULL, method = c("exact", "approx"), 
B.start = NULL, maxit.out = 100, maxit.in = 500, 
tol.out = 0.001, tol.in = 1e-04, kfold = 10, silent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{An \eqn{n} by \eqn{p} matrix of the prediction variables, where \eqn{n}
	is the number of observations and \eqn{p} is the number of predictors.  Should be uncentered
	if this function is used to cross-validate, and should be column-centered otherwise.}
  \item{Y}{An \eqn{n} by \eqn{q} matrix of the responses, where \eqn{q} is the number
          of responses.  Should be uncentered if this function is used to cross-validate, and should be column-centered otherwise.}
  \item{lam1}{Either a scalar or a vector of tuning parameter values for the penalty on the inverse error covariance.  
	If \code{lam1} is a vector with length at least 2, \code{kfold}-Fold CV will be performed}
  \item{lam2}{Either a scalar or a vector of tuning parameter values for the penalty on the regression coefficient matrix.  
	If \code{lam2} is a vector with length at least 2, \code{kfold}-Fold CV will be performed}
  \item{lam0}{Used for the approximate MRCE method (\code{method="approx"}).  
	Either a scalar or a vector of tuning parameter values for the penalty on the regression coefficient matrix
	when performing separate lasso regressions for each of the \eqn{q} responses all with the same tuning parameter.
  Used in the first step of the approximate MRCE algorithm.  
	If \code{lam0} is a vector with length at least 2, \code{kfold}-Fold CV will be performed for the first step
	of the approximate MRCE algorithm}
  \item{method}{\code{method="exact"} computes the exact MRCE solution, \code{method="approx"} computes the
	approximate MRCE solution.  The default is \code{method="exact"}}
  \item{B.start}{If \code{method="approx"}, the user can optionally specify a regression coefficient matrix to replace
	the first step of the approximate MRCE algorithm.  In this case, \code{lam0} is not used}
  \item{maxit.out}{The maximum allowable outer iterations of the exact MRCE algorithm.}
  \item{maxit.in}{The maximum allowable iterations of the algorithm that minimizes with respect to the regression coefficient matrix}
  \item{tol.out}{Convergence tolerance for outer-loop of the exact MRCE algorithm}
  \item{tol.in}{Convergence tolerance for the algorithm that minimizes with respect to the regression coefficient matrix}
  \item{kfold}{If cross validation is to be performed, the number of folds to use.}
  \item{silent}{Either \code{TRUE} or \code{FALSE}.  If \code{silent=FALSE} this function will display progress updates to the screen}
}
\details{No intercepts are to be included.  If the function is used to cross-validate for the tuning parameters, \code{X} and \code{Y}
should be uncentered.  For a description of the objective function see Rothman, Levina, and Zhu (2010).}
\value{
A list containing:
\item{Bhat}{Regression coefficient matrix estimate.  If the function performed cross-validation, then \code{Bhat} is
computed with the CV selected tuning parameters using column-centered versions of \code{X} and \code{Y}.}
\item{omega}{Inverse error covariance matrix estimate}
\item{lambda1}{The optimal tuing parameter for the penalty on the inverse error covariance, 
selected with cross validation.  If cross-validation is not performed, \code{lambda1} will be \code{NULL}}
\item{lambda2}{The optimal tuing parameter for the penalty on the regression coefficient matrix, 
selected with cross validation. If cross-validation is not performed, \code{lambda2} will be \code{NULL}}
\item{lambda0}{The optimal tuing parameter for the penalty on the regression coefficient matrix for the
first step in the approximation algorithm, selected with cross validation. If the approximation
of cross-validation are not performed, \code{lambda0} will be \code{NULL}}
\item{cv.err}{Cross validation error matrix.  Will be \code{NULL} if cross validation over \code{lam1} and \code{lam2} is 
not performed}
\item{cv.err.0}{Cross validation error vector for the first step of the approximation.
Will be \code{NULL} if cross validation over \code{lam0} is not performed}
}
\references{
Rothman, A.J., Levina, E., and Zhu, J. (2010) Sparse multivariate regression with covariance estimation.
  Journal of Computational and Graphical Statistics. 19: 947--962.
}
\author{
Adam J. Rothman
}
\note{
Certain values of the tuning parameters could make \code{method="exact"} run very slowly.  Use \code{silent=FALSE} to 
check if the algorithm is converging before the total iterations exceeds \code{maxit.out}.
}
\examples{
set.seed(5)
n=50
p=20
q=20

true.B=matrix( rnorm(p*q), nrow=p, ncol=q)
indi = matrix( rbinom(p*q, 1, 0.1), nrow=p, ncol=q)
true.B = true.B*indi
	
X.tr=matrix(rnorm(n*p), nrow=n, ncol=p)
E.tr=matrix(rnorm(n*q), nrow=n, ncol=q)
Y.tr=X.tr\%*\%true.B+E.tr

mx = apply(X.tr, 2, mean)
my = apply(Y.tr, 2, mean)
X.tr=scale(X.tr, center=mx, scale=FALSE)
Y.tr=scale(Y.tr, center=my, scale=FALSE)
	
mrce.out=mrce(X=X.tr,Y=Y.tr, lam1=0.1, lam2=0.1, method="exact")
mrce.out$Bhat


}
 
