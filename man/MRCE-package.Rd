\name{MRCE-package}
\alias{MRCE-package}
\alias{MRCE}
\docType{package}
\title{Multivariate regression with covariance estimation}
\description{
Computes the MRCE covariance estimators (Rothman, Levina, and Zhu, 2010) 
and has the dataset \code{stock04} used in Rothman, Levina, and Zhu (2010),
originally analyzed in Yuan et. al (2007)
}
\details{The primary function is \code{mrce} (type \code{help(mrce)} for more information). The dataset is called \code{stock04}
(type \code{help(stock04)} for additional information).}
\author{Adam J. Rothman

Maintainer: Adam J. Rothman <arothman@umn.edu>
}
\references{
Rothman, A.J., Levina, E., and Zhu, J. (2010) Sparse multivariate regression with covariance estimation.
  Journal of Computational and Graphical Statistics.  19: 974--962.
}
 
\keyword{ package }
 
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
