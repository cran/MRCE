compute.mrce.approx <-
function(X,Y, lam1, lam2, B.start=NULL, lam.vec.0=NULL, 
                            maxit=500, tol=1e-4, kfold=5, silent=TRUE)
{
  lambda=NULL
	cv.err=NULL
	if(is.null(B.start))
	{
	  if(!silent) cat("cross-validating for the initial B matrix\n")
	  multi.lasso.out=multi.lasso.cv(X=X,Y=Y, lam.vec=lam.vec.0, kfold=kfold)
	  old.B=multi.lasso.out$Bhat
		lambda=multi.lasso.out$lambda
		cv.err=multi.lasso.out$cv.err
		X=scale(X, scale=FALSE)
	  Y=scale(Y, scale=FALSE)
	}else
	{
	  old.B = B.start
	}

  n=dim(Y)[1]
	p=dim(X)[2]
	nlam=n*lam2
  xty=crossprod(X,Y)
  xtx=crossprod(X)
	
	
	R = Y-X%*%old.B
	samp.cov=crossprod(R)/n
  mab = sum(abs(qr.solve(xtx+nlam*diag(p), xty)))
  if(!silent) cat("Solving for omega\n")	
	g.out=glasso(s=samp.cov, rho=lam1, penalize.diagonal=FALSE, thr=1e-3)
	om=g.out$wi
  xtyom=crossprod(X,Y)%*%om 
	xtx=crossprod(X)
	if(!silent) cat("Solving for B\n")
  Bhat=rblasso(s=xtx, m=xtyom, om=om, nlam=(n*lam2), tol=tol, sbols=mab, maxit=maxit, warm=0, B0=old.B)	
	
  return(list(Bhat=Bhat, omega=om, cv.err=cv.err, lambda=lambda))
}

