mrce <-
function(X,Y, lam1, lam2, lam0=NULL, method=c("exact", "approx"), B.start=NULL,
                          maxit.out=100, maxit.in=500, tol.out=1e-3, tol.in=1e-4, kfold=10, silent=TRUE)
{
  method=match.arg(method)
  doCV=(length(lam1)+length(lam2) > 2)
	doCVexact = doCV & (method=="exact")
  doCVapprox= doCV & (method=="approx")
	
	noCV=(length(lam1)+length(lam2) == 2)
	noCVexact = noCV & (method=="exact")
  noCVapprox= noCV & (method=="approx")
	
	Bhat=NULL
	omega=NULL
	lambda1=NULL
	lambda2=NULL
	lambda0=NULL
	cv.err=NULL
	cv.err.0=NULL
	
	if(noCVexact)
	{
	  tmp.out=compute.mrce(X=X,Y=Y, lam1=lam1, lam2=lam2, tol.out=tol.out, tol.in=tol.in,
                         maxit.out=maxit.out, maxit.in=maxit.in,silent=silent)
		 Bhat=tmp.out$Bhat
     omega=tmp.out$omega		 
	} else
	if(noCVapprox)
	{
	  tmp.out=compute.mrce.approx(X=X,Y=Y, lam1=lam1, lam2=lam2, B.start=B.start, lam.vec.0=lam0, 
                            maxit=maxit.in, tol=tol.in, kfold=kfold, silent=silent)
		Bhat=tmp.out$Bhat
    omega=tmp.out$omega
		cv.err.0=tmp.out$cv.err
    lambda0=tmp.out$lambda		
	} else
  if(doCVexact)
  {	
	  tmp.out=mrce.cv(X=X,Y=Y, lam.vec.1=lam1, lam.vec.2=lam2, maxit.out=maxit.out, maxit.in=maxit.in, 
		        tol.out=tol.out, tol.in=tol.in, kfold=kfold, silent=silent)
		Bhat=tmp.out$Bhat
    omega=tmp.out$omega
		lambda1=tmp.out$lambda1		
		lambda2=tmp.out$lambda2
    cv.err=tmp.out$cv.err			
  } else
	if(doCVapprox)
	{
	  tmp.out=mrce.approx.cv(X=X,Y=Y, lam.vec.0=lam0, lam.vec.1=lam1, lam.vec.2=lam2, 
                     maxit.in=maxit.in,  tol.in=tol.in, kfold=kfold, silent=silent)	
    Bhat=tmp.out$Bhat
    omega=tmp.out$omega
		lambda0=tmp.out$lambda0
		lambda1=tmp.out$lambda1		
		lambda2=tmp.out$lambda2
    cv.err=tmp.out$cv.err			
		cv.err.0=tmp.out$cv.err.0								 
	}
	
  return(list(Bhat=Bhat, omega=omega, lambda1=lambda1, lambda2=lambda2, lambda0=lambda0,
	   cv.err=cv.err, cv.err.0=cv.err.0))
}

