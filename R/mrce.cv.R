mrce.cv <-
function(X,Y, lam.vec.1, lam.vec.2,  
                          maxit.out=100, maxit.in=500, tol.out=1e-4, tol.in=1e-4, kfold=5, silent=FALSE)
{
  n=dim(Y)[1]
  err.mat = matrix(0, nrow=length(lam.vec.1), ncol=length(lam.vec.2))
  ind=sample(1:n)
	
	for (k in 1:kfold)
  {
	  if(!silent) cat("Staring fold : ", k, "\n")
    foldind = ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
	  X.tr=X[-foldind, ]
	  Y.tr=Y[-foldind, ]
	  X.va=X[foldind, ]
		Y.va=Y[foldind, ]
 
	  mtrx=apply(X.tr, 2, mean)
	  X.tr=scale(X.tr, scale=FALSE, center=mtrx)
	  X.va=scale(X.va, scale=FALSE, center=mtrx)
	
	  mtry=apply(Y.tr, 2, mean)
	  Y.tr=scale(Y.tr, scale=FALSE, center=mtry)
	  Y.va=scale(Y.va, scale=FALSE, center=mtry)
	  for(i in 1:length(lam.vec.1))
	  {
	    for(j in 1:length(lam.vec.2))
    	{
	      lam1=lam.vec.1[i]
	      lam2=lam.vec.2[j] 
	      tmp.out=compute.mrce(X=X.tr,Y=Y.tr, lam1=lam1, lam2=lam2, 
		                        tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, maxit.in=maxit.in,  silent=silent)
	      err.mat[i,j]=err.mat[i,j]+mean((Y.va-X.va%*%tmp.out$Bhat)^2)				
	    }
    }
	}
	
  ## find the (i,j) for the minimum of err.mat
	tmp = which.min(err.mat) %% (dim(err.mat)[1])
  tmp = (tmp != 0)*tmp + (tmp == 0)*(dim(err.mat)[1])
  best.i=tmp
  best.j=which.min(err.mat[tmp,])
 
  X=scale(X, scale=FALSE)
	Y=scale(Y, scale=FALSE)
	 
  tmp.out=compute.mrce(X=X,Y=Y, lam1=lam.vec.1[best.i], lam2=lam.vec.2[best.j],  
	                        tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, maxit.in=maxit.in, silent=silent) 
	 
	best.B = tmp.out$Bhat
  best.omega = tmp.out$omega
	
  return(list(Bhat=best.B, omega=best.omega, lambda1=lam.vec.1[best.i], lambda2=lam.vec.2[best.j],
	   cv.err=err.mat))
}

