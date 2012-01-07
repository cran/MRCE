mrce.approx.cv <-
function(X,Y, lam.vec.0, lam.vec.1, lam.vec.2, 
                            maxit.in=500,  tol.in=1e-4, kfold=5, silent=FALSE)
{
  n=dim(Y)[1]
  err.mat = matrix(0, nrow=length(lam.vec.1), ncol=length(lam.vec.2))
  ind=sample(1:n)
  
	multi.lasso.out=multi.lasso.cv(X=X,Y=Y, lam.vec=lam.vec.0)
	B.step1=multi.lasso.out$Bhat
	lambda0=multi.lasso.out$lambda
	cv.err.0=multi.lasso.out$cv.err

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
	  
	  B.start=multi.lasso(X=X.tr,Y=Y.tr, lam=lambda0)
		for(i in 1:length(lam.vec.1))
	  {
	    for(j in 1:length(lam.vec.2))
    	{
	      lam1=lam.vec.1[i]
	      lam2=lam.vec.2[j]	
				tmp.out=compute.mrce.approx(X=X.tr,Y=Y.tr, lam1=lam1, lam2=lam2, B.start=B.start,   
                            maxit=maxit.in, tol=tol.in, silent=silent)
				
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
	
	tmp.out=compute.mrce.approx(X=X,Y=Y, lam1=lam.vec.1[best.i], lam2=lam.vec.2[best.j], 
		                       B.start=B.step1, maxit=maxit.in, tol=tol.in, silent=silent)
	
	Bhat = tmp.out$Bhat
  omega = tmp.out$omega
	
  return(list(Bhat=Bhat, omega=omega, lambda1=lam.vec.1[best.i], lambda2=lam.vec.2[best.j], lambda0=lambda0,
	   cv.err=err.mat, cv.err.0=cv.err.0))
}

