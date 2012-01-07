compute.mrce <-
function(X,Y, lam1, lam2, tol.out=1e-3, tol.in=1e-4,
                         maxit.out=100, maxit.in=500,silent=FALSE)
{
  n=dim(X)[1]
  p=dim(X)[2]
   
  nlam=n*lam2
  xty=crossprod(X,Y)
  xtx=crossprod(X)

  old.B = qr.solve(xtx+nlam*diag(p), xty)
  k=0
  mab = sum(abs(old.B))
  while(1)
  {
    k=k+1
	  R = Y-X%*%old.B
	  samp.cov=crossprod(R)/n
	  g.out=glasso(s=samp.cov, rho=lam1, penalize.diagonal=FALSE, thr=1e-3)
    om=g.out$wi
	  om.i=g.out$w  
	  xtyom=xty%*%om 
     
	  warmstart=1
	  if(k == 1) warmstart=0
    B=rblasso(s=xtx, m=xtyom, om=om, nlam=nlam, tol=tol.in, sbols=mab, maxit=maxit.in, warm=warmstart, B0=old.B)		
    bdist = sum(abs(B-old.B))
    old.B=B
    if( (bdist < tol.out*mab) | (k > maxit.out))
      break	
  }
  if(silent ==FALSE) cat("Total outer iterations for MRCE : ", k, "\n")
  return(list(Bhat=old.B, omega=om))
}

