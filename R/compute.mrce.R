compute.mrce=function(X,Y, lam1, lam2, tol.out, tol.in, maxit.out, maxit.in, silent,
                      cov.tol, cov.maxit,informed=NULL, eps, standardize)
{
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  ## no diagonal penalty on Omega if n > p. 
  pdiag=(n <= p)
  if(is.null(informed))
  {
    if(!is.matrix(lam2))
      nlam=matrix(n*lam2, nrow=p, ncol=q) else nlam=n*lam2
    mx=apply(X, 2, mean)
    my=apply(Y, 2, mean)
    X=scale(X, center=mx, scale=FALSE)
    Y=scale(Y, center=my, scale=FALSE)
    if(standardize)
    {
      ## X is already centered so only need to scale
      sx=sqrt(apply(X^2, 2, sum)/n)
      X=scale(X, center=FALSE, scale=sx)
    }else sx=rep(1,p)
    yty=crossprod(Y)
    xty=crossprod(X,Y)
    xtx=crossprod(X)
    old.B=matrix(0, nrow=p, ncol=q)
    tolmult=sum(diag(yty)/n)
    tout=tol.out*tolmult
    residual.cov = yty/n
    sigma=diag(diag(residual.cov))
    om=diag(1/diag(residual.cov))
    omoff=om
    diag(omoff)=0
    old.obj=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
  } else
  {
    nlam=matrix(n*lam2, nrow=p, ncol=q)
    mx=informed$mx
    my=informed$my
    xtx=informed$xtx
    xty=informed$xty
    yty=informed$yty
    old.B=informed$Bhat
    sx=informed$sx
    if(standardize)
    {
      ## use the final iterate on the standardized
      ## scale as the initial iterate here
      old.B=old.B*informed$sx
      ## if informed and standardize then
      ## X is already standardized
    }
    om=informed$omega
    sigma=informed$sigma
    tolmult=sum(diag(yty)/n)
    tout=tol.out*tolmult
    residual.cov = crossprod(Y-X%*%old.B)/n
    omoff=om
    diag(omoff)=0
    old.obj=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
  }
  k=0  
  iterating=TRUE
  while(iterating)
  {
    k=k+1  
    if(min(diag(residual.cov)) < eps)
    {
      cat("A perfect fit occured. Terminated early.\n")
      break
    }   
    cov.out=glasso(s=residual.cov, rho=lam1, thr=cov.tol, maxit=cov.maxit, penalize.diagonal=pdiag)
    om=cov.out$wi   
    tolinmult=sum(yty*om)/n
    omoff=om
    diag(omoff)=0
    obj.after.omega=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
    if(!silent) cat("k =", k, "obj. fn. val. after Omega update is", obj.after.omega, "\n")  
    xtyom=xty%*%om  
    soft=xtyom - xtx%*%old.B%*%om + old.B*tcrossprod(diag(xtx), diag(om)) 
    outlasso=rblasso(s=xtx, m=xtyom, om=om, nlam=nlam, n=n,B0=old.B, soft=soft, objective=obj.after.omega, tol=(tol.in*tolinmult), maxit=maxit.in, quiet=silent)		
    old.B=outlasso$B
    residual.cov = crossprod(Y-X%*%old.B)/n
    new.obj=outlasso$obj
    bdist = old.obj-new.obj
    iterating= (bdist > tout) & (k <= maxit.out)
    old.obj=new.obj
    if(!silent) cat("k =", k, "obj. fn. val. after B update is", new.obj, "\n")
  }
  if(!silent) cat("Total outer iterations for MRCE : ", k, "\n")
  if(standardize)  ## compute old.B[,j] = old.B[,j]/sx using R's / operator:
    old.B=old.B/sx
  muhat=as.numeric(my - crossprod(old.B, mx))
  return(list(Bhat=old.B, muhat=muhat, omega=om, sigma=cov.out$w, mx=mx, my=my, sx=sx))
}
