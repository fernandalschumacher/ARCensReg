
### AUXILIARY FUNCTIONS FOR AR NORMAL MODEL ###

## Covariance matrix (MatArp * sig2)
MatArp = function(phi,n) {
  p = length(phi)
  if (n==1) Rn = 1
  else Rn = toeplitz(ARMAacf(ar=phi, ma=0, lag.max = n-1))
  rhos = ARMAacf(ar=phi, ma=0, lag.max = p)[(1:p)+1]
  return(Rn/(1-sum(rhos*phi)))
}


## Transformation function: pi to phi
estphit = function(pit) {
  p = length(pit)
  Phi = matrix(0, ncol=p, nrow=p)
  if (p>1) {
    diag(Phi) = pit
    for (j in 2:p) {
      for (k in 1:(j-1)) {
        Phi[j,k] = Phi[j-1,k] - pit[j]*Phi[j-1,j-k]
      }
    }
    return(Phi[p,])
  }
  else return(pit)
}


## Transformation function: phi to pi
tphitopi = function(phit) {
  p = length(phit)
  Phi = matrix(0, ncol=p, nrow=p)
  Phi[p,] = phit
  if (p>1) {
    for (k in p:2) {
      for (i in 1:(k-1)) {
        Phi[k-1,i] = (Phi[k,i] + Phi[k,k]*Phi[k,k-i])/(1-Phi[k,k]^2)
      }
    }
    return(diag(Phi))
  }
  else return(phit)
}


## Estimate pi - case with censoring
lc = function(pi,D,n) {
  phi = estphit(pi)
  p = length(phi)
  lambda = matrix(c(-1,phi))
  spi = t(lambda)%*%D%*%lambda
  gp = 1
  for (i in 1:p) gp = gp*((1-pi[i]^2)^(-i))
  sig2hat = spi/n
  l = as.numeric(-n/2*log(sig2hat) -1/2*log(gp))
  return(-l)
}


## Estimate pi - case without censoring (Normal)
lcc = function(pi,y,x) {
  n = length(y)
  phi = estphit(pi)
  p = length(phi)
  lambda = matrix(c(-1,phi))
  betahat = solve(t(x)%*%solve(MatArp(phi,n))%*%x)%*%t(x)%*%solve(MatArp(phi,n))%*%y
  spi = t(lambda)%*%Dbeta(betahat,y,x,p)%*%lambda
  gp = 1
  for (i in 1:p) gp = gp*((1-pi[i]^2)^(-i))
  l = as.numeric(-n/2*log(spi) -1/2*log(gp))
  return(-l)
}


## Log-likelihood Normal model
LogVerosCens = function(cc, y, media, Psi, LI, LS){
  m = length(cc)
  gammai = media

  if(sum(cc)==0){
    ver = log(dmvnorm(as.vector(y), as.vector(gammai), Psi))
  }
  if(sum(cc)>0){
    if(sum(cc)==m){
      ver = log(pmvnorm(lower=c(LI),upper=c(LS),mean=c(media),sigma=Psi) + .Machine$double.xmin)
      
      } else {
      vero = numeric(2)
      vero[1] = dmvnorm(y[cc==0],gammai[cc==0,],Psi[cc==0,cc==0])
      inverse00 = solve(Psi[cc==0,cc==0])
      muc = gammai[cc==1,]+ Psi[cc==1,cc==0]%*%inverse00%*%(y[cc==0]-gammai[cc==0,])
      Sc = Psi[cc==1,cc==1]-Psi[cc==1,cc==0]%*%inverse00%*%Psi[cc==0,cc==1]
      vero[2] = pmvnorm(lower=c(LI[cc==1]),upper=c(LS[cc==1]),mean=c(muc),sigma=Sc)
      if(length(which(vero == 0)) > 0) vero[which(vero == 0)] = .Machine$double.xmin
      ver = sum(log(vero))
    }
  }
  obj.out = list(ver = ver)
  return(obj.out)
}


## Gibbs sampler (normal model)
amostradordegibbs = function(M, M0, nj, t1, cc1, y1, media, Gama, LI, LS){

    draws = matrix(NA, nrow=M, ncol=nj)
    draws[1,1:nj] = t1
    gammai = media

    if(sum(cc1)==0){
      for(i in 2:M){
        t1 = y1
        draws[i,1:nj] = t1
      }
    } else {
      
      if(sum(cc1)>0 & sum(cc1)==nj){
        g = as.vector(gammai)
        for(i in 2:M){
          t1 = as.vector(rtmvnorm(1, mean=g, sigma=(Gama), lower=LI, upper=LS, algorithm="gibbs", thinning=2))
          draws[i,1:nj] = t1
        }
      }
      if(sum(cc1)>0 & sum(cc1)<nj){
        if(sum(cc1)==1){
          g = gammai
          t1[cc1==0] = y1[cc1==0]
          inverse0 = solve(Gama[cc1==0,cc1==0])
          muc = g[cc1==1]+Gama[cc1==1,cc1==0]%*%inverse0%*%(y1[cc1==0]-g[cc1==0])
          muc = as.vector(muc)
          Sc = Gama[cc1==1,cc1==1]-Gama[cc1==1,cc1==0]%*%inverse0%*%Gama[cc1==0,cc1==1]
          Sc = as.numeric(Sc)
          for(i in 2:M){
            y_r = rtnorm(1, mean=muc, sd=(sqrt(Sc)), lower=LI[cc1==1], upper=LS[cc1==1])
            t1[cc1==1] = y_r
            draws[i,1:nj] = t1
          }
        } else{
          g = gammai
          t1[cc1==0] = y1[cc1==0]
          inverse0 = solve(Gama[cc1==0,cc1==0])
          muc = g[cc1==1]+Gama[cc1==1,cc1==0]%*%inverse0%*%(y1[cc1==0]-g[cc1==0])
          muc = as.vector(muc)
          Sc = Gama[cc1==1,cc1==1]-Gama[cc1==1,cc1==0]%*%inverse0%*%Gama[cc1==0,cc1==1]
          for(i in 2:M){
            y_r = rtmvnorm(1, mean=muc, sigma=(Sc), lower=LI[cc1==1], upper=LS[cc1==1], algorithm="gibbs", thinning=2)
            t1[cc1==1] = y_r
            draws[i,1:nj] = t1
          }
        }
      }
    }
  # Sample with burn-in (M0)
  amostragibbs = draws[(M0+1):M,]
  obj.out = list(amostragibbs = amostragibbs)
  return(obj.out)
}


## Derivatives
#################################
aphi = function(phi) ifelse(length(phi)==1,log(MatArp(phi,length(phi))),
                            log(det(MatArp(phi,length(phi)))))

Dbeta = function(beta,y,x,p) {
  n = length(y)
  D = matrix(0,p+1,p+1)
  for (ii in 1:(p+1)) {
    for (jj in 1:(p+1)) {
      D[ii,jj] = sum((y-x%*%beta)[ii:(n+1-jj)]*(y-x%*%beta)[jj:(n+1-ii)])
    }
  }
  return(D)
}
Dphi1 = function(beta,y,xx,p) matrix(Dbeta(beta,y,xx,p)[2:(p+1),1])
Dphiphi2 = function(beta,phi,y,xx,p) (Dbeta(beta,y,xx,p)[2:(p+1),2:(p+1)])%*%phi
einvM = function(phi,e) {
  n = length(e)
  invM = solve(MatArp(phi,n))
  return(t(e)%*%invM)
}

M1 = function(phi,yy){
  n= nrow(yy)
  return(sum(diag(yy%*%solve(MatArp(phi,n)))))
}

M1i = function(phi,yy,Di){
  n= nrow(yy)
  return(sum(diag(yy%*%solve(MatArp(phi,n))%*%Di)))
}

M2 = function(phi,vec1,vec2){
  n= nrow(vec2)
  return(vec1%*%solve(MatArp(phi,n))%*%vec2)
}

M3 = function(phi,vec1){
  n= length(vec1)
  return(t(vec1)%*%solve(MatArp(phi,n)))
}

Jt = function(theta,y,x) {
  l = ncol(x)
  n = length(y)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):length(theta)]
  p = length(phi)
  Mn = MatArp(phi,n)
  lambda = matrix(c(-1,phi))
  spi = t(lambda)%*%Dbeta(beta,y,x,p)%*%lambda
  invMn = solve(Mn)
  dbeta = 1/sig2*(t(x)%*%invMn%*%y - t(x)%*%invMn%*%x%*%beta)
  dsig2 = -n/2/sig2 +1/2/sig2^2*spi
  da = matrix(jacobian(aphi,phi))
  dphi = -1/sig2*(-Dphi1(beta,y,x,p) + Dphiphi2(beta,phi,y,x,p))-1/2*da
  return(rbind(dbeta,dsig2,dphi))
}

Ht = function(theta,y,x) {
  l = ncol(x)
  n = length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  lambda = matrix(c(-1,phi))
  spi = t(lambda)%*%Dbeta(beta,y,x,p)%*%lambda
  invMn = solve(Mn)
  dbetabeta = -1/sig2*(t(x)%*%invMn%*%x)
  dsig2sig2 = n/2/sig2^2 - 1*spi/sig2^3
  daa = (hessian(aphi,phi))
  dphiphi = -1/sig2*Dbeta(beta,y,x,p)[2:(p+1),2:(p+1)] - 1/2*daa
  dbetasig2 = -1/sig2^2*(t(x)%*%invMn%*%y - t(x)%*%invMn%*%x%*%beta )
  dD1beta = jacobian(Dphi1,beta,y=y,xx=x,p=p)
  dDphibeta = jacobian(Dphiphi2,beta,phi=phi,y=y,xx=x,p=p)
  dbetaphi = 1/sig2*(dD1beta - dDphibeta)
  dphisig = 1/sig2^2*(-Dphi1(beta,y,x,p)+ Dphiphi2(beta,phi,y,x,p))
  H = matrix(0,r,r)
  H[1:l,1:l] = dbetabeta
  H[l+1,l+1] = dsig2sig2
  H[(l+2):r,(l+2):r] = dphiphi
  H[l+1,1:l] = H[1:l,l+1] = dbetasig2
  H[(l+2):r,1:l] = dbetaphi
  H[1:l,(l+2):r] = t(dbetaphi)
  H[l+1,(l+2):r] = H[(l+2):r,l+1] = dphisig
  return(H)
}

Qt = function(theta,y,yy,x) {
  l = ncol(x)
  n = length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  invMn = solve(Mn)
  delta = sum(diag(yy%*%invMn)) - 2*t(y)%*%invMn%*%x%*%beta+t(beta)%*%t(x)%*%invMn%*%x%*%beta
  dbetabeta = -1/sig2*(t(x)%*%invMn%*%x)
  dsig2sig2 = n/2/sig2^2 - 1*delta/sig2^3
  daa = (hessian(aphi,phi)) #hessiana do log(det(Mn))
  dphiphi = - 1/2*daa -1/2/sig2*(hessian(M1,phi,yy=yy)+hessian(M2,phi,vec1=t(-2*y+x%*%beta),vec2=x%*%beta))
  dbetasig2 = -1/sig2^2*(t(x)%*%invMn%*%y - t(x)%*%invMn%*%x%*%beta )
  dbetaphi = 1/sig2*t(jacobian(M2,phi,vec1=t(y-x%*%beta),vec2=x))
  dphisig = 1/2/sig2^2*(jacobian(M1,phi,yy=yy)+ jacobian(M2,phi,vec1=t(-2*y+x%*%beta),vec2=x%*%beta))
  H = matrix(0,r,r)
  H[1:l,1:l] = dbetabeta
  H[l+1,l+1] = dsig2sig2
  H[(l+2):r,(l+2):r] = dphiphi
  H[l+1,1:l] = H[1:l,l+1] = dbetasig2
  H[(l+2):r,1:l] = dbetaphi
  H[1:l,(l+2):r] = t(dbetaphi)
  H[l+1,(l+2):r] = H[(l+2):r,l+1] = dphisig
  return(H)
}

######################################################
#local influence
######################################################
#ytil = y+w
deltaw = function(theta,y,x) {
  l = ncol(x)
  n = length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  invMn = solve(Mn)
  e = y-x%*%beta
  dbeta = 1/sig2*t(x)%*%invMn
  dsig2 = 1/sig2^2*t(e)%*%invMn
  dphi = -1/sig2*t(jacobian(einvM,phi,e=e))
  ddelta = rbind(dbeta,dsig2,dphi)
  return(ddelta)
}

################################################
#scheme 2 ----> Sigmatil = D(w)*Sigma
deltaSigi = function(theta,y,yy,x,i) {
  l = ncol(x)
  n = length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  invMn = solve(Mn)
  e = y-x%*%beta
  vec = rep(0,n);vec[i]=1
  Di = diag(vec)
  invMDi = invMn%*%Di
  dbeta = -1/sig2*t(x)%*%invMDi%*%e
  dsig2 = -1/2/sig2^2*(sum(diag(yy%*%invMDi))-2*t(y)%*%invMDi%*%x%*%beta+t(x%*%beta)%*%invMDi%*%x%*%beta)
  d1 = jacobian(M1i,phi,yy=yy,Di=Di)
  dphi = 1/2/sig2*t((d1)-2*jacobian(M2,phi,vec1 = t(y),vec2=Di%*%x%*%beta)+jacobian(M2,phi,vec1 = t(x%*%beta),vec2=Di%*%x%*%beta))
  ddelta = rbind(dbeta,dsig2,dphi)
  return(ddelta)
}

deltaSigi = function(theta,y,yy,x,i) {
  l = ncol(x)
  n = length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  invMn = solve(Mn)
  vec = rep(0,n);vec[i]=1
  Di = diag(vec)
  invMDi = invMn%*%Di
  DiinvM= Di%*%invMn
  dbeta = -1/2/sig2*t(x)%*%(2*DiinvM%*%y-(invMDi+DiinvM)%*%x%*%beta)
  dsig2 = -1/2/sig2^2*(sum(diag(yy%*%invMDi))-2*t(y)%*%invMDi%*%x%*%beta+t(x%*%beta)%*%invMDi%*%x%*%beta)
  d1 = jacobian(M1i,phi,yy=yy,Di=Di)
  dphi = 1/2/sig2*t((d1)-2*jacobian(M2,phi,vec1 = t(y),vec2=Di%*%x%*%beta)+jacobian(M2,phi,vec1 = t(x%*%beta),vec2=Di%*%x%*%beta))
  ddelta = rbind(dbeta,dsig2,dphi)
  return(ddelta)
}

#scheme 3 ----> x(w)=x+w*t(1)
deltaxp = function(theta,y,x,indp) {
  indp = matrix(indp,ncol=1)
  l = ncol(x)
  n = length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  invMn = solve(Mn)
  dbetaw = 1/sig2*(indp%*%t(y-x%*%beta)-as.numeric(t(indp)%*%beta)*t(x))%*%invMn
  dsigw = -as.numeric(t(indp)%*%beta)/sig2^2*(t(y-x%*%beta)%*%invMn)
  dphiw = as.numeric(t(indp)%*%beta)/sig2*t(jacobian(M3,phi,vec1=(y-x%*%beta)))
  ddelta = rbind(dbetaw,dsigw,dphiw)
  return(ddelta)
}
