
##############################################################
# Parameter estimation - Normal innovations #
##############################################################
#cc must be 1 for miss

SAEM = function(cc, LI, LS, y, x, p, M, perc, MaxIter, pc, miss, tol, show_ep, quiet){
  
  if (!quiet) pb = txtProgressBar(min = 0, max = MaxIter, style = 3)
  
  # Initial values
  m = length(y)
  if (length(miss)==0) {
    beta1 = solve(t(x)%*%x)%*%t(x)%*%y
    l = length(beta1)
    pi1 = as.numeric(pacf((y - x%*%beta1), lag.max = p, plot=F)$acf)
    pi1 = optim(pi1, lcc, y=y, x=x, lower=rep(-.999,p), upper=rep(0.999,p), method='L-BFGS-B')$par
    phi1 = estphit(pi1)
    lambda = matrix(c(-1,phi1))
    beta1 = solve(t(x)%*%solve(MatArp(phi1,m))%*%x)%*%t(x)%*%solve(MatArp(phi1,m))%*%y
    spi = t(lambda)%*%Dbeta(beta1,y,x,p)%*%lambda
    sigmae = as.numeric(spi/m)
  } else {
    beta1 = solve(t(x[-miss,])%*%x[-miss,])%*%t(x[-miss,])%*%y[-miss]
    l = length(beta1)
    pi1 = as.numeric(pacf((y - x%*%beta1)[-miss],lag.max=p,plot=F)$acf)
    pi1 = optim(pi1,lcc,y=y[-miss],x=as.matrix(x[-miss,]),lower = rep(-.999,p), upper = rep(0.999,p),method='L-BFGS-B')$par
    phi1 = estphit(pi1)
    lambda = matrix(c(-1,phi1))
    beta1 = solve(t(as.matrix(x[-miss,]))%*%solve(MatArp(phi1,m-length(miss)))%*%as.matrix(x[-miss,]))%*%
      t(as.matrix(x[-miss,]))%*%solve(MatArp(phi1,m-length(miss)))%*%as.matrix(y[-miss])
    spi = t(lambda)%*%Dbeta(beta1,y[-miss],as.matrix(x[-miss,]),p)%*%lambda
    sigmae = as.numeric(spi/(m-length(miss)))
  }
  teta = c(beta1, sigmae, phi1)
  r = length(teta)
  
  tempoi = Sys.time()
  if ((sum(cc)==0) & (length(miss)==0)) {
    teta1 = teta
    media = x%*%beta1
    V = MatArp(phi1, m)
    Psi = sigmae*V
    #
    H = -Ht(teta1, y, x)
    ep_par = sqrt(diag(solve(H)))
    logver = (LogVerosCens(cc, y, media, Psi, LI, LS)$ver)
    npar = length(c(teta1))
    loglik = logver
    AICc = -2*loglik +2*npar
    AICcorr = AICc + ((2*npar*(npar+1))/(m-npar-1))
    BICc = -2*loglik +log(m)*npar
    tempof = Sys.time()
    dift = tempof - tempoi
    
    if (!quiet) setTxtProgressBar(pb, MaxIter)
    if (show_ep) resultados = list(beta=c(beta1), sigma2=sigmae, phi=phi1, pi1=pi1, theta=teta1, SE=ep_par)
    else resultados = list(beta=c(beta1), sigma2=sigmae, phi=phi1, pi1=pi1, theta=teta1)
    resultados$loglik=loglik 
    resultados$AIC=AICc
    resultados$BIC=BICc
    resultados$AICcorr=AICcorr
    resultados$yest=y
    resultados$yyest=y%*%t(y)
    resultados$x=x
    obj.out = list(res=resultados, time=dift)

  } else {
    ## SAEM algorithm ################################
    MG = round(M/(1-perc),0)
    M0 = MG - M #burn in
    Theta = matrix(NA,MaxIter,length(teta))

    #criterio <- 10000
    critval = critval2 = 1
    delta1 = 0.001
    delta2 = tol

    count = 0
    media = x%*%beta1
    V = MatArp(phi1,m)
    Psi = sigmae*V

    # Sequence of decreasing positive numbers: smoothing parameter
    if (pc==1){
      seqq = rep(1,pc*MaxIter)
    } else {
      seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
      seqq = c(rep(1,MaxIter-length(seqq)),seqq)
    }

    SAEM_ss = array(data=0,dim=c(MaxIter+1))
    SAEM_xx = array(data=0,dim=c(MaxIter+1,l,l))
    SAEM_xy = array(data=0,dim=c(MaxIter+1,l))
    SAEM_y = array(data=0,dim=c(MaxIter+1,sum(cc)))
    SAEM_yy = array(data=0,dim=c(MaxIter+1,m,m))
    SAEM_D = array(data=0,dim=c(MaxIter+1,p+1,p+1))
    if (show_ep & (sum(cc)!=0)) {
      SAEM_delta = array(data=0,dim=c(MaxIter+1,r))
      SAEM_G = array(data=0,dim=c(MaxIter+1,r,r))
      SAEM_H = array(data=0,dim=c(MaxIter+1,r,r))
    }
    tyi = y
    
    while(critval2 < 3) {
      count = count + 1
      if (!quiet) setTxtProgressBar(pb, count)

      # E-1: Sampling step
      t1 = tyi
      gibbs = amostradordegibbs(MG,M0,m,t1,cc,y,media,Psi,LI,LS)
      amostragibbs = gibbs$amostragibbs
      uyi = matrix(amostragibbs[,1:m],nrow=M,ncol=m)

      # E-2: Stochastic approximation
      somaD  = matrix(0, p+1, p+1)
      somass = 0
      somaxx = matrix(0, l, l)
      somaxy = matrix(0, l, 1)
      somay  = matrix(0, sum(cc), 1)
      somayy = matrix(0, m, m)
      somadelta = matrix(0, r, 1)
      somaG = matrix(0, r, r)
      invV = solve(V)
      for (k in 1:M) {
        yi = matrix(uyi[k,], nrow=m, ncol=1)
        somass = somass + (t(yi-media)%*%invV%*%(yi-media))
        somaxx = somaxx + (t(x)%*%invV%*%x)
        somaxy = somaxy + t(x)%*%invV%*%(yi)
        somay  = somay + as.matrix(yi[cc==1,])
        somayy = somayy + yi%*%t(yi)
        if (show_ep & (sum(cc)!=0)) {
          J = Jt(teta,yi,x)
          H = Ht(teta,yi,x)
          somadelta = somadelta + J
          somaG = somaG + (-H - (J)%*%t(J))
        }
        somaD = somaD + Dbeta(beta1,yi,x,p)
      }
      E_D  = 1/M*somaD
      E_ss = (1/M)*somass
      E_xx = (1/M)*somaxx
      E_xy = (1/M)*somaxy
      E_y  = (1/M)*somay
      E_yy = (1/M)*somayy
      if (show_ep & (sum(cc)!=0)) {
        E_delta = 1/M * somadelta
        E_G = 1/M*somaG
      }
      
      # E-2: Stochastic approximation
      SAEM_D[count+1,,] = SAEM_D[count,,] + seqq[count]*(E_D - SAEM_D[count,,])
      SAEM_ss[count+1] = SAEM_ss[count] + seqq[count]*(E_ss - SAEM_ss[count])
      SAEM_xx[count+1,,] = SAEM_xx[count,,] + seqq[count]*(E_xx - SAEM_xx[count,,])
      SAEM_xy[count+1,] = SAEM_xy[count,] + seqq[count]*(E_xy - SAEM_xy[count,])
      SAEM_y[count+1,] = SAEM_y[count,] + seqq[count]*(E_y - SAEM_y[count,])
      SAEM_yy[count+1,,] = SAEM_yy[count,,] + seqq[count]*(E_yy - SAEM_yy[count,,])
      if (show_ep & (sum(cc)!=0)) {
        SAEM_delta[count+1,] = SAEM_delta[count,] + seqq[count]*(E_delta - SAEM_delta[count,])
        SAEM_G[count+1,,] = SAEM_G[count,,] + seqq[count]*(E_G - SAEM_G[count,,])
        SAEM_H[count+1,,] = SAEM_G[count+1,,] - (SAEM_delta[count+1,])%*%t(SAEM_delta[count+1,])
      }
      tyi = y;   tyi[cc==1] = SAEM_y[count+1,]
      
      # CM: Conditional maximization
      beta1 = solve(SAEM_xx[count+1,,])%*%SAEM_xy[count+1,]
      media = x%*%beta1

      sigmae = (1/m)*(SAEM_ss[count+1])
      sigmae = as.numeric(sigmae)

      pi1 = optim(pi1,lc,lower = rep(-.999,p), upper = rep(0.999,p), n=m, D = SAEM_D[count+1,,],method='L-BFGS-B')$par
      phi1 = estphit(pi1)
      teta1 = c(beta1, sigmae, phi1)

      V = MatArp(phi1,m)
      Psi = sigmae*V

      criterio2 = sqrt((teta1/teta-1)%*%(teta1/teta-1))
      if(max(criterio2) < delta2){critval2 = critval2+1}else{critval2 = 0}
      if(count == MaxIter){critval2 = 10}

      Theta[count,] = teta1
      teta = teta1
    } # End while
    
    if (!quiet) setTxtProgressBar(pb, MaxIter)
    Theta = Theta[1:count,]
    logver = LogVerosCens(cc, y, media, Psi, LI, LS)$ver
    tempof = Sys.time()
    dift = tempof - tempoi
    npar = length(c(teta1))
    loglik = logver
    AICc = -2*loglik +2*npar
    AICcorr = AICc + ((2*npar*(npar+1))/(m-npar-1))
    BICc = -2*loglik +log(m)*npar
    
    if (show_ep & (sum(cc)!=0)) {
      vartheta = solve(SAEM_H[count+1,,])
      ep_par = sqrt(diag(vartheta))
    }
    
    yest = numeric(m)
    yest[cc==0] = y[cc==0]
    yest[cc==1] = SAEM_y[count+1,]
    
    if (show_ep) resultados = list(beta=c(beta1), sigma2=sigmae, phi=phi1, pi1=pi1, theta=teta1, SE=ep_par)
    else resultados = list(beta=c(beta1), sigma2=sigmae, phi=phi1, pi1=pi1, theta=teta1)
    resultados$loglik = loglik 
    resultados$AIC = AICc
    resultados$BIC = BICc
    resultados$AICcorr = AICcorr
    resultados$yest = yest
    resultados$yyest = SAEM_yy[count+1,,]
    resultados$x = x
    resultados$iter = count
    resultados$criteria=criterio2
    obj.out = list(res=resultados, time=dift, Theta=Theta)
  } # End else
  
  return(obj.out)
}
