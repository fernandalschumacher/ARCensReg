## The autoregressive AR(p) Student-t model ##

##############################################################
# Random numbers from the multivariate normal distribution #
##############################################################

rtmvnormal = function(n, mean, sigma, lower, upper, thinning){
  if (!is.positive.definite(sigma)){ stop("The variance matrix could be positive definite.",call.=FALSE) }
  x0 = ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  
  if (length(mean)==1){ ret = rtnorm(n, mean=mean, sd=sqrt(sigma), lower=lower, upper=upper)
  } else {   ret = rtmvnorm(n, mean=mean, sigma=sigma, lower=lower, upper=upper, algorithm="gibbs",
                             burn.in.samples=0, start.value=x0, thinning=thinning) }
  return (ret)
}


##############################################################
# Function to maximize nu #
##############################################################

Maxnu = function(e.nu, SAEMu, SAEMlu){
  ni = length(SAEMu) # n-p
  fmax = as.numeric()
  fmax = 0.5*ni*(e.nu*log(0.5*e.nu) - 2*log(gamma(0.5*e.nu))) + 0.5*e.nu*(SAEMlu - sum(SAEMu))
  return(fmax)  # Function to be maximized
}


##############################################################
# Gibbs sampler #
##############################################################

Gibbs_samplerCLRT = function(X, yobs, SAEMu, cc, LI, LS, phi, beta, sigma2, nu, M1, M0, fixnu, show_ep){
 nj = length(yobs) # n
 pj = length(phi)  # p
 diff1 = nj - pj   # n - p
 drawsY = vector("numeric", length=nj)      # Sum Y
 drawsU = vector("numeric", length=diff1)   # Sum u
 drawslU  = 0                               # Sum i=1^m (sum j=p+1^n log u_j)
 drawsUY2 = 0                              # Sum i=1^m (sum j=p+1^n u_jy_j^2)
 drawsUYi = vector("numeric", length=diff1) # Sum u*Y
 drawsUYZ = vector("numeric", length=pj)    # Sum i=1^m (sum j=p+1^n u_jy_ijz_j-1) 
 drawsUZi = matrix(0,nrow=pj, ncol=diff1)   # Sum u_i z_i-1
 drawsUZ2 = matrix(0,nrow=pj, ncol=pj)      # Sum i=1^m (sum j=p+1^n u_jz_i-1t(z_i-1))
 sampler.y = yobs
 if (show_ep){
   if (!fixnu){ deltaM = matrix(0, nrow=(pj+length(beta)+2), ncol=(pj+length(beta)+2))     # The second part of the information matrix
   } else {     deltaM = matrix(0, nrow=(pj+length(beta)+1), ncol=(pj+length(beta)+1)) }   # The second part of the information matrix 
 }
 
 if (sum(cc)==0){
   diff.media = sampler.y - X%*%beta
   
   for (k in 1:M1){
     sampler.u  = vector("numeric", length=diff1)
     samplerUYZ = vector("numeric", length=pj)    
     samplerUZi = matrix(0, nrow=pj, ncol=diff1)     
     samplerUZ2 = matrix(0, nrow=pj, ncol=pj)        
     
     # Sampling from u|Y,theta
     for (i in 1:diff1){
       bi = nu + (diff.media[pj+i] - sum(phi*diff.media[(pj+i-1):i]))^2/sigma2
       sampler.u[i] = rgamma(1, shape=0.5*(nu+1), rate=0.5*bi)
       
       samplerUYZ = samplerUYZ + sampler.u[i]*sampler.y[pj+i]*sampler.y[(pj+i-1):i]
       samplerUZi[,i] = sampler.u[i]*sampler.y[(pj+i-1):i]
       samplerUZ2 = samplerUZ2 + sampler.u[i]*(sampler.y[(pj+i-1):i]%*%t(sampler.y[(pj+i-1):i]))
      } # End for
     
     if (k > M0){
       drawsU   = drawsU + sampler.u
       drawslU  = drawslU + sum(log(sampler.u))
       drawsUY2 = drawsUY2 + sum(sampler.u*(sampler.y[-(1:pj)])^2)
       drawsUYi = drawsUYi + sampler.u*sampler.y[-(1:pj)]
       drawsUYZ = drawsUYZ + samplerUYZ
       drawsUZi = drawsUZi + samplerUZi
       drawsUZ2 = drawsUZ2 + samplerUZ2
       if (show_ep){
         score  = Gradient(X, sampler.u, sampler.y, phi, beta, sigma2, nu, fixnu)
         deltaM = deltaM + score%*%t(score)
       }
     } # End if
   } # End for
   EY = sampler.y
 } # Sampling without censored observations
 
 if (sum(cc)>=1 & sum(cc)<=nj){
   mean.vc = ComputeMean(yobs, X, beta, phi)
   mean1   = mean.vc$mvc
   cc1 = cc[-(1:pj)]; Li = LI[-c(1:pj)]; Ls = LS[-c(1:pj)]
   sampler.u = SAEMu
   
   for (k in 1:M1){
     samplerUYZ = vector("numeric", length=pj)    
     samplerUZi = matrix(0, nrow=pj, ncol=diff1)     
     samplerUZ2 = matrix(0, nrow=pj, ncol=pj) 
     
     # Sampling from Yc|u,Yo,theta
     var1 = ComputeVar(sigma2, mean.vc$Pphi, phi, sampler.u)
     y.np = sampler.y[-(1:pj)]
     inversa.vc = var1[cc1==1,cc1==0]%*%inversa(var1[cc1==0,cc1==0])
     mean.c = mean1[cc1==1] + inversa.vc%*%(y.np[cc1==0]-mean1[cc1==0])
     mean.c = as.vector(mean.c)   # Mean
     var.c  = var1[cc1==1,cc1==1] - inversa.vc%*%var1[cc1==0,cc1==1]
     var.c  = (var.c + t(var.c))/2 # Variance matrix
     random.y = rtmvnormal(1, mean=mean.c, sigma=var.c, lower=Li[cc1==1], upper=Ls[cc1==1], thinning=2)
     y.np[cc1==1] = random.y
     sampler.y[-c(1:pj)] = y.np
     
     # Sampling from u|Y,theta
     diff.media = sampler.y - X%*%beta
     for (i in 1:diff1){
       bi = nu + (diff.media[pj+i] - sum(phi*diff.media[(pj+i-1):i]))^2/sigma2
       sampler.u[i] = rgamma(1,shape=0.5*(nu+1),rate=0.5*bi)
       
       samplerUYZ = samplerUYZ + sampler.u[i]*sampler.y[pj+i]*sampler.y[(pj+i-1):i]
       samplerUZi[,i] = sampler.u[i]*sampler.y[(pj+i-1):i]
       samplerUZ2 = samplerUZ2 + sampler.u[i]*(sampler.y[(pj+i-1):i]%*%t(sampler.y[(pj+i-1):i]))
      } # End for
     
     if (k > M0){
       drawsY   = drawsY + sampler.y
       drawsU   = drawsU + sampler.u
       drawslU  = drawslU + sum(log(sampler.u))
       drawsUY2 = drawsUY2 + sum(sampler.u*(sampler.y[-(1:pj)])^2)
       drawsUYi = drawsUYi + sampler.u*sampler.y[-(1:pj)]
       drawsUYZ = drawsUYZ + samplerUYZ
       drawsUZi = drawsUZi + samplerUZi
       drawsUZ2 = drawsUZ2 + samplerUZ2
       if (show_ep){
         score  = Gradient(X, sampler.u, sampler.y, phi, beta, sigma2, nu, fixnu)
         deltaM = deltaM + score%*%t(score)
       }
     } # End if
   } # End for
   EY = drawsY/(M1 - M0)
 } # Sampling censored variables
 
 # Estimating the expectations
 m = M1 - M0
 EU  = drawsU/m
 ElU = drawslU/m
 EUY2 = drawsUY2/m
 EUYi = drawsUYi/m
 EUYZ = drawsUYZ/m
 EUZi = drawsUZi/m
 EUZ2 = drawsUZ2/m
 if (show_ep){
   EdeltaM = deltaM/m
   return (list(EY=EY, EU=EU, ElU=ElU, EUY2=EUY2, EUYi=EUYi, EUYZ=EUYZ, EUZi=EUZi, EUZ2=EUZ2, Edelta=EdeltaM))
 } else {
   return (list(EY=EY, EU=EU, ElU=ElU, EUY2=EUY2, EUYi=EUYi, EUYZ=EUYZ, EUZi=EUZi, EUZ2=EUZ2))
 }
}

##############################################################
# Parameter estimation #
##############################################################

SAEM_temporalT = function(cens, LI, LS, y, x, p, x_pred, tol, M, perc, MaxIter, pc, 
                           nufix, show_ep, quiet){
  if (!quiet){ pb = txtProgressBar(min = 0, max = MaxIter, style = 3) }
  yobs = y
  n = length(yobs)
  X = x
  p = p
  q = ncol(X)
  model0 = arima(yobs, order=c(p=p,d=0,q=0), xreg=X, method="ML", include.mean=FALSE)
  beta   = c((model0$coef)[-(1:p)])
  phi    = c((model0$coef)[1:p])
  sigma2 = model0$sigma2
  if (is.null(nufix)) { 
    nu = 10;  fixed.nu = FALSE 
    theta = c(beta, sigma2, phi, nu)
  } else {  
    nu = nufix;  fixed.nu = TRUE
    theta = c(beta, sigma2, phi)
  }
  Theta = theta
  
  # Stop criterion
  criterio = criterio2 = 10
  count = 0
  
  ## SAEM algorithm ################################
  MG = round(M/(1 - perc),0) # Number of samples to generate
  M0 = MG - M                # Number of burn samples
  # Sequence of decreasing positive numbers: smoothing parameter
  if (pc==1){
    seqq = rep(1,MaxIter)
  } else {
    seqq = c(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter)))
    seqq = c(rep(1,MaxIter-length(seqq)),seqq)
  }
  SAEM.Y = vector("numeric", length=n)
  SAEM.U = runif(n-p, 0, 1)
  SAEM.lU  = 0
  SAEM.UY2 = 0
  SAEM.UYi = vector("numeric", length=n-p)
  SAEM.UYZ = vector("numeric", length=p)
  SAEM.UZi = matrix(0, nrow=p, ncol=n-p)
  SAEM.UZ2 = matrix(0, nrow=p, ncol=p)
  if (show_ep){
    if (!fixed.nu){ SAEM.delta = matrix(0,nrow=(p+q+2),ncol=(p+q+2)) 
    } else { SAEM.delta = matrix(0,nrow=(p+q+1),ncol=(p+q+1)) }
  }
  
  initime = Sys.time()
  while (criterio > tol){
    count = count + 1
    if (!quiet){ setTxtProgressBar(pb, count) }
    
    # E-1: Sampling step
    amostras = Gibbs_samplerCLRT(X, yobs, SAEM.U, cens, LI, LS, phi, beta, sigma2, nu, MG, M0, fixed.nu, show_ep)
    
    # E-2: Stochastic approximation
    SAEM.Y = SAEM.Y + seqq[count]*(amostras$EY - SAEM.Y)
    SAEM.U = SAEM.U + seqq[count]*(amostras$EU - SAEM.U)
    SAEM.lU  = SAEM.lU + seqq[count]*(amostras$ElU - SAEM.lU)
    SAEM.UY2 = SAEM.UY2 + seqq[count]*(amostras$EUY2 - SAEM.UY2)
    SAEM.UYi = SAEM.UYi + seqq[count]*(amostras$EUYi - SAEM.UYi)
    SAEM.UYZ = SAEM.UYZ + seqq[count]*(amostras$EUYZ - SAEM.UYZ)
    SAEM.UZi = SAEM.UZi + seqq[count]*(amostras$EUZi - SAEM.UZi)
    SAEM.UZ2 = SAEM.UZ2 + seqq[count]*(amostras$EUZ2 - SAEM.UZ2)
    if (show_ep){ SAEM.delta = SAEM.delta + seqq[count]*(amostras$Edelta - SAEM.delta) }
    
    # CM: Conditional maximization
    media1 = X%*%beta
    part1 = vector("numeric", length=p)
    part2 = matrix(0, ncol=p, nrow=p)
    for (i in 1:(n-p)){
      part1 = part1 + SAEM.UYi[i]*media1[(p+i-1):i,] + media1[p+i,]*matrix(SAEM.UZi[,i]) - SAEM.U[i]*media1[p+i,]*media1[(p+i-1):i,]
      part2 = part2 + SAEM.UZi[,i]%*%t(media1[(p+i-1):i,]) + media1[(p+i-1):i,]%*%t(SAEM.UZi[,i]) - SAEM.U[i]*(media1[(p+i-1):i,]%*%t(media1[(p+i-1):i,]))
    }
    uy2 = SAEM.UY2 - 2*sum(SAEM.UYi*media1[-(1:p)]) + sum(SAEM.U*(media1[-(1:p)])^2)
    uyw = SAEM.UYZ - part1
    uw2 = SAEM.UZ2 - part2
    
    phi = as.vector(solve(uw2)%*%uyw)   # Update phi
    sigma2 = (uy2 - t(phi)%*%uyw - t(uyw)%*%phi + t(phi)%*%uw2%*%phi)/(n-p)
    sigma2 = as.numeric(sigma2)  # Update sigma2
    
    ciclos = ciclobeta(X,phi,SAEM.U,SAEM.UYi,SAEM.UZi)
    Ai = ciclos$Ai
    Bi = ciclos$Bi
    beta = as.vector(solve(Ai)%*%Bi)   # Update beta
    
    if (!fixed.nu){
    nu = optimize(f=Maxnu, lower=2, upper=200, maximum=TRUE, SAEMu=SAEM.U, SAEMlu=SAEM.lU)$maximum
    nu = as.numeric(nu)     # Update nu
    theta1 = c(beta, sigma2, phi, nu)
    } else { theta1 = c(beta, sigma2, phi) }
    
    # Stopping criteria
    criterio2 = sqrt((theta1/theta-1)%*%(theta1/theta-1))
    criterio = criterio2
    if (count==MaxIter){criterio = 1e-12}
    theta = theta1
    Theta = rbind(Theta, theta)
  } # End SAEM algorithm

  if (!quiet) setTxtProgressBar(pb, MaxIter)
  endtime   = Sys.time()
  timediffe = endtime-initime
  
  # Observed information matrix
  if (show_ep){
    score1   = GradientExp(SAEM.lU, SAEM.U, uy2, uyw, uw2, Ai, Bi, phi, beta, sigma2, nu, fixed.nu)
    hessian1 = HessianExp(SAEM.U, SAEM.UZi, SAEM.UYi, uy2, uyw, uw2, Ai, Bi, media1, X, phi, beta,
                           sigma2, nu, fixed.nu)
    ObsInfM  = score1%*%t(score1) - hessian1 - SAEM.delta
    variancias = diag(inversa(ObsInfM)) # variance error approximation
    se.app     = sqrt(variancias)       # standard error approximation
  }
  
  # Prediction
  if (!is.null(x_pred)){
    m = nrow(x_pred)
    meanDiff = SAEM.Y - X%*%beta
    media.pre = as.matrix(x_pred)%*%beta
    y_pred = matrix(0, ncol=1, nrow=m)
    for (k in 1:m){
      y_pred[k,] = media.pre[k] + 
        if(k==1){t(phi)%*%c(meanDiff[n:(n-p+k)])}else{0} +
        if(1<k & k<=p){t(phi)%*%c((y_pred[(k-1):1] - media.pre[(k-1):1]),(meanDiff[n:(n-p+k)]))}else{0} +
        if(k>p){t(phi)%*%c(y_pred[(k-1):(k-p)] - media.pre[(k-1):(k-p)])}else{0}
    }
    if (show_ep) resultados = list(beta=beta, sigma2=sigma2, phi=phi, nu=nu, theta=theta, SE=se.app,
                                    pred=y_pred, criteria=criterio2)
    else resultados = list(beta=beta, sigma2=sigma2, phi=phi, nu=nu, theta=theta, pred=y_pred, criteria=criterio2)
  } else {
    if (show_ep) resultados = list(beta=beta, sigma2=sigma2, phi=phi, nu=nu, theta=theta, SE=se.app, criteria=criterio2)
    else resultados = list(beta=beta, sigma2=sigma2, phi=phi, nu=nu, theta=theta, criteria=criterio2)
  }
  return (list(res=resultados, SAEMy=SAEM.Y, SAEMu=SAEM.U, iter=count, time=timediffe, Theta=Theta))
}
