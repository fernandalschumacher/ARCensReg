#source('princfunction.R')

InfDiagys=function(theta,yest,yyest,x,k=3,plots=T,indpar=rep(1,length(theta))){
  n = length(yest)
  hes = Qt(theta,yest,yyest,x)
  if (plots) par(mfrow=c(1,1),mar=c(4, 4, 3, 2) + 0.1)
  ############################### y(w) = y+w
  delta = deltaw(theta,yest,x)
  npar = length(theta)
  Q = solve(-hes)
  if (sum(indpar)<npar) {
    indt = (1:npar)[indpar==0]
    b22 = matrix(0,ncol=npar,nrow=npar)
    b22[indt,indt] = solve(hes[indt,indt])
    Q = Q + b22
  }
  matF = 2*t(delta)%*%Q%*%delta
  ci = diag(matF)
  auto = eigen(matF)
  autoval = auto$values[abs(auto$values)>.0001]
  autovalp = autoval/sum(autoval)
  autovec = auto$vectors[,abs(auto$values)>.0001]
  Mautovalp = matrix(autovalp,n,length(autovalp),byrow=T)
  M0y = apply(Mautovalp*autovec^2,1,sum)
  bm = mean(M0y)+k*sd(M0y)
  if (plots) {
    plot(M0y,pch=ifelse(M0y>bm,16,1),ylab='Response perturbation')
    abline(h=bm,lty='dashed')
    if (length(which(M0y>bm))>0) text(which(M0y>bm)+length(M0y)/10,M0y[M0y>bm],labels=which(M0y>bm))
  }
  return(M0y)
}
InfDiagSigs=function(theta,yest,yyest,x,k=3,plots=T,indpar=rep(1,length(theta))){
  n = length(yest)
  hes = Qt(theta,yest,yyest,x)
  ############################### Sigmatil=Sigma D(w)
  delta = NULL
  for (i in 1:n) delta = cbind(delta,deltaSigi(theta,yest,yyest,x,i))
  npar = length(theta)
  Q = solve(-hes)
  if (sum(indpar)<npar) {
    indt = (1:npar)[indpar==0]
    b22 = matrix(0,ncol=npar,nrow=npar)
    b22[indt,indt] = solve(hes[indt,indt])
    Q = Q + b22
  }
  matF = 2*t(delta)%*%Q%*%delta
  auto = eigen(matF)
  autoval = auto$values[abs(auto$values)>.0001]
  autovalp = autoval/sum(autoval)
  autovec = auto$vectors[,abs(auto$values)>.0001]
  Mautovalp = matrix(autovalp,n,length(autovalp),byrow=T)
  M0sig = apply(Mautovalp*autovec^2,1,sum)
  bm = mean(M0sig)+k*sd(M0sig)
  if (plots) {
    plot(M0sig,pch=ifelse(M0sig>bm,16,1),ylab='Scale matrix perturbation')#,type='h')
    abline(h=bm,lty='dashed')
    if (length(which(M0sig>bm))>0) text(which(M0sig>bm)+length(M0sig)/10,M0sig[M0sig>bm],labels=which(M0sig>bm))
  }
  ###
  return(M0sig)
}

InfDiagxps = function(theta,yest,yyest,x,k=3,indp=rep(1,ncol(x)),plots=T,indpar=rep(1,length(theta))){
  n = length(yest)
  hes = Qt(theta,yest,yyest,x)
  ###############################
  delta = deltaxp(theta,yest,x,indp)
  npar = length(theta)
  Q = solve(-hes)
  if (sum(indpar)<npar) {
    indt = (1:npar)[indpar==0]
    b22 = matrix(0,ncol=npar,nrow=npar)
    b22[indt,indt] = solve(hes[indt,indt])
    Q = Q + b22
  }
  matF = 2*t(delta)%*%Q%*%delta
  auto = eigen(matF)
  autoval = auto$values[abs(auto$values)>.0001]
  autovalp = autoval/sum(autoval)
  autovec = auto$vectors[,abs(auto$values)>.0001]
  Mautovalp = matrix(autovalp,n,length(autovalp),byrow=T)
  M0xp = apply(Mautovalp*autovec^2,1,sum)
  bm = mean(M0xp)+k*sd(M0xp)
  if (plots) {
    plot(M0xp,pch=ifelse(M0xp>bm,16,1),ylab='x perturbation')
    abline(h=bm,lty='dashed')
    if (length(which(M0xp>bm))>0) text(which(M0xp>bm)+length(M0xp)/20,M0xp[M0xp>bm],labels=which(M0xp>bm))
  }
  return(M0xp)
}
