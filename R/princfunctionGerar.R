#source('auxfunctions.r')

################gerando os dados
#rARCens
gerarARCens = function(n,beta,pit,sig2,x=rep(1,n),cens='left',pcens=.1) {
  phi = estphit(pit);p = length(phi)
  beta<- as.matrix(beta)
  x=as.matrix(x)
  erro = as.matrix(arima.sim(n,model=list(ar=phi),sd=sqrt(sig2)))
  resp<-x%*%beta+erro
  if (cens=='left') {
    cte<-as.numeric(quantile(resp,probs=pcens))
    cc<-(resp<cte)+0
  }
  else {
    cte<-as.numeric(quantile(resp,probs=1-pcens))
    cc<-(resp>cte)+0
  }
  y<-resp*(1-cc)+cte*cc
  return(data.frame(y=y,cc=cc))
}
