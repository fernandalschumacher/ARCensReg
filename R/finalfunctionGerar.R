
rARCens = function(n, beta, pit, sig2=1, x=rep(1,n), cens='left', pcens=0.1, innov="norm", nu=NULL)
{
  if ((!is.numeric(n)) | (length(n)!=1) | n<=0) stop("n must be a positive integer number")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x = as.matrix(x)
  if (det(t(x)%*%x) == 0) stop("The columns of x must be linearly independent")
  if (!is.numeric(beta)) stop("beta must be a numeric vector")
  if (!is.numeric(pit)) stop("pit must be a numeric vector")
  if (!is.numeric(sig2)) stop("sig2 must be a number")
  if (!is.numeric(pcens)) stop("pcens must be a number")
  
  ## Verify error at parameters specification
  #No data
  if ((length(x) == 0) | (length(beta) == 0) | (length(pit) == 0)) stop("All parameters must be provided.")

  #Validating if exists NA's
  if (sum(is.na(x)) > 0) stop("There are some NA values in x")

  #Validating dims data set
  if (length(beta) != ncol(as.matrix(x))) stop("The length of beta must be equal to the number of columns of x")
  if (nrow(x) != n) stop("The number of rows of x must be equal to n")

  #Validating supports
  if (length(sig2)!=1) stop("sig2 must be a positive value")
  if (sig2 <= 0) stop("sig2 must be a positive value")
  if (!is.numeric(pcens)) stop("pcens must be a real number in [0,1]")
  if (pcens > 1 | pcens < 0) stop("pcens must be a real number in [0,1]")
  if (cens!='left' & cens!='right') stop('cens must be left or right')
  if (any(pit>=1) | any(pit<=-1)) stop('pit must be a vector with elements in (-1,1)')
  
  if (innov!='norm' & innov!='t') stop('innov must be norm or t')
  if (innov == 't'){
    if (is.null(nu)) stop('nu must be provided for Student-t innovations')
    if (!is.numeric(nu)) stop('nu must be a positive number (greater than 2)')
    if (nu <= 2) stop('nu must be a positive number (greater than 2)')
  }
  #Load required libraries

  #Running the algorithm
  out = list()
  phi = estphit(pit)
  out$data = gerarARCens(n=n, beta=beta, phi=phi, sig2=sig2, x=x, cens=cens, pcens=pcens,
                         innov=innov, nu=nu)
  if (innov == "norm") out$param = c(beta, sig2, phi)
  if (innov == "t") out$param = c(beta, sig2, phi, nu)
  
  class(out)  = 'rARCens'
  return(out)
}


################################
# Simulate a dataset #
################################

gerarARCens = function(n, beta, phi, sig2, x, cens, pcens, innov, nu) {
  p = length(phi)
  beta = as.matrix(beta)
  x = as.matrix(x)
  sigma = sqrt(sig2)
  
  if (innov == 'norm'){
    erro = as.matrix(arima.sim(n=n, model=list(ar=phi), sd=sigma))
  } else {
    erro = as.matrix(arima.sim(model=list(ar=phi), n=n, 
                               rand.gen=function(n,...){sigma*rt(n,nu)}))
  }
  resp = x%*%beta + erro
  
  if (pcens == 0){
    return (data.frame(y=resp, cc=rep(0, n)))
    
  } else {
    if (innov == 't'){ prob2 = ifelse(pcens*n/(n-p)>1, 1, pcens*n/(n-p)) }
    
    if (cens=='left') {
      if (innov == 'norm'){
        cte = as.numeric(quantile(resp, probs=pcens))
        cc = (resp<cte) + 0
      } else {
        cte = as.numeric(quantile(resp[-(1:p)], probs=prob2))
        cc = rep(0, n);   cc[-(1:p)] = (resp[-(1:p)]<cte) + 0
      }
      y = resp*(1-cc) + cte*cc
      LI = rep(-Inf, n)
      LS = rep(cte, n)
      if (innov == 't'){ LS[1:p] = rep(min(y[1:p], cte), p) }

    } else {
      
      if (innov == 'norm'){
        cte = as.numeric(quantile(resp, probs=1-pcens))
        cc = (resp>cte) + 0
      } else {
        cte = as.numeric(quantile(resp[-(1:p)], probs=1-prob2))
        cc = rep(0, n);  cc[-(1:p)] = (resp[-(1:p)]>cte) + 0
      }
      y = resp*(1-cc) + cte*cc
      LI = rep(cte, n)
      LS = rep(Inf, n)
      if (innov == 't'){ LI[1:p] = rep(max(y[1:p], cte), p) }
    }
    return(data.frame(y=y, cc=cc, lcl=LI, ucl=LS))    
  } # End else
}
