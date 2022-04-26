#source('princfunctionGerar.R')

rARCens = function(n,beta,pit,sig2=1,x=rep(1,n),cens='left',pcens=.1)
{
  if ((!is.numeric(n))|(length(n)!=1)|n<=0) stop("n must be a positive integer number")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x=as.matrix(x)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
  if (!is.numeric(beta)) stop("beta must be a numeric vector")
  if (!is.numeric(pit)) stop("pit must be a numeric vector")
  if (!is.numeric(sig2)) stop("sig2 must be a number")
  if (!is.numeric(pcens)) stop("pcens must be a number")
  ## Verify error at parameters specification

  #No data
  if( (length(x) == 0) | (length(beta) == 0) | (length(pit) == 0)) stop("All parameters must be provided.")

  #Validating if exists NA's
  if(sum(is.na(x)) > 0) stop("There are some NA values in x")

  #Validating dims data set
  if( length(beta) != ncol(as.matrix(x)) ) stop("the length of beta must the equal to the number of columns of x")
  if (nrow(x) !=n) stop("the number of rows of x must be equal to n")

  #Validating supports
  if(length(sig2)!=1) stop("sig2 must be a positive value")
  if(!is.numeric(pcens)) stop("pcens must be a real number in [0,1]")
  if(pcens > 1 | pcens < 0) stop("pcens must be a real number in [0,1]")
  if (cens!='left'& cens!='right') stop('cens must be left or right')
  if (any(pit>=1)|any(pit<=-1)) stop('pit must be a vector with elements in (-1,1)')
  #Load required libraries

  #Running the algorithm
  out = list()
  out$data = gerarARCens(n=n,beta=beta,pit=pit,sig2=sig2,x=x,cens=cens,pcens=pcens)
  phi = estphit(pit)
  out$param = c(beta,sig2,phi)
  return(out)
}
