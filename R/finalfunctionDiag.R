
InfDiag = function(theta, yest, yyest, x, k=3, plots=T, indpar=rep(1,length(theta)),
                   perturbation='y', indcolx=rep(1,ncol(x)))
{
  n = length(yest)

  if ((!is.numeric(yest))|any(is.na(yest))) stop("yest must be a numeric vector")
  if (!is.numeric(yyest)|any(is.na(yyest))|any(dim(yyest)!=c(n,n))) stop("yyest must be a numeric matrix of dimension n x n")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x=as.matrix(x)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
  if (!is.numeric(k)) stop("k must be a number")
  if (!is.numeric(theta)) stop("theta must be a numeric vector")
  ## Verify error at parameters specification

  #No data
  if( (length(x) == 0) | (length(yest) == 0) | (length(yyest) == 0)| (length(theta) == 0)) stop("All parameters must be provided.")

  #Validating if exists NA's
  if(sum(is.na(x)) > 0) stop("There are some NA values in x")
  if(perturbation=='x') if(sum(indcolx)==0) stop("indcolx must have at least one element equal to 1")
  if(sum(indpar)==0) stop("indpar must have at least one element equal to 1")

  #Validating dims data set
  if (ncol(as.matrix(yest)) > 1) stop("yest must have just one column")
  if( length(yest) != nrow(as.matrix(x)) ) stop("x does not have the same number of lines than yest")
  if( nrow(as.matrix(yyest)) != nrow(as.matrix(x)) ) stop("x does not have the same number of lines than yyest")
  if (length(indpar)!= length(theta)) stop("indpar does not have the same length than theta")
  if(perturbation=='x') if (length(indcolx)!= ncol(x)) stop("the lenght of indcolx must be equal to the number of columns of x")

  #Validating supports
  if(length(k)!=1) stop("k must be an integer value")
  if(!is.logical(plots)) stop("plots must be TRUE or FALSE")
  if(!any(perturbation==c('y','Sigma','x'))) stop("perturbation must be 'y', 'Sigma' or 'x'")
  if(perturbation=='x') if(any(!(indcolx %in% c(0,1)))) stop("indcolx must a vector with elements 0 or 1")
  if(any(!(indpar %in% c(0,1)))) stop("indpar must a vector with elements 0 or 1")
  #Load required libraries

  #Running the algorithm
  if (perturbation=='y') M0 = InfDiagys(theta=theta, yest=yest, yyest=yyest, x=x, k=k,
                                        plots=plots, indpar=indpar)
  if (perturbation=='Sigma') M0 = InfDiagSigs(theta=theta, yest=yest, yyest=yyest, x=x, k=k,
                                              plots=plots, indpar=indpar)
  if (perturbation=='x') M0 = InfDiagxps(theta=theta, yest=yest, yyest=yyest, x=x, k=k,
                                    plots=plots, indp=indcolx, indpar=indpar)

  cat('\n')
  cat("Perturbation scheme:",perturbation)
  cat('\n')
  cat('Benchmark =',round(mean(M0)+k*sd(M0),3))
  cat('\n')
  detecpoints= which(M0>mean(M0)+k*sd(M0))
  if (length(detecpoints)>0) cat('Detected points:',detecpoints,'\n')
  else cat('Detected points:',0,'\n')
  return(M0)
}
