
InfDiag = function(object, k=3, indpar=rep(1,length(object$theta)), indcolx=rep(1,ncol(object$x)), 
                   perturbation='y')
{
  ## Verify error at parameters specification
  if (!is(object, "ARpCRM")) stop("object must be of class 'ARpCRM'")
  if (!is.numeric(k)) stop("k must be a number")
  
  theta = object$theta
  yest  = object$yest
  yyest = object$yyest 
  x     = object$x

  #Validating if exists NA's
  if(perturbation=='x') if(sum(indcolx)==0) stop("indcolx must have at least one element equal to 1")
  if(sum(indpar)==0) stop("indpar must have at least one element equal to 1")

  #Validating dims data set
  if (length(indpar)!= length(theta)) stop("indpar does not have the same length than theta")
  if(perturbation=='x') if (length(indcolx)!= ncol(x)) stop("the lenght of indcolx must be equal to the number of columns of x")

  #Validating supports
  if(length(k)!=1) stop("k must be an integer value")
  if(!any(perturbation==c('y','Sigma','x'))) stop("perturbation must be 'y', 'Sigma' or 'x'")
  if(perturbation=='x') if(any(!(indcolx %in% c(0,1)))) stop("indcolx must a vector with elements 0 or 1")
  if(any(!(indpar %in% c(0,1)))) stop("indpar must a vector with elements 0 or 1")

  #Running the algorithm
  if (perturbation=='y') M0 = InfDiagys(theta=theta, yest=yest, yyest=yyest, x=x, 
                                        k=k, indpar=indpar)
  if (perturbation=='Sigma') M0 = InfDiagSigs(theta=theta, yest=yest, yyest=yyest, 
                                              x=x, k=k, indpar=indpar)
  if (perturbation=='x') M0 = InfDiagxps(theta=theta, yest=yest, yyest=yyest, x=x, 
                                         k=k, indp=indcolx, indpar=indpar)

  bench = mean(M0) + k*sd(M0)
  cat('\n')
  cat("Perturbation scheme:",perturbation)
  cat('\n')
  cat('Benchmark =',round(bench,3))
  cat('\n')
  detecpoints= which(M0>bench)
  if (length(detecpoints)>0) cat('Detected points:',detecpoints,'\n')
  else cat('Detected points:',0,'\n')
  
  M1 = list(M0=M0, perturbation=perturbation, benchmark=bench)
  class(M1) = "DiagARpCRM"
  return(M1)
}


#' @export
plot.DiagARpCRM = function(x, ...) {
  M0 = x$M0
  n = nrow(matrix(M0))
  M0y2 = data.frame(M0 = M0) 
  rownames(M0y2) = 1:n
  bm = x$benchmark
  if (x$perturbation == "y") text = 'Response perturbation'
  if (x$perturbation == "Sigma") text = 'Scale matrix perturbation'
  if (x$perturbation == "x") text = 'x perturbation'
  
  ggplot(M0y2, aes(x=1:n, y=M0)) + geom_point(shape=ifelse(M0>bm,16,21), size=1.9) +
    geom_hline(yintercept=bm, linetype="dashed") +
    geom_text(aes(label=ifelse(M0>bm, rownames(M0y2), '')), vjust=1.5) +
    labs(x="Index", y=text) + theme_bw()
}
