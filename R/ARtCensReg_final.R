globalVariables(c("z","lag","..density.."))

ARtCensReg = function(cc, lcl=NULL, ucl=NULL, y, x, p=1, x_pred=NULL, tol=0.0001, M=10,
                       perc=0.25, MaxIter=400, pc=0.18, nufix=NULL, show_se=TRUE, quiet=FALSE){
  m = length(y)
  
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x = as.matrix(x)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
  
  ## Verify error at parameters specification
  #No data
  if ((length(x) == 0) | (length(y) == 0) | (length(cc) == 0)) stop("All parameters must be provided")

  #Validating if exists NA's
  if (sum(cc[1:p]) > 0) stop("The first p values in y must be completely observed")
  if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")
  if (sum(is.na(x)) > 0) stop("There are some NA values in x")
  if (sum(is.na(cc)) > 0) stop("There are some NA values in cc")
  miss = which(is.na(y))
  if (sum(cc[miss]) != length(miss)) stop ("NA values in y must be specified through arguments cc, lcl, and ucl")
  
  #Validating dims data set
  if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if (ncol(as.matrix(cc)) > 1) stop("cc must have just one column")
  if (nrow(as.matrix(x)) != m) stop("x does not have the same number of lines than y")
  if (length(cc) != m) stop("cc does not have the same length than y")

  if (sum(cc) > 0){
    if (is.null(lcl) | is.null(ucl)) stop("lcl and ucl must be provided for censored data")
    if (!is.numeric(lcl) | !is.numeric(ucl)) stop("lcl and ucl must be numeric vectors")
    if (length(miss)>0){
      censor = (cc==1 & !is.na(y))
      if (any(is.infinite(lcl[censor])) & any(is.infinite(ucl[censor]))) stop("lcl or ucl must be finite for censored data")
    } else {
      if (any(is.infinite(lcl[cc==1])) & any(is.infinite(ucl[cc==1]))) stop("lcl or ucl must be finite for censored data") 
    }
    if (length(lcl) != m) stop("lcl does not have the same length than y")
    if (length(ucl) != m) stop("ucl does not have the same length than y")
    if (ncol(as.matrix(lcl)) > 1) stop("lcl must have just one column")
    if (ncol(as.matrix(ucl)) > 1) stop("ucl must have just one column")
    if (sum(is.na(lcl))>0 | sum(is.na(ucl))>0) stop("There are some NA values in lcl or ucl")
    if (!all(lcl[cc==1]<ucl[cc==1])) stop ("lcl must be smaller than ucl")
  }
  
  if (!is.null(x_pred)) {
    x_pred = as.matrix(x_pred)
    if (ncol(x_pred)!=ncol(as.matrix(x))) stop("x_pred must have the same number of columns than x")
    if (sum(is.na(x_pred))>0) stop("There are some NA values in x_pred")
    if (!is.numeric(x_pred)) stop("x_pred must be a numeric matrix")
  }
  
  #Validating supports
  if (!is.null(nufix)){ 
    if (length(c(nufix)) != 1) stop("nufix must be a positive value or 'NULL'")
    if (!is.numeric(nufix)) stop("nufix must be a positive value")
    if (nufix <= 2) stop("nufix must be a positive value (greater than 2)")
  }
  if (length(p) != 1) stop("p must be a positive integer value")
  if (!is.numeric(p)) stop("p must be a positive integer value")
  if (p!=round(p) | p<=0) stop("p must be a positive integer value")
  if (tol <= 0) stop("tolerance must be a positive value (suggested to be small)")
  if (!is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
  if (MaxIter<=0 | MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  if (!is.numeric(M)) stop("M must be a positive integer value")
  if (M<=1 | M%%1!=0) stop("M must be a positive integer value (greater than 1)")
  if (!is.numeric(pc)) stop("pc must be a real number in [0,1]")
  if (pc>1 | pc<0) stop("pc must be a real number in [0,1]")
  if (!is.numeric(perc)) stop("perc must be a real number in [0,1)")
  if (perc>=1 | perc<0) stop("perc must be a real number in [0,1)")
  if (!is.logical(show_se)) stop("show_se must be TRUE or FALSE")
  if (!is.logical(quiet)) stop("quiet must be TRUE or FALSE")
  
  #Running the algorithm
  if (!quiet) {
    cat('\n')
    call <- match.call()
    cat("Call:\n")
    print(call)
    cat('\n')
  }
  out = suppressWarnings(SAEM_temporalT(cc, lcl, ucl, y, x, p, x_pred, tol, M, perc, MaxIter, pc, nufix, show_se, quiet))
  q = ncol(x)
  if (is.null(nufix)){ lab = numeric(p+q+2); lab[p+q+2] = 'nu' } else { lab = numeric(p+q+1) }
  if (sum(abs(x[,1])) == nrow(x)){ for (i in 1:q) lab[i] = paste('beta',i-1,sep='') 
  } else { for (i in 1:q) lab[i] = paste('beta',i,sep='') }
  lab[q+1] = 'sigma2'
  for (i in ((q+2):(p+q+1))) lab[i] = paste('phi',i-q-1,sep='')
  if (show_se) {
    tab = round(rbind(out$res$theta, out$res$SE),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(out$res$theta),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  obj.out = out$res
  obj.out$call = match.call()
  obj.out$tab = tab
  if (sum(cc) == 0){ cens = "no-censoring" 
  } else {
    if (sum(cc) == length(miss)){ cens = "missing"
    } else {
      if (all(is.infinite(lcl)) & any(is.finite(ucl))){ cens = "left" }
      if (all(is.infinite(ucl)) & any(is.finite(lcl))){ cens = "right" }
      if (any(is.finite(ucl-lcl))){ cens = "interval" }
    }
  }
  obj.out$cens = cens
  obj.out$nmiss = length(miss)
  obj.out$ncens = sum(cc)
  obj.out$converge = (out$res$iter < MaxIter)
  obj.out$MaxIter = MaxIter
  obj.out$M = M
  obj.out$pc = pc
  obj.out$time = out$time
  #plot
  obj.out$plot$cpl = pc*MaxIter
  obj.out$plot$npar = length(out$res$theta)
  obj.out$plot$labels = list()
  if (sum(abs(x[,1]))==nrow(x)) { for(i in 1:q){obj.out$plot$labels[[i]] = bquote(beta[.(i-1)])} 
    } else { for(i in 1:q){obj.out$plot$labels[[i]] = bquote(beta[.(i)])} }
  obj.out$plot$labels[[q+1]] = bquote(sigma^2)
  for(i in 1:p){obj.out$plot$labels[[i+q+1]] = bquote(phi[.(i)])}
  if (obj.out$plot$npar == (p+q+2)) obj.out$plot$labels[[p+q+2]] = bquote(nu)
  obj.out$plot$Theta = out$Theta
  #class
  class(obj.out)  = 'ARtpCRM'
  invisible(obj.out)
}


#' @export
print.ARtpCRM = function(x, ...){
  cat('---------------------------------------------------\n')
  cat('  Censored Linear Regression Model with AR Errors \n')
  cat('---------------------------------------------------\n')
  cat("Call:\n")
  print(x$call)
  cat('\n')
  cat('Estimated parameters:\n')
  print(x$tab)
  cat('\n')
  cat('Details:\n')
  cat('Type of censoring:', x$cens, '\n')
  if (x$ncens > 0){ cat('Number of missing values:', x$nmiss, '\n') }
  cat("Convergence reached?:", x$converge, '\n')
  cat('Iterations:', x$iter,"/",x$MaxIter, '\n')
  cat('MC sample:', x$M, '\n')
  cat('Cut point:', x$pc, '\n')
  cat("Processing time:", x$time, units(x$time), '\n')
}


#' @export
summary.ARtpCRM = function(object, ...){
  cat('---------------------------------------------------\n')
  cat('  Censored Linear Regression Model with AR Errors \n')
  cat('---------------------------------------------------\n')
  cat("Call:\n")
  print(object$call)
  cat('\n')
  cat('Estimated parameters:\n')
  print(object$tab)
  cat('\n')
  cat('Details:\n')
  cat('Type of censoring:', object$cens, '\n')
  if (object$ncens > 0) { cat('Number of missing values:', object$nmiss, '\n') }
  cat("Convergence reached?:", object$converge, '\n')
  cat('Iterations:', object$iter,"/",object$MaxIter, '\n')
  cat('MC sample:', object$M, '\n')
  cat('Cut point:', object$pc, '\n')
  cat("Processing time:", object$time, units(object$time), '\n')
}


#' @export
plot.ARtpCRM = function(x, ...) {
  count = x$iter
  npar  = x$plot$npar
  label = x$plot$labels
  myplot = vector("list", npar)
  
  for (i in 1:npar){
    data1 = data.frame(z=x$plot$Theta[,i])
    myplot[[i]] = ggplot(data1, aes(x=seq(0,count), y=z)) + geom_line() +
      geom_vline(xintercept=x$plot$cpl, color="red", linetype="twodash") +
      labs(x="Iteration", y=label[[i]]) + theme_bw()
  }
  nrows = ifelse(npar%%3==0, npar%/%3, (npar%/%3)+1)
  grid.arrange(grobs=myplot, nrow=nrows, ncol=3)
}


#' @export
residuals.ARtpCRM = function(object, ...) {
  
  x = object$x
  p = length(object$phi)
  m = nrow(x)
  residuals = numeric(m)
  residuals[1:p] = 0
  
  res = object$yest - x%*%object$beta
  for (i in (p+1):m) residuals[i] = res[i] - sum(object$phi*res[(i-1):(i-p)])
  #
  quant = residuals/sqrt(object$sigma2)
  quant = qnorm(pt(quant, object$nu))
  
  resid = list(residuals=residuals, quantile.resid=quant)
  class(resid) = "residARpCRM"
  return(resid)
}
