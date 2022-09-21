
ARCensReg = function(cc, lcl=NULL, ucl=NULL, y, x, p=1, M=10, perc=0.25, MaxIter=400, 
                     pc=0.18, tol=0.0001, show_se=TRUE, quiet=FALSE){
  m = length(y)
  
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x = as.matrix(x)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
  
  ## Verify error at parameters specification
  #No data
  if ( (length(x) == 0) | (length(y) == 0) | (length(cc) == 0)) stop("All parameters must be provided")
  
  #Validating if exists NA's
  if (sum(cc%in%c(0,1)) < length(cc)) stop("The elements of the vector cc must be 0 or 1")
  if (sum(is.na(x)) > 0) stop("There are some NA values in x")
  if (sum(is.na(cc)) > 0) stop("There are some NA values in cc")
  miss = which(is.na(y))
  if (length(miss)>0) { if (sum(cc[miss]) != length(miss)) stop ("NA values in y must be specified through arguments cc, lcl, and ucl")
  } else { miss = NULL }
  
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
      if (any(is.infinite(lcl[censor]) & is.infinite(ucl[censor]))) stop("lcl or ucl must be finite for censored data")
    } else { 
      if (any(is.infinite(lcl[cc==1]) & is.infinite(ucl[cc==1]))) stop("lcl or ucl must be finite for censored data") 
    }
    if (length(lcl) != m) stop("lcl does not have the same length than y")
    if (length(ucl) != m) stop("ucl does not have the same length than y")
    if (ncol(as.matrix(lcl)) > 1) stop("lcl must have just one column")
    if (ncol(as.matrix(ucl)) > 1) stop("ucl must have just one column")
    if (sum(is.na(lcl))>0 | sum(is.na(ucl))>0) stop("There are some NA values in lcl or ucl")
    if (!all(lcl[cc==1]<ucl[cc==1])) stop ("lcl must be smaller than ucl")
  }
  
  #Validating supports
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

  #Load required libraries

  #Running the algorithm
  if (!quiet) {
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  }
  out = suppressWarnings(SAEM(cc, lcl, ucl, y, x, p, M, perc, MaxIter, pc, miss, tol, show_se, quiet))
  l = ncol(x)
  lab = numeric(p + l + 1)
  if (sum(abs(x[,1])) == nrow(x)){ for (i in 1:ncol(x)) lab[i] = paste('beta',i-1,sep='') 
  } else { for (i in 1:ncol(x)) lab[i] = paste('beta',i,sep='') }
  lab[l+1] = 'sigma2'
  for (i in ((l+2):length(lab))) lab[i] = paste('phi',i-l-1,sep='')
  if (show_se) {
    tab = round(rbind(out$res$theta, out$res$SE), 4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(out$res$theta), 4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  critFin = c(out$res$loglik, out$res$AIC, out$res$BIC, out$res$AICcorr)
  critFin = round(t(as.matrix(critFin)), digits=3)
  dimnames(critFin) = list(c("Value"),c("Loglik", "AIC", "BIC","AICcorr"))
  
  obj.out = out$res
  obj.out$call = match.call()
  obj.out$tab = tab
  obj.out$critFin = critFin
  if (sum(cc) == 0){ cens = "no censoring" 
  } else {
    if (sum(cc) == length(miss)){ cens = "missing"
    } else {
      if (all(is.infinite(lcl)) & any(is.finite(ucl))){ cens = "left" }
      if (all(is.infinite(ucl)) & any(is.finite(lcl))){ cens = "right" }
      if (any(is.finite(ucl-lcl))){ cens = "interval" }
    }
  }
  obj.out$cens = cens
  obj.out$nmiss = ifelse(is.null(miss),0,length(miss))
  obj.out$ncens = sum(cc)
  if (sum(cc)>0){
    obj.out$converge = (out$iter < MaxIter)
    obj.out$MaxIter = MaxIter
    obj.out$M = M
    obj.out$pc = pc
  }
  obj.out$time = out$time
  #plot
  if (sum(cc) > 0){
    obj.out$plot$cpl = pc*MaxIter
    obj.out$plot$npar = l+1+p
    obj.out$plot$labels = list()
    if (sum(abs(x[,1]))==nrow(x)) { for(i in 1:l){obj.out$plot$labels[[i]] = bquote(beta[.(i-1)])} 
    } else { for(i in 1:l){obj.out$plot$labels[[i]] = bquote(beta[.(i)])} }
    obj.out$plot$labels[[l+1]] = bquote(sigma^2)
    for(i in 1:p){obj.out$plot$labels[[i+l+1]] = bquote(phi[.(i)])}
    obj.out$plot$Theta = out$Theta
  }
  #class
  class(obj.out)  = 'ARpCRM'
  invisible(obj.out)
}


#' @export
print.ARpCRM = function(x, ...){
  cat('---------------------------------------------------\n')
  cat('  Censored Linear Regression Model with AR Errors \n')
  cat('---------------------------------------------------\n')
  cat("Call:\n")
  print(x$call)
  cat('\n')
  cat('Estimated parameters:\n')
  print(x$tab)
  cat('\n')
  cat('Model selection criteria:\n')
  print(x$critFin)
  cat('\n')
  cat('Details:\n')
  cat('Type of censoring:', x$cens, '\n')
  if (x$ncens > 0) {
    cat('Number of missing values:', x$nmiss, '\n')
    cat("Convergence reached?:", x$converge, '\n')
    cat('Iterations:', x$iter,"/",x$MaxIter, '\n')
    cat('MC sample:', x$M, '\n')
    cat('Cut point:', x$pc, '\n')
  }
  cat("Processing time:", x$time, units(x$time), '\n')
}


#' @export
summary.ARpCRM = function(object, ...){
  cat('---------------------------------------------------\n')
  cat('  Censored Linear Regression Model with AR Errors \n')
  cat('---------------------------------------------------\n')
  cat("Call:\n")
  print(object$call)
  cat('\n')
  cat('Estimated parameters:\n')
  print(object$tab)
  cat('\n')
  cat('Model selection criteria:\n')
  print(object$critFin)
  cat('\n')
  cat('Details:\n')
  cat('Type of censoring:', object$cens, '\n')
  if (object$ncens>0) {
    cat('Number of missing values:', object$nmiss, '\n')
    cat("Convergence reached?:", object$converge, '\n')
    cat('Iterations:', object$iter,"/",object$MaxIter, '\n')
    cat('MC sample:', object$M, '\n')
    cat('Cut point:', object$pc, '\n')
  }
  cat("Processing time:", object$time, units(object$time), '\n')
}


#' @export
plot.ARpCRM = function(x, ...) {
  if (x$ncens == 0) stop("plot only defined for cases with censoring")
  
  count = x$iter
  npar  = x$plot$npar
  label = x$plot$labels
  myplot = vector("list", npar)
  
  for (i in 1:npar){
    data1 = data.frame(z=x$plot$Theta[,i])
    myplot[[i]] = ggplot(data1, aes(x=seq(1,count), y=z)) + geom_line() +
      geom_vline(xintercept=x$plot$cpl, color="red", linetype="twodash") +
      labs(x="Iteration", y=label[[i]]) + theme_bw()
  }
  nrows = ifelse(npar%%3==0, npar%/%3, (npar%/%3)+1)
  grid.arrange(grobs=myplot, nrow=nrows, ncol=3)
}


#' @export
predict.ARpCRM = function(object, x_pred, ...){

  # validation
  x_pred = as.matrix(x_pred)
  if (ncol(x_pred)!=ncol(as.matrix(object$x))) stop("x_pred must have the same number of columns than x")
  if (sum(is.na(x_pred))>0) stop("There are some NA values in x_pred")
  if (!is.numeric(x_pred)) stop("x_pred must be a numeric matrix")
  
  y = object$yest
  x = object$x
  beta1 = object$beta
  sigmae = object$sigma2
  phi1 = object$phi
  m = length(c(y))
  h = nrow(x_pred)
  sig_pred = MatArp(phi1, m+h)*sigmae
  pred = x_pred%*%beta1 + sig_pred[m+1:h, 1:m]%*%solve(sig_pred[1:m, 1:m])%*%(y - x%*%beta1)
  
  return (pred)
}


#' @export
residuals.ARpCRM = function(object, ...) {
  
  x = object$x
  p = length(object$phi)
  m = nrow(x)
  residuals = numeric(m)
  residuals[1:p] = 0
  
  res = object$yest - x%*%object$beta
  for (i in (p+1):m) residuals[i] = res[i] - sum(object$phi*res[(i-1):(i-p)])
  
  resid = list(residuals=residuals[-(1:p)], quantile.resid=(residuals[-(1:p)])/sqrt(object$sigma2))
  class(resid) = "residARpCRM"
  return(resid)
}


#' @export
plot.residARpCRM = function(x, ...) {
  
  resid = data.frame(resid=x$quantile.resid)
  replot = list(4)
  m = nrow(resid)
  #
  replot[[1]] = ggplot(resid, aes(x=seq(1,m),y=resid)) + geom_line() + labs(x="Time", y="Quantile Residual") + 
    geom_hline(yintercept=c(-2,0,2), color="red", linetype="twodash") + theme_bw()
  #
  replot[[2]] = ggplot(resid, aes(sample=resid)) + stat_qq_band(distribution="norm", identity=TRUE) + 
    stat_qq_line(distribution="norm", color="red", linetype="twodash", identity=TRUE) +
    stat_qq_point(distribution="norm", identity=TRUE, size=1, alpha=0.5) + 
    labs(x="Theoretical Quantiles", y="Sample Quantiles") + theme_bw()
  #
  replot[[3]] = ggplot(resid, aes(x=resid)) + geom_histogram(aes(y=..density..), fill="grey", color="black", bins=15) +
    stat_function(fun=dnorm, col="red", linetype="twodash") + labs(x="Quantile Residual",y="Density") + theme_bw()
  #
  bacfdf = with(acf(resid, plot=FALSE), data.frame(lag, acf))
  replot[[4]] = ggplot(data=bacfdf, aes(x=lag, y=acf)) + geom_hline(aes(yintercept=0)) + theme_bw() +
    geom_segment(aes(xend=lag, yend=0)) + labs(x="Lag", y="ACF") +
    geom_hline(yintercept=c(qnorm(0.975)/sqrt(m),-qnorm(0.975)/sqrt(m)), colour="red", linetype="twodash")
  #
  grid.arrange(grobs=replot, widths=c(1, 1, 1), layout_matrix = rbind(c(1, 1, 1), c(2, 3, 4)))
}
