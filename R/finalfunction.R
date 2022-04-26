#source('princfunction.R')

ARCensReg = function(cc,y,x,p=1,cens='left',x_pred=NULL,miss=NULL,tol=0.0001,
                     M=10,perc=0.25,MaxIter=400,pc=0.18,show_se=TRUE, quiet=FALSE)
{
  m = length(y)

  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x=as.matrix(x)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
  ## Verify error at parameters specification

  if (cens!='left'& cens!='right') stop('cens must be left or right')
  #No data
  if( (length(x) == 0) | (length(y) == 0) | (length(cc) == 0)) stop("All parameters must be provided.")

  #Validating if exists NA's
  if (length(miss)>0) {
    cc[miss] = 0
    if (sum(is.na(y[-miss]))>0) stop("NA values in y must be specified in the argument miss")
    if (sum(is.na(y[miss]))!=length(miss)) warning("y[miss] is not NA")
  }
  else {
    if (sum(is.na(y))>0) stop("NA values in y must be specified in the argument miss")
  }
  if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")
  if(sum(is.na(x)) > 0) stop("There are some NA values in x")
  if (sum(is.na(cc)) > 0) stop("There are some NA values in cc")


  #Validating dims data set
  if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if (ncol(as.matrix(cc)) > 1) stop("cc must have just one column")
  if( length(y) != nrow(as.matrix(x)) ) stop("x does not have the same number of lines than y")
  if( length(y) != length(cc) ) stop("cc does not have the same length than y")

  if (!is.null(x_pred)) {
    x_pred = as.matrix(x_pred)
    if (ncol(x_pred)!=ncol(as.matrix(x))) stop("x_pred must have the same number of columns than x")
    if (sum(is.na(x_pred))>0) stop("There are some NA values in x_pred")
    if (!is.numeric(x_pred)) stop("x_pred must be a numeric matrix")
  }

  if (sum(miss %in% 1:m)<length(miss) ) stop("miss must indicate the index of the missing data on y")

  #Validating supports
  if(length(p)!=1) stop("p must be a positive integer value")
  if(!is.numeric(p)) stop("p must be a positive integer value")
  if(p!= round(p)|p<=0) stop("p must be a positive integer value")
  if(tol <= 0) stop("tolerance must be a positive value (suggested to be small)")
  if(!is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
  if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  if(!is.numeric(M)) stop("M must be a positive integer value")
  if(M <= 1 |M%%1!=0) stop("M must be a positive integer value (greater than 1)")
  if(!is.numeric(pc)) stop("pc must be a real number in [0,1]")
  if(pc > 1 | pc < 0) stop("pc must be a real number in [0,1]")
  if(!is.numeric(perc)) stop("perc must be a real number in [0,1)")
  if(perc >= 1 | perc < 0) stop("perc must be a real number in [0,1)")
  if(!is.logical(show_se)) stop("show_se must be TRUE or FALSE.")

  #Load required libraries

  #Running the algorithm
  if (!quiet) {
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  }
  out <-suppressWarnings(SAEM(cc,y,x,p,M=M,cens=cens,perc=perc,MaxIter=MaxIter,
                              pc = pc,x_pred = x_pred,miss = miss,tol=tol,
                              show_ep=show_se, quiet = quiet))
  l = ncol(x)
  lab = numeric(p +l +1)
  for (i in 1:ncol(x)) lab[i] = paste('beta',i-1,sep='')
  lab[l+1] = 'sigma2'
  for (i in ((l+2):length(lab))) lab[i] = paste('phi',i-l-1,sep='')
  if (show_se) {
    tab = round(rbind(out$theta,out$ep),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  }
  else {
    tab = round(rbind(out$theta),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  critFin <- c(out$loglik, out$AIC, out$BIC, out$AICcorr)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","AICcorr"))


  if (!is.null(x_pred)) res = list(beta = out$beta,sigma2= out$sigmae,phi = out$phi1,pi1=out$pi1,theta =out$theta,SE=out$ep,
                                   loglik=out$loglik,AIC=out$AIC,BIC=out$BIC,AICcorr=out$AICcorr,time = out$timediff,pred=out$pred,criteria = out$criteria)
  else res = list(beta = out$beta,sigma= out$sigmae,phi = out$phi1,pi1=out$pi1,theta =out$theta,SE=out$ep,
            loglik=out$loglik,AIC=out$AIC,BIC=out$BIC,AICcorr=out$AICcorr,time = out$timediff,criteria = out$criteria)
  if (sum(cc)==0) {
    obj.out = list(res = res)
  }
  else obj.out = list(res = res,yest=out$yest,yyest=out$yyest,iter = out$iter)#,conv=out$Theta
  obj.out$call <- match.call()
  obj.out$tab <- tab
  obj.out$critFin <- critFin
  obj.out$cens <- cens
  obj.out$nmiss<- ifelse(is.null(miss),0,length(miss))
  obj.out$ncens <- sum(cc)
  obj.out$converge <- (out$iter < MaxIter)
  obj.out$MaxIter <- MaxIter
  obj.out$M <- M
  obj.out$pc <- pc
  obj.out$timediff <- out$timediff
  #plot
  obj.out$plot$cpl = pc*MaxIter
  obj.out$plot$npar   = l+1+p
  obj.out$plot$labels = list()
  for(i in 1:l){obj.out$plot$labels[[i]] = bquote(beta[.(i-1)])}
  obj.out$plot$labels[[l+1]] = bquote(sigma^2)
  for(i in 1:p){obj.out$plot$labels[[i+l+1]] = bquote(phi[.(i)])}
  obj.out$plot$Theta <-out$Theta
  #class
  class(obj.out)  = 'ARpCRM' #ifelse(sum(cc)==0,'ARp-LRM','ARp-CRM')
  invisible(obj.out)
}


print.ARpCRM <- function(x, ...){
  cat('\n')
  cat("Call:\n")
  print(x$call)
  cat('\n')
  cat('---------------------------------------------------\n')
  cat('  Censored Linear Regression Model with AR Errors \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  print(x$tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  print(x$critFin)
  cat('\n')
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat('Type of censoring =',x$cens)
  cat('\n')
  cat('Number of missing values =',x$nmiss)
  cat('\n')
  if (x$ncens>0) {
    cat("Convergence reached? =",x$converge)
    cat('\n')
    cat('Iterations =',x$iter,"/",x$MaxIter)
    cat('\n')
    cat('MC sample =',x$M)
    cat('\n')
    cat('Cut point =',x$pc)
    cat('\n')
  }
  cat("Processing time =",x$timediff,units(x$timediff))
  cat('\n','\n')
}

plot.ARpCRM<- function(x, ...) {
  if (x$ncens == 0) stop("plot only defined for cases with censoring")
  par(mar=c(4, 4.5, 1, 0.5))
  op <- suppressWarnings(par(mfrow=c(ifelse(x$plot$npar%%3==0,x$plot$npar%/%3,
                                            (x$plot$npar%/%3)+1),3)))
  for(i in 1:x$plot$npar)
  {
    suppressWarnings(plot.ts(x$plot$Theta[,i],xlab="Iteration",
                             ylab=x$plot$labels[[i]]))
    abline(v=x$plot$cpl,lty=2)
  }
  par(mfrow=c(1,1))
}
