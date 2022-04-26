\name{rARCens}
\alias{rARCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Censored Autoregressive Data

}
\description{
It simulates a censored response variable with autoregressive errors of order p, with a established censoring rate.
}
\usage{
rARCens(n,beta,pit,sig2,x,cens,pcens)
}

\arguments{
 \item{n}{
Length of the desired time serie.
}
\item{beta}{
Vector of theoretical regression parameters of length \code{l}.
}
\item{pit}{
Vector of theoretical partial autocorrelations. Each element must be in (-1,1).
}
\item{sig2}{
Theoretical variance of the error.
}
\item{x}{
Matrix of covariates of dimension \code{n x l} (in models that include an intercept \code{x} should contain a column of ones).
}
\item{cens}{
"left" for left censoring, "right" for right censoring.
}
\item{pcens}{
Desired censoring rate.
}
}
\value{
\item{data}{Generated response (y) and censoring indicator (cc).}
\item{param}{Theoretical parameters (beta, sig2, phi).}
}

\author{Fernanda L. Schumacher <\email{fernandalschumacher@gmail.com}>,
Victor H. Lachos
<\email{hlachos@ime.unicamp.br}> and Christian E. Galarza <\email{cgalarza88@gmail.com}>

Maintainer: Fernanda L. Schumacher <\email{fernandalschumacher@gmail.com}>
}

\examples{
#generating a sample
dat = rARCens(n=100,beta = c(1,-1),pit = c(.4,-.2),
     sig2=.5,x=cbind(1,runif(100)),cens='left',pcens=.05)
#
plot.ts(dat$data$y)
table(dat$data$cc)

dat$param
#[1]  1.00 -1.00  0.50  0.48 -0.20
}


\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}

\seealso{\code{\link[ARCensReg]{ARCensReg}}}
