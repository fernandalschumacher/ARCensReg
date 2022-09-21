\name{summary}
\alias{summary.ARpCRM} \alias{summary.ARtpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of an ARpCRM or ARtpCRM object
  %%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "ARpCRM" or "ARtpCRM".
}

\usage{
  \method{summary}{ARpCRM}(object, \dots)
  
  \method{summary}{ARtpCRM}(object, \dots)
}

\arguments{
  \item{object}{An object inheriting from class \code{ARpCRM} or \code{ARtpCRM},  
  representing a fitted censored autoregressive model of order \eqn{p}.}
  \item{\dots}{Additional arguments.}
}

\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza, and Larissa A. Matos
}


\seealso{
  \code{\link{ARCensReg}}, \code{\link{ARtCensReg}}, \code{\link{print}}, \code{\link{plot}}
}
\examples{
n = 80; x = rep(1, n)
dat = rARCens(n=n, beta=2, phi=.6, sig2=.3, x=x, cens='right', pcens=.05, 
              innov='t', nu=4)

fit = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                 M=5, pc=.12, tol=0.001)
summary(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
