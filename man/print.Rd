\name{print}
\alias{print.ARpCRM} \alias{print.ARtpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print an ARpCRM or ARtpCRM object
  %%  ~~function to do ... ~~
}
\description{
  Print an ARpCRM or ARtpCRM object.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  \method{print}{ARpCRM}(x, \dots)

  \method{print}{ARtpCRM}(x, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{ARpCRM} or \code{ARtpCRM}, 
  representing a fitted censored autoregressive model of order \eqn{p}.}
  \item{\dots}{Additional print arguments.}
}

\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza, and Larissa A. Matos
}


\seealso{
  \code{\link{ARCensReg}}, \code{\link{ARtCensReg}}, \code{\link{summary}}, \code{\link{plot}}
}
\examples{
n = 50; x = rep(1, n)
dat = rARCens(n=n, beta=2, phi=.5, sig2=.3, x=x, cens='left', pcens=.1)

fit = ARCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                M=5, pc=.12, tol=0.001, show_se=FALSE)
fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
