\name{CloudCeiling}
\alias{CloudCeiling}
\docType{data}
\title{
Cloud Ceiling Height
}
\description{
The cloud ceiling heights, collected by the National Center for Atmospheric Research (NCAR), were observed hourly in San Francisco during the month of March 1989, consisting of n=716 observations (Park et al., 2007).
}

\usage{data(CloudCeiling)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{y}}{logarithm of the cloud ceiling heights.}
      \item{\code{cc}}{right censoring indicator (1 if the observation is right censored and 0 otherwise).}
  }
}
\source{
Park, J. W., Genton, M. G. & Ghosh, S. K. (2007). Censored time series analysis with autoregressive moving average models. Journal of Computational and Graphical Statistics, 18(4), 797-817.
}

\seealso{\code{\link{ARCensReg}}}

\examples{
\dontrun{
data(CloudCeiling)

plot.ts(CloudCeiling$y)

#Proportion of censoring
prop.table(table(CloudCeiling$cc))

#A censored regression model
##this may take a long time due to the number of censored observations.
##For other examples see help(ARCensReg).

x = as.matrix(rep(1,length(CloudCeiling$y)))
miss =  which(is.na(CloudCeiling$y))
AR_reg = ARCensReg(CloudCeiling$cc,CloudCeiling$y,x,cens='right',miss=miss,p=1,tol=.001)
}
}
\keyword{datasets}
