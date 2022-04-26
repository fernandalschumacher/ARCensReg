\name{phosphorus}
\alias{phosphorus}
\docType{data}
\title{
Phosphorus Concentration data
}
\description{
The phosphorus concentration (P) data of West Fork Cedar River at Finchford, Iowa, USA, collected under the ambient water quality program conducted by the Iowa Department of Natural Resources (Iowa DNR) were observed monthly from 10/1998 to 10/2013 (n=181) and the dataset was first available in the R package carx.
The phosphorus concentration measurement was subject to a limit of detection (lcl), thereby the P data are left censored.

The water discharge dataset were obtained from the website of U.S. Geological Survey (site number 05458900), and it is measured in cubic feet per second.
}

\usage{data(phosphorus)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{lP}}{logarithm of the phosphorus concentration.}
    \item{\code{cc}}{left censoring indicator (1 if the observation is left censored and 0 otherwise).}
    \item{\code{lQ}}{logarithm of the water discharge.}
    \item{\code{lcl}}{lower censoring limit.}
    \item{\code{time}}{Year-Month.}

  }
}
\source{
\url{http://waterdata.usgs.gov/ia/nwis/monthly/}

\url{https://cran.r-project.org/web/packages/carx/carx.pdf}
}

\seealso{\code{\link{ARCensReg}}}

\examples{
\dontrun{
data(phosphorus)

plot.ts(phosphorus$lP)
lines(phosphorus$lcl,col=2,lty=2)

#Proportion of censoring
prop.table(table(phosphorus$cc))

#A censored regression model
x = cbind(1,phosphorus$lQ)
miss =  which(is.na(phosphorus$lP))
AR_reg = ARCensReg(phosphorus$cc,phosphorus$lP,x,cens='left',miss=miss,p=1,tol=.001)
}
}
\keyword{datasets}
