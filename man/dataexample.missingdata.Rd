
\name{dataexample.missingdata}
\alias{dataexample.missingdata}

\title{
[Deprecated] Data for examples with missing data
}
\description{
[\code{dataexample.missingdata} is deprecated and will be removed in the next version of the package].

Simulated cohort and case-cohort for examples.
The case-cohort is a stratified phase-three sample, because of missing covariate
information in the stratified phase-two data.
}
%\details{ }
\examples{
 data(dataexample.missingdata, package="CaseCohortCoxSurvival")

 # Display some of the data
 dataexample.missingdata$cohort[1:5, ]
}
\keyword{data}
\seealso{\code{\link{dataexample.missingdata.stratified}}, \code{\link{dataexample.missingdata.unstratified}}}



