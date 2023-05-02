
\name{dataexample.missingdata}
\alias{dataexample.missingdata}

\title{
  Data for examples with missing data
}
\description{
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
%\seealso{\code{\link{caseCohortCoxSurvival}}}



