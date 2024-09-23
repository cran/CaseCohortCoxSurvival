
\name{dataexample}
\alias{dataexample}

\title{
[Deprecated] Data for examples
}
\description{
[\code{dataexample} is deprecated and will be removed in the next version of the package].

Simulated cohort, case-cohort and set of auxiliary variables for examples.
The case-cohort is a stratified phase-two sample with no missing covariate data.
}
%\details{ }
\examples{
 data(dataexample, package="CaseCohortCoxSurvival")

 # Display some of the data
 dataexample$cohort[1:5, ]

 dataexample$A[1:5, ] # auxiliary variable values in the cohort
}
\keyword{data}
\seealso{\code{\link{dataexample.stratified}}, \code{\link{dataexample.unstratified}}}



