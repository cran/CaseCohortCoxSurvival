\name{influences.RH.missingdata}
\alias{influences.RH.missingdata}
\title{influences.RH.missingdata}

\description{
Computes the influences on the log-relative hazard, when covariate data is 
missing for certain individuals in the phase-two data.
}

\usage{
influences.RH.missingdata(mod, riskmat.phase2, dNt.phase2 = NULL, 
status.phase2 = NULL, estimated.weights = FALSE, B.phase2 = NULL)
}

\arguments{

\item{mod}{a cox model object, result of function coxph.}

\item{riskmat.phase2}{at risk matrix for the phase-two data at all of the cases 
event times, even those with missing covariate data.}

\item{dNt.phase2}{counting process matrix for failures in the phase-two data. 
Needs to be provided if \code{status.phase2 = NULL}. }

\item{status.phase2}{vector indicating the case status in the phase-two data. 
Needs to be provided if \code{dNt.phase2 = NULL}.}

\item{estimated.weights}{are the weights for the third phase of sampling (due to 
missingness) estimated? If \code{estimated.weights = TRUE}, the argument below 
needs to beprovided. Default is \code{FALSE}.}

\item{B.phase2}{matrix for the phase-two data, with phase-three sampling strata 
indicators. It should have as many columns as phase-three strata (\eqn{J^{(3)}}), 
with one 1 per row, to indicate the phase-three stratum position. Needs to be 
provided if \code{estimated.weights = TRUE}.}
}

\value{

\code{infl.beta}: matrix with the overall influences on the log-relative hazard estimates.

\code{infl2.beta}: matrix with the phase-two influences on the log-relative hazard estimates.

\code{infl3.beta}: matrix with the phase-three influences on the log-relative hazard estimates.

\code{beta.hat}: vector of length \eqn{p} with log-relative hazard estimates.
}

\details{

\code{influences.RH.missingdata} works for estimation from a case-cohort with design 
weights and when covariate data was missing for certain individuals in the 
phase-two data (i.e., case-cohort obtained from three phases of sampling and 
consisting of individuals in the phase-two data without missing covariate information).

If there are no missing covariates in the phase- two sample, use \code{influences.RH}
with either design weights or calibrated weights.

\code{influences.RH.missingdata} uses the influence formulas provided in Etievant 
and Gail (2023). More precisely, as in Section 5.4 if 
\code{estimated.weights = TRUE}, and as in Web Appendix H if 
\code{estimated.weights = FALSE}.
} 

\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{
\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{estimation.PR}}, 
\code{\link{influences.missingdata}}, \code{\link{influences.CumBH.missingdata}},
\code{\link{influences.PR.missingdata}}, \code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}}, 
\code{\link{influences.PR}}, \code{\link{robustvariance}} and \code{\link{variance}}.
}

\examples{
data(dataexample.missingdata, package="CaseCohortCoxSurvival")

cohort          <- dataexample.missingdata$cohort # a simulated cohort
casecohort      <- dataexample.missingdata$casecohort # a simulated stratified case-cohort
# phase-two data: dataexample.missingdata$casecohort.phase2 
riskmat.phase2  <- dataexample.missingdata$riskmat.phase2
dNt.phase2      <- dataexample.missingdata$dNt.phase2
B.phase2        <- dataexample.missingdata$B.phase2

# Estimation using the stratified case-cohort with true known design weights 

mod <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
             weight = weights.true, id = id, robust = TRUE)

estimation <- influences.RH.missingdata(mod = mod, 
                                           riskmat.phase2 = riskmat.phase2, 
                                           dNt.phase2 = dNt.phase2)

# print the influences on the log-relative hazard estimates
#estimation$infl.beta

# print the phase-two influences on the log-relative hazard estimates
#estimation$infl2.beta

# print the phase-three influences on the log-relative hazard estimates
#estimation$infl3.beta

# Estimation using the stratified case-cohort with estimated weights, and
# accounting for the estimation through the influences

mod.est <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                 weight = weights.est, id = id, robust = TRUE)

estimation.est  <- influences.RH.missingdata(mod.est, 
                                                riskmat.phase2 = riskmat.phase2, 
                                                dNt.phase2 = dNt.phase2, 
                                                estimated.weights = TRUE,
                                                B.phase2 = B.phase2)
}
