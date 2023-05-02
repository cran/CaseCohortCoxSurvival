\name{influences.CumBH.missingdata}
\alias{influences.CumBH.missingdata}
\title{influences.CumBH.missingdata}

\description{
Computes the influences on the log-relative hazard, baseline hazards at each 
unique event time, and on the cumulative baseline hazard in a given time 
interval [Tau1, Tau2], when covariate data is missing for certain individuals 
in the phase-two data.
}

\usage{
  influences.CumBH.missingdata(mod, riskmat.phase2, dNt.phase2 = NULL, 
status.phase2 = NULL, Tau1 = NULL, Tau2 = NULL, estimated.weights = FALSE, 
B.phase2 = NULL)
}

\arguments{

\item{mod}{a cox model object, result of function coxph.}

\item{riskmat.phase2}{at risk matrix for the phase-two data at all of the cases 
event times, even those with missing covariate data.}

\item{dNt.phase2}{counting process matrix for failures in the phase-two data. 
Needs to be provided if \code{status.phase2 = NULL}.} 

\item{status.phase2}{vector indicating the case status in the phase-two data. 
Needs to be provided if \code{dNt.phase2 = NULL}.}

\item{Tau1}{left bound of the time interval considered for the cumulative 
baseline hazard and pure risk. Default is the first event time.}

\item{Tau2}{right bound of the time interval considered for the cumulative 
baseline hazard and pure risk. Default is the last event time.}

\item{estimated.weights}{are the weights for the third phase of sampling (due to 
missingness) estimated? If \code{estimated.weights = TRUE}, the argument below 
needs to be provided. Default is \code{FALSE}.}

\item{B.phase2}{matrix for the phase-two data, with phase-three sampling strata 
indicators. It should have as many columns as phase-three strata (\eqn{J^{(3)}}), 
with one 1 per row, to indicate the phase-three stratum position. Needs to be 
provided if \code{estimated.weights = TRUE}.}
}

\value{

\code{infl.beta}: matrix with the overall influences on the log-relative hazard estimates.

\code{infl.lambda0.t}: matrix with the overall influences on the baseline hazards estimates at each unique event time.

\code{infl.Lambda0.Tau1Tau2.hat}: vector with the overall influences on the cumulative baseline hazard estimate in [Tau1, Tau2].

\code{infl2.beta}: matrix with the phase-two influences on the log-relative hazard estimates.

\code{infl2.lambda0.t}: matrix with the phase-two influences on the baseline hazards estimates at each unique event time.

\code{infl2.Lambda0.Tau1Tau2.hat}: vector with the phase-two influences on the cumulative baseline hazard estimate in [Tau1, Tau2]. 

\code{infl3.beta}: matrix with the phase-three influences on the log-relative hazard estimates.

\code{infl3.lambda0.t}: matrix with the phase-three influences on the baseline hazards estimates at each unique event time.

\code{infl3.Lambda0.Tau1Tau2.hat}: vector with the phase-three influences on the cumulative baseline hazard estimate in [Tau1, Tau2]. 

\code{beta.hat}: vector of length \eqn{p} with log-relative hazard estimates.

\code{lambda0.t.hat}: vector with baseline hazards estimates at each unique event time.

\code{Lambda0.Tau1Tau2.hat}: cumulative baseline hazard estimate in [Tau1, Tau2].
}

\details{

\code{influences.CumBH.missingdata} works for estimation from a case-cohort with design
weights and when covariate data was missing for certain individuals in the
phase-two data (i.e., case-cohort obtained from three phases of sampling).

If there are no missing covariates in the phase-two sample, use \code{\link{influences.CumBH}}
with either design weights or calibrated weights.

\code{influences.CumBH.missingdata} uses the influence formulas provided in Etievant
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
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},
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

estimation <- influences.CumBH.missingdata(mod = mod, 
                                           riskmat.phase2 = riskmat.phase2, 
                                           dNt.phase2 = dNt.phase2, Tau1 = 0, 
                                           Tau2 = 8)

# print the influences on the log-relative hazard estimates
#estimation$infl.beta

# print the influences on the cumulative baseline hazard estimate
#estimation$infl.Lambda0.Tau1Tau2

# print the phase-two influences on the log-relative hazard estimates
#estimation$infl2.beta

# print the phase-two influences on the cumulative baseline hazard estimate
#estimation$infl2.Lambda0.Tau1Tau2

# print the phase-three influences on the log-relative hazard estimates
#estimation$infl3.beta

# print the phase-three influences on the cumulative baseline hazard estimate
#estimation$infl3.Lambda0.Tau1Tau2

# Estimation using the stratified case-cohort with estimated weights, and
# accounting for the estimation through the influences 

mod.est <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                 weight = weights.est, id = id, robust = TRUE)

estimation.est  <- influences.CumBH.missingdata(mod.est, 
                                                riskmat.phase2 = riskmat.phase2, 
                                                dNt.phase2 = dNt.phase2, 
                                                estimated.weights = TRUE,
                                                B.phase2 = B.phase2, 
                                                Tau1 = 0, Tau2 = 8)
}
