\name{estimation}
\alias{estimation}
\title{estimation}
\description{
Estimates the log-relative hazard, baseline hazards at each unique event time, 
cumulative baseline hazard in a given time interval [Tau1, Tau2] and pure risk 
in [Tau1, Tau2] and for a given covariate profile x.
}

\usage{
estimation(mod, Tau1 = NULL, Tau2 = NULL, x = NULL, missing.data = NULL, 
riskmat.phase2 = NULL, dNt.phase2 = NULL, status.phase2 = NULL)
}

\arguments{
\item{mod}{a Cox model object, result of function \code{coxph}.}
\item{Tau1}{left bound of the time interval considered for the cumulative 
baseline hazard and pure risk. Default is the first event time.}
\item{Tau2}{right bound of the time interval considered for the cumulative 
baseline hazard and pure risk. Default is the last event time.}
\item{x}{vector of length \eqn{p}, specifying the covariate profile considered for the
pure risk. Default is (0,...,0).}
\item{missing.data}{was data on the \eqn{p} covariates missing for certain individuals 
in the phase-two data (i.e., was a third phase of sampling performed)? If 
\code{missing.data = TRUE}, the arguments below need to be provided. Default is \code{FALSE}.}
\item{riskmat.phase2}{at risk matrix for the phase-two data at all of the case 
event times, even those with missing covariate data. Needs to be provided if
\code{missing.data = TRUE}.}
\item{dNt.phase2}{counting process matrix for failures in the phase-two data. 
Needs to be provided if \code{missing.data = TRUE} and \code{status.phase2 = NULL}.}
\item{status.phase2}{vector indicating the case status in the phase-two data. 
Needs to be provided if \code{missing.data = TRUE} and \code{dNt.phase2 = NULL}.}
}

\value{
\code{beta.hat}: vector of length \eqn{p} with log-relative hazard estimates.

\code{lambda0.t.hat}: vector with baseline hazards estimates at each unique event time.

\code{Lambda0.Tau1Tau2.hat}: cumulative baseline hazard estimate in [Tau1, Tau2].

\code{Pi.x.Tau1Tau2.hat}: pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}.
}

\details{

  \code{estimation} returns the log-relative hazard estimates provided by 
  \code{mod}, and estimates the baseline hazard point mass at any event time 
non-parametrically.

  \code{estimation} works for estimation from a case-cohort with design weights
or calibrated weights, when the case-cohort consists of the subcohort and cases 
not in the subcohort (i.e., case-cohort obtained from two phases of sampling), 
as well as with design weights when covariate data was missing for certain 
individuals in the phase-two data (i.e., case-cohort obtained from three phases 
of sampling).
}

\references{
Breslow, N. (1974). Covariance Analysis of Censored Survival Data. Biometrics, 30, 89-99.

Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{
\code{\link{estimation.CumBH}}, \code{\link{estimation.PR}}, \code{\link{influences}}, \code{\link{influences.RH}},
 \code{\link{influences.CumBH}}, \code{\link{influences.PR}}, 
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, \code{\link{influences.CumBH.missingdata}},
and \code{\link{influences.PR.missingdata}}.
}

\examples{
data(dataexample.missingdata, package="CaseCohortCoxSurvival")
cohort      <- dataexample.missingdata$cohort # a simulated cohort
casecohort  <- dataexample.missingdata$casecohort # a simulated stratified case-cohort
riskmat.phase2  <- dataexample.missingdata$riskmat.phase2
dNt.phase2      <- dataexample.missingdata$dNt.phase2

Tau1    <- 0 # given time interval for the pure risk
Tau2    <- 8
x       <- c(-1, 1, -0.6) # given covariate profile for the pure risk

# Estimation using the whole cohort

mod.cohort <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = cohort, 
                    robust = TRUE)
estimation.cohort <- estimation(mod = mod.cohort, Tau1 = Tau1, Tau2 = Tau2, 
                                x = x)

# print the vector with log-relative hazard estimates
estimation.cohort$beta.hat

# print the cumulative baseline hazard estimate
estimation.cohort$Lambda0.Tau1Tau2.hat

# print the pure risk estimate
estimation.cohort$Pi.x.Tau1Tau2.hat

# Estimation using the stratified case-cohort with known design weights

mod <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
             weight = weights.true, id = id, robust = TRUE)
estimation.casecohort <- estimation(mod = mod, Tau1 = Tau1, Tau2 = Tau2, x = x,
                                    missing.data = TRUE, 
                                    riskmat.phase2 = riskmat.phase2, 
                                    dNt.phase2 = dNt.phase2)

# print the vector with log-relative hazard estimates
estimation.casecohort$beta.hat

# print the cumulative baseline hazard estimate
estimation.casecohort$Lambda0.Tau1Tau2.hat

# print the pure risk estimate
estimation.casecohort$Pi.x.Tau1Tau2.hat
}
