\name{influences}
\alias{influences}
\title{influences}
\description{
Computes the influences on the log-relative hazard, baseline hazards at each 
unique event time, cumulative baseline hazard in a given time interval 
[Tau1, Tau2] and on the pure risk in [Tau1, Tau2] and for a given covariate 
profile x. Can take calibration of the design weights into account.
}

\usage{
influences(mod, Tau1 = NULL, Tau2 = NULL, x = NULL, calibrated = NULL, 
A = NULL)
}

\arguments{

\item{mod}{a cox model object, result of function coxph.}
\item{Tau1}{left bound of the time interval considered for the cumulative 
baseline hazard and pure risk. Default is the first event time.}
\item{Tau2}{right bound of the time interval considered for the cumulative 
baseline hazard and pure risk. Default is the last event time.}
\item{x}{vector of length \eqn{p}, specifying the covariate profile considered for the
pure risk. Default is (0,...,0).}
\item{calibrated}{are calibrated weights used for the estimation of the
parameters? If \code{calibrated = TRUE}, the argument below needs to be
provided. Default is \code{FALSE}.}
\item{A}{\eqn{n \times q} matrix with the values of the auxiliary variables used for
the calibration of the weights in the whole cohort. Needs to be provided if
\code{calibrated = TRUE}.}
}


\value{

\code{infl.beta}: matrix with the overall influences on the log-relative hazard estimates.

\code{infl.lambda0.t}: matrix with the overall influences on the baseline hazards estimates at each unique event time.

\code{infl.Lambda0.Tau1Tau2.hat}: vector with the overall influences on the cumulative baseline hazard estimate in [Tau1, Tau2].

\code{infl.Pi.x.Tau1Tau2.hat}: vector with the overall influences on the pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}.

\code{infl2.beta}: matrix with the phase-two influences on the log-relative hazard estimates. Returned if \code{calibrated = TRUE}.

\code{infl2.lambda0.t}: matrix with the phase-two influences on the baseline hazards estimates at each unique event time. Returned if \code{calibrated = TRUE}.

\code{infl2.Lambda0.Tau1Tau2.hat}: vector with the phase-two influences on the cumulative baseline hazard estimate in [Tau1, Tau2]. Returned if \code{calibrated = TRUE}.

\code{infl2.Pi.x.Tau1Tau2.hat}: vector with the phase-two influences on the pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}. Returned if \code{calibrated = TRUE}.

\code{beta.hat}: vector of length \eqn{p} with log-relative hazard estimates.

\code{lambda0.t.hat}: vector with baseline hazards estimates at each unique event time.

\code{Lambda0.Tau1Tau2.hat}: cumulative baseline hazard estimate in [Tau1, Tau2].

\code{Pi.x.Tau1Tau2.hat}: pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}.
}

\details{

\code{influences} works for estimation from a case-cohort with design weights
or calibrated weights (case-cohort consisting of the subcohort and cases not in 
the subcohort, i.e., case-cohort obtained from two phases of sampling).

If covariate information is missing for certain individuals in the phase-two data 
(i.e., case-cohort obtained from three phases of sampling), use \code{\link{influences.missingdata}}.

\code{influences} uses the influence formulas provided in Etievant and Gail 
(2023). More precisely, as in Section 3.2 if \code{calibrated = FALSE}, and as in
Section 4.3 if \code{calibrated = TRUE}. 

If \code{calibrated = FALSE}, the infuences are only provided for the individuals
in the case-cohort. If \code{calibrated = TRUE}, the influences are provided for
all the individuals in the cohort.
}

\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{
\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{estimation.PR}},
\code{\link{influences.RH}}, \code{\link{influences.CumBH}}, \code{\link{influences.PR}}, 
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, 
\code{\link{influences.CumBH.missingdata}}, 
\code{\link{influences.PR.missingdata}}, \code{\link{robustvariance}} and  \code{\link{variance}}.
}


\examples{
data(dataexample, package="CaseCohortCoxSurvival")

cohort      <- dataexample$cohort # a simulated cohort
casecohort  <- dataexample$casecohort # a simulated stratified case-cohort
A           <- dataexample$A # matrix with auxiliary variables values in the cohort

Tau1    <- 0 # given time interval for the pure risk
Tau2    <- 8
x       <- c(-1, 1, -0.6) # given covariate profile for the pure risk

# Estimation using the whole cohort 

mod.cohort <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = cohort, 
                    robust = TRUE)
estimation.cohort <- influences(mod.cohort, Tau1 = Tau1, Tau2 = Tau2, x = x)

# print the vector with log-relative hazard estimates
#estimation.cohort$beta.hat

# print the cumulative baseline hazard estimate
#estimation.cohort$Lambda0.Tau1Tau2.hat

# print the pure risk estimate
#estimation.cohort$Pi.x.Tau1Tau2.hat

# print the influences on the log-relative hazard estimates
#estimation.cohort$infl.beta

# print the influences on the cumulative baseline hazard estimate
#estimation.cohort$infl.Lambda0.Tau1Tau2

# print the influences on the pure risk estimate
#estimation.cohort$infl.Pi.x.Tau1Tau2

# Estimation using the stratified case-cohort with weights calibrated on A

mod.calib <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                   weight = weights.calib, id = id, robust = TRUE)
estimation.calib    <- influences(mod.calib, A = A, calibrated = TRUE, 
                                  Tau1 = Tau1, Tau2 = Tau2, x = x)

# print the influences on the log-relative hazard estimates
#estimation.calib$infl.beta

# print the influences on the cumulative baseline hazard estimate
#estimation.calib$infl.Lambda0.Tau1Tau2

# print the influences on the pure risk estimate
#estimation.calib$infl.Pi.x.Tau1Tau2

# print the phase-two influences on the log-relative hazard estimates
#estimation.calib$infl2.beta

# print the phase-two influences on the cumulative baseline hazard estimate
#estimation.calib$infl2.Lambda0.Tau1Tau2

# print the phase-two influences on the pure risk estimate
#estimation.calib$infl2.Pi.x.Tau1Tau2
}

