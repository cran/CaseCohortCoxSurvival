\name{influences.PR}
\alias{influences.PR}
\title{influences.PR}
\description{
Computes the influences on the pure risk in the time interval [Tau1, Tau2] and 
for a given covariate profile x, from that on the log-relative hazard and 
cumulative baseline hazard. Can take calibration of the design weights into 
account.
}

\usage{
influences.PR(beta, Lambda0.Tau1Tau2, x = NULL, infl.beta, 
infl.Lambda0.Tau1Tau2, calibrated = NULL, infl2.beta = NULL, 
infl2.Lambda0.Tau1Tau2 = NULL)
}

\arguments{
\item{beta}{vector of length \eqn{p} with log-relative hazard values.}
\item{Lambda0.Tau1Tau2}{cumulative baseline hazard in [Tau1, Tau2].}
\item{x}{vector of length \eqn{p}, specifying the covariate profile considered for the
pure risk. Default is (0,...,0).}
\item{infl.beta}{matrix with the overall influences on the log-relative hazard 
estimates.}
\item{infl.Lambda0.Tau1Tau2}{vector with the overall influences on the 
cumulative baseline hazard estimate in [Tau1, Tau2].}
\item{calibrated}{are calibrated weights used for the estimation of the
parameters? If \code{calibrated = TRUE}, the arguments below need to be
provided. Default is \code{FALSE}.}
\item{infl2.beta}{matrix with the phase-two influences on the log-relative 
hazard estimates. Needs to be provided if \code{missing.data = TRUE}.}
\item{infl2.Lambda0.Tau1Tau2}{vector with the phase-two influences on the 
cumulative baseline hazard estimate in [Tau1, Tau2]. Needs to be provided 
if \code{missing.data = TRUE}.}
}

\value{

\code{infl.Pi.x.Tau1Tau2.hat}: vector with the overall influences on the pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}.

\code{infl2.Pi.x.Tau1Tau2.hat}: vector with the phase-two influences on the pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}. 
   Returned if \code{calibrated = TRUE}.

\code{Pi.x.Tau1Tau2.hat}: pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}.

}

\details{

\code{influences.PR} works for estimation from a case-cohort with design weights
or calibrated weights (case-cohort consisting of the subcohort and cases not in 
the subcohort, i.e., case-cohort obtained from two phases of sampling).

If covariate information is missing for certain individuals in the phase-two data 
(i.e., case-cohort obtained from three phases of sampling), use \code{\link{influences.PR.missingdata}}.

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
\code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}}, 
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, 
\code{\link{influences.CumBH.missingdata}}, \code{\link{influences.PR.missingdata}}, 
\code{\link{robustvariance}} and  \code{\link{variance}}.

}

\examples{
data(dataexample, package="CaseCohortCoxSurvival")
cohort      <- dataexample$cohort # a simulated cohort
casecohort  <- dataexample$casecohort # a simulated stratified case-cohort
A           <- dataexample$A # matrix with auxiliary variables values in the cohort

Tau1    <- 0 # given time interval for the pure risk
Tau2    <- 8
x       <- c(-1, 1, -0.6) # given covariate profile for the pure risk
v       <- c(1, 1, 0.6) # over covariate profile

# Estimation using the whole cohort

mod.cohort <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = cohort, 
                    robust = TRUE)
estimation.cohort     <- influences(mod.cohort, Tau1 = Tau1, Tau2 = Tau2, x = x)
beta.hat.cohort       <- estimation.cohort$beta.hat
Lambda0.hat.cohort    <- estimation.cohort$Lambda0.Tau1Tau2.hat
infl.beta.cohort      <- estimation.cohort$infl.beta
infl.Lambda0.cohort   <- estimation.cohort$infl.Lambda0.Tau1Tau2

estimation.cohort2    <- influences.PR(beta = beta.hat.cohort, 
                                       Lambda0.Tau1Tau2 = Lambda0.hat.cohort, 
                                       x = v, infl.beta = infl.beta.cohort, 
                                       infl.Lambda0.Tau1Tau2 = infl.Lambda0.cohort)

# print the influences on the pure risk estimate for profile v
#estimation.cohort2$infl.Pi.x.Tau1Tau2

# Estimation using the stratified case-cohort with weights calibrated on A

mod.calib <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                   weight = weights.calib, id = id, robust = TRUE)
estimation.calib    <- influences(mod.calib, A = A, calibrated = TRUE, 
                                  Tau1 = Tau1, Tau2 = Tau2, x = x)
beta.hat.calib      <- estimation.calib$beta.hat
Lambda0.hat.calib   <- estimation.calib$Lambda0.Tau1Tau2.hat
infl.beta.calib     <- estimation.calib$infl.beta
infl.Lambda0.calib  <- estimation.calib$infl.Lambda0.Tau1Tau2
infl2.beta.calib    <- estimation.calib$infl2.beta
infl2.Lambda0.calib <- estimation.calib$infl2.Lambda0.Tau1Tau2

estimation.calib2   <- influences.PR(beta = beta.hat.calib,
                                     Lambda0.Tau1Tau2 = Lambda0.hat.calib, 
                                     x = v, infl.beta = infl.beta.calib,
                                     infl.Lambda0.Tau1Tau2 = infl.Lambda0.calib,
                                     calibrated = TRUE, 
                                     infl2.beta = infl2.beta.calib,
                                     infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.calib)

# print the influences on the pure risk estimate for profile v
#estimation.calib2$infl.Pi.x.Tau1Tau2

# print the phase-two influences on the pure risk estimate for profile v
#estimation.calib2$infl2.Pi.x.Tau1Tau2
}

