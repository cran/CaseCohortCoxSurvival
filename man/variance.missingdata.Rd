\name{variance.missingdata}
\alias{variance.missingdata}
\title{variance.missingdata}
\description{
Computes the variance estimate that follows the complete variance decomposition, 
for a parameter such as log-relative hazard, cumulative baseline hazard or 
covariate specific pure-risk, when covariate information is missing for 
individuals in the phase-two sample. 
}

\usage{
variance.missingdata(n, casecohort, casecohort.phase2, weights, 
weights.phase2, weights.p2.phase2, infl2, infl3, stratified.p2 = NULL, 
estimated.weights = NULL)
}

\arguments{
\item{n}{number of individuals in the whole cohort.}

\item{casecohort}{If \code{stratified = TRUE}, data frame with \code{W} 
 (the \eqn{J} phase-two strata), \code{strata.m} (vector of length \eqn{J} with the numbers of 
sampled individuals in the strata in the second phase of sampling) and 
\code{strata.n} (vector of length \eqn{J} with the strata sizes in the cohort), for 
each individual in the stratified case cohort data. If \code{stratified = FALSE},
data frame with \code{m} (number of sampled individuals in the second phase of 
sampling) and \code{n} (cohort size), for each individual in the unstratified 
case cohort data.}

\item{casecohort.phase2}{If \code{stratified = TRUE}, data frame with \code{W} 
(the \eqn{J} phase-two strata), \code{strata.m} (vector of length \eqn{J} with the 
numbers of sampled individuals in the strata in the second phase of sampling), 
\code{strata.n} (vector of length \eqn{J} with the strata sizes in the cohort) and 
\code{phase3} (phase-three sampling indicator), for each individual in the 
phase-two sample. If \code{stratified = FALSE}, data frame with  \code{m} (number 
of sampled individuals in the second phase of sampling), \code{n} (cohort size)
and unstrat.phase3 (phase-three sampling indicator), for each individual in the 
phase-two sample.}

\item{weights}{vector with design weights for the individuals
in the case cohort data.}

\item{weights.phase2}{vector with design weights for the individuals in the 
phase-two sample.}

\item{weights.p2.phase2}{vector with phase-two design weights for the individuals
in the phase-two sample.}

\item{infl2}{matrix with the phase-two influences on the parameter.}

\item{infl3}{matrix with the phase-three influences on the parameter.}

\item{stratified.p2}{was the second phase of sampling stratified on \code{W}?
Default is \code{FALSE}.}

\item{estimated.weights}{were the phase-three weights estimated? Default is \code{FALSE}.}
}

\value{
 \code{variance}: variance estimate.
}

\details{

\code{variance.missingdata} works for estimation from a case cohort with design 
weights and when covariate information was missing for certain individuals in the 
phase-two data (i.e., case cohort obtained from three phases of sampling and 
consisting of individuals in the phase-two data without missing covariate 
information).

If there are no missing covariates in the phase- two sample, use \code{\link{variance}} 
with either design weights or calibrated weights.

\code{variance.missingdata} uses the variance formulas provided in Etievant and 
Gail (2023). More precisely, as in Section 5.4 if \code{estimated.weights = TRUE}, and as in
Web Appendix H.2 if \code{estimated.weights = FALSE}.
} 

\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{

\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, \code{\link{influences.CumBH.missingdata}}, 

\code{\link{influences.PR.missingdata}}, \code{\link{robustvariance}} and \code{\link{variance}}.
}

\examples{

data(dataexample.missingdata, package="CaseCohortCoxSurvival")

cohort            <- dataexample.missingdata$cohort # a simulated cohort
n                 <- nrow(cohort)
casecohort        <- dataexample.missingdata$casecohort # a simulated stratified case cohort
casecohort.phase2 <- dataexample.missingdata$casecohort.phase2 
riskmat.phase2    <- dataexample.missingdata$riskmat.phase2
dNt.phase2        <- dataexample.missingdata$dNt.phase2
B.phase2          <- dataexample.missingdata$B.phase2

Tau1    <- 0 # given time interval for the pure risk
Tau2    <- 8
x       <- c(-1, 1, -0.6) # given covariate profile for the pure risk

# Estimation using the stratified case cohort with true known design weights

mod <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
             weight = weights.true, id = id, robust = TRUE)

estimation <- influences.missingdata(mod = mod, riskmat.phase2 = riskmat.phase2, 
                                     dNt.phase2 = dNt.phase2, Tau1 = Tau1, 
                                     Tau2 = Tau2, x = x)
infl.beta     <- estimation$infl.beta
infl.Lambda0  <- estimation$infl.Lambda0.Tau1Tau2
infl.Pi.x    <- estimation$infl.Pi.x.Tau1Tau2
infl2.beta    <- estimation$infl2.beta
infl2.Lambda0 <- estimation$infl2.Lambda0.Tau1Tau2
infl2.Pi.x   <- estimation$infl2.Pi.x.Tau1Tau2
infl3.beta    <- estimation$infl3.beta
infl3.Lambda0 <- estimation$infl3.Lambda0.Tau1Tau2
infl3.Pi.x   <- estimation$infl3.Pi.x.Tau1Tau2

# variance estimate for the log-relative hazard
variance.missingdata(n = n, casecohort = casecohort, 
                                 casecohort.phase2 = casecohort.phase2, 
                                 weights = casecohort$weights.true, 
                                 weights.phase2 = casecohort.phase2$weights.true, 
                                 weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                                 infl2 = infl2.beta, infl3 = infl3.beta, 
                                 stratified.p2 = TRUE)

# variance estimate for the cumulative baseline hazard estimate
variance.missingdata(n = n, casecohort = casecohort, 
                     casecohort.phase2 = casecohort.phase2, 
                     weights = casecohort$weights.true, 
                     weights.phase2 = casecohort.phase2$weights.true,
                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                     infl2 = infl2.Lambda0, infl3 = infl3.Lambda0, 
                     stratified.p2 = TRUE)

# variance estimate for the pure risk estimate
variance.missingdata(n = n, casecohort = casecohort, 
                     casecohort.phase2 = casecohort.phase2, 
                     weights = casecohort$weights.true, 
                     weights.phase2 = casecohort.phase2$weights.true, 
                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                     infl2 = infl2.Pi.x, infl3 = infl3.Pi.x, 
                     stratified.p2 = TRUE)


# Estimation using the stratified case cohort with estimated weights, and
# accounting for the estimation through the influences 

mod.est <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                 weight = weights.est, id = id, robust = TRUE)

estimation.est  <- influences.missingdata(mod.est, 
                                          riskmat.phase2 = riskmat.phase2, 
                                          dNt.phase2 = dNt.phase2, 
                                          estimated.weights = TRUE,
                                          B.phase2 = B.phase2, Tau1 = Tau1, 
                                          Tau2 = Tau2, x = x)
infl.beta.est     <- estimation.est$infl.beta
infl.Lambda0.est  <- estimation.est$infl.Lambda0.Tau1Tau2
infl.Pi.x.est     <- estimation.est$infl.Pi.x.Tau1Tau2
infl2.beta.est    <- estimation.est$infl2.beta
infl2.Lambda0.est <- estimation.est$infl2.Lambda0.Tau1Tau2
infl2.Pi.x.est    <- estimation.est$infl2.Pi.x.Tau1Tau2
infl3.beta.est    <- estimation.est$infl3.beta
infl3.Lambda0.est <- estimation.est$infl3.Lambda0.Tau1Tau2
infl3.Pi.x.est    <- estimation.est$infl3.Pi.x.Tau1Tau2

# variance estimate for the log-relative hazard
variance.missingdata(n = n, casecohort = casecohort, 
                     casecohort.phase2 = casecohort.phase2, 
                     weights = casecohort$weights.est, 
                     weights.phase2 = casecohort.phase2$weights.est, 
                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                     infl2 = infl2.beta.est, infl3 = infl3.beta.est, 
                     stratified.p2 = TRUE, estimated.weights = TRUE)

# variance estimate for the cumulative baseline hazard estimate
variance.missingdata(n = n, casecohort = casecohort,
                     casecohort.phase2 = casecohort.phase2, 
                     weights = casecohort$weights.est, 
                     weights.phase2 = casecohort.phase2$weights.est, 
                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                     infl2 = infl2.Lambda0.est, infl3 = infl3.Lambda0.est, 
                     stratified.p2 = TRUE, estimated.weights = TRUE)

# variance estimate for the pure risk estimate
variance.missingdata(n = n, casecohort = casecohort, 
                     casecohort.phase2 = casecohort.phase2, 
                     weights = casecohort$weights.est, 
                     weights.phase2 = casecohort.phase2$weights.est, 
                     weights.p2.phase2 = casecohort.phase2$weights.p2.true,
                     infl2 = infl2.Pi.x.est, infl3 = infl3.Pi.x.est, 
                     stratified.p2 = TRUE, estimated.weights = TRUE)
  
}
