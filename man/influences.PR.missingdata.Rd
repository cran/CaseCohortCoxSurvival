\name{influences.PR.missingdata}
\alias{influences.PR.missingdata}
\title{influences.PR.missingdata}

\description{
Computes the influences on the pure risk in the time interval [Tau1, Tau2] and 
for a given covariate profile x, from that on the log-relative hazard and 
cumulative baseline hazard, when covariate data is missing for certain individuals 
in the phase-two data.
}

\usage{
  influences.PR.missingdata(beta, Lambda0.Tau1Tau2, x = NULL, infl2.beta, 
  infl2.Lambda0.Tau1Tau2, infl3.beta, infl3.Lambda0.Tau1Tau2)
}

\arguments{

\item{beta}{vector of length \eqn{p} with log-relative hazard values.}

\item{Lambda0.Tau1Tau2}{cumulative baseline hazard in [Tau1, Tau2].}

\item{x}{vector of length \eqn{p}, specifying the covariate profile considered for 
the pure risk. Default is (0,...,0).}

\item{infl2.beta}{matrix with the overall influences on the log-relative hazard 
estimates.}

\item{infl2.Lambda0.Tau1Tau2}{vector with the overall influences on the 
cumulative baseline hazard estimate in [Tau1, Tau2].}

\item{infl3.beta}{matrix with the phase-three influences on the log-relative 
hazard estimates. }

\item{infl3.Lambda0.Tau1Tau2}{vector with the phase-three influences on the 
cumulative baseline hazard estimate in [Tau1, Tau2]. }
}


\value{

\code{infl.Pi.x.Tau1Tau2.hat}: vector with the overall influences on the pure risk estimate 
  in [Tau1, Tau2] and for covariate profile \code{x}.

\code{infl2.Pi.x.Tau1Tau2.hat}: vector with the phase-two influences on the pure risk estimate 
  in [Tau1, Tau2] and for covariate profile \code{x}. 

\code{infl3.Pi.x.Tau1Tau2.hat}: vector with the phase-three influences on the pure risk estimate 
 in [Tau1, Tau2] and for covariate profile \code{x}. 

 \code{Pi.x.Tau1Tau2.hat}: pure risk estimate in [Tau1, Tau2] and for covariate profile \code{x}.
}

\details{

\code{influences.PR.missingdata} works for estimation from a case-cohort with design 
weights and when covariate data was missing for certain individuals in the 
phase-two data (i.e., case-cohort obtained from three phases of sampling).

If there are no missing covariates in the phase- two sample, use \code{\link{influences.PR}}
with either design weights or calibrated weights.

\code{influences.PR.missingdata} uses the influence formulas provided in Etievant 
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
\code{\link{influences.CumBH.missingdata}}, \code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}}, 
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

Tau1    <- 0 # given time interval for the pure risk
Tau2    <- 8
x       <- c(-1, 1, -0.6) # given covariate profile for the pure risk
v       <- c(1, 1, 0.6) # over covariate profile

# Estimation using the stratified case-cohort with true known design weights

mod <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
             weight = weights.true, id = id, robust = TRUE)

estimation <- influences.missingdata(mod = mod, riskmat.phase2 = riskmat.phase2, 
                                     dNt.phase2 = dNt.phase2, Tau1 = Tau1, 
                                     Tau2 = Tau2, x = x)
beta.hat      <- estimation$beta.hat
Lambda0.hat   <- estimation$Lambda0.Tau1Tau2.hat
infl.beta     <- estimation$infl.beta
infl.Lambda0  <- estimation$infl.Lambda0.Tau1Tau2
infl2.beta    <- estimation$infl2.beta
infl2.Lambda0 <- estimation$infl2.Lambda0.Tau1Tau2
infl3.beta    <- estimation$infl3.beta
infl3.Lambda0 <- estimation$infl3.Lambda0.Tau1Tau2

estimation2   <- influences.PR.missingdata(beta = beta.hat, 
                                           Lambda0.Tau1Tau2 = Lambda0.hat, 
                                           x = v, infl2.beta = infl2.beta, 
                                           infl2.Lambda0.Tau1Tau2 = infl2.Lambda0,
                                           infl3.beta = infl3.beta, 
                                           infl3.Lambda0.Tau1Tau2 = infl3.Lambda0)

# print the overall influences on the pure risk estimate for profile v
#estimation2$infl.Pi.x.Tau1Tau2

# print the phase-two influences on the pure risk estimate for profile v
#estimation2$infl2.Pi.x.Tau1Tau2

# print the phase-three influences on the pure risk estimate for profile v
#estimation2$infl3.Pi.x.Tau1Tau2

# Estimation using the stratified case-cohort with estimated weights, and
# accounting for the estimation through the influences

mod.est <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                 weight = weights.est, id = id, robust = TRUE)

estimation.est  <- influences.missingdata(mod.est, 
                                          riskmat.phase2 = riskmat.phase2, 
                                          dNt.phase2 = dNt.phase2, 
                                          estimated.weights = TRUE,
                                          B.phase2 = B.phase2, Tau1 = Tau1, 
                                          Tau2 = Tau2, x = x)
beta.hat.est      <- estimation.est$beta.hat
Lambda0.hat.est   <- estimation.est$Lambda0.Tau1Tau2.hat
infl.beta.est     <- estimation.est$infl.beta
infl.Lambda0.est  <- estimation.est$infl.Lambda0.Tau1Tau2
infl2.beta.est    <- estimation.est$infl2.beta
infl2.Lambda0.est <- estimation.est$infl2.Lambda0.Tau1Tau2
infl3.beta.est    <- estimation.est$infl3.beta
infl3.Lambda0.est <- estimation.est$infl3.Lambda0.Tau1Tau2
estimation.est2   <- influences.PR.missingdata(beta = beta.hat.est,
                                               Lambda0.Tau1Tau2 = Lambda0.hat.est, 
                                               x = v, infl2.beta = infl2.beta.est,
                                               infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.est,
                                               infl3.beta = infl3.beta.est,
                                               infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.est)


}
