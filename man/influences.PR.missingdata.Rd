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
and Gail (2024).
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{
\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{estimation.PR}},
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},
\code{\link{influences.CumBH.missingdata}}, \code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}},
\code{\link{influences.PR}}, \code{\link{robustvariance}} and \code{\link{variance}}.
}

\examples{
\donttest{
  data(dataexample.missingdata.stratified, package="CaseCohortCoxSurvival")

  cohort <- dataexample.missingdata.stratified$cohort
  phase2 <- cohort[which(cohort$phase2 == 1),] # the phase-two sample
  casecohort <- cohort[which(cohort$phase3 == 1),] # the stratified case-cohort

  B.phase2 <- cbind(1 * (phase2$W3 == 0), 1 * (phase2$W3 == 1))
  rownames(B.phase2)  <- cohort[cohort$phase2 == 1, "id"]
  B.phase3 <- cbind(1 * (casecohort$W3 == 0), 1 * (casecohort$W3 == 1))
  rownames(B.phase3)  <- cohort[cohort$phase3 == 1, "id"]
  total.B.phase2 <- colSums(B.phase2)
  J3 <- ncol(B.phase3)
  n <- nrow(cohort)

  # Quantities needed for estimation of the cumulative baseline hazard when
  # covariate data is missing
  mod.cohort <- coxph(Surv(event.time, status) ~ X2, data = cohort,
                      robust = TRUE) # X2 is available on all cohort members
  mod.cohort.detail <- coxph.detail(mod.cohort, riskmat = TRUE)

  riskmat.phase2 <- with(cohort, mod.cohort.detail$riskmat[phase2 == 1,])
  rownames(riskmat.phase2) <- cohort[cohort$phase2 == 1, "id"]
  observed.times.phase2 <- apply(riskmat.phase2, 1,
                                 function(v) {which.max(cumsum(v))})
  dNt.phase2 <- matrix(0, nrow(riskmat.phase2), ncol(riskmat.phase2))
  dNt.phase2[cbind(1:nrow(riskmat.phase2), observed.times.phase2)] <- 1
  dNt.phase2 <- sweep(dNt.phase2, 1, phase2$status, "*")
  colnames(dNt.phase2) <- colnames(riskmat.phase2)
  rownames(dNt.phase2) <- rownames(riskmat.phase2)

  Tau1 <- 0 # given time interval for the pure risk
  Tau2 <- 8
  x <- c(-1, 1, -0.6) # given covariate profile for the pure risk
  v <- c(1, -1, 0.6) # over covariate profile

  # Estimation using the stratified case cohort with true known design weights
  mod.true <- coxph(Surv(event.time, status) ~ X1 + X2 + X3, data = casecohort,
                    weight = weight.true, id = id, robust = TRUE)

  est.true <- influences.missingdata(mod = mod.true,
                                     riskmat.phase2 = riskmat.phase2,
                                     dNt.phase2 = dNt.phase2, Tau1 = Tau1,
                                     Tau2 = Tau2, x = x)

  beta.true <- est.true$beta.hat
  Lambda0.true <- est.true$Lambda0.Tau1Tau2.hat
  infl2.beta.true <- est.true$infl2.beta
  infl2.Lambda0.true <- est.true$infl2.Lambda0.Tau1Tau2
  infl3.beta.true <- est.true$infl3.beta
  infl3.Lambda0.true <- est.true$infl3.Lambda0.Tau1Tau2

  est.PR2.true <- influences.PR.missingdata(beta = beta.true,
                                            Lambda0.Tau1Tau2 = Lambda0.true,
                                            x = v,
                                            infl2.beta = infl2.beta.true,
                                            infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.true,
                                            infl3.beta = infl3.beta.true,
                                            infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.true)

  # print the influences on the pure risk estimate
  # est.PR2.true$infl.Pi.x.Tau1Tau2
  # print the phase-two influences on the pure risk estimate
  # est.PR2.true$infl2.Pi.x.Tau1Tau2
  # print the phase-three influences on the pure risk estimate
  # est.PR2.true$infl3.Pi.x.Tau1Tau2

  # Estimation using the stratified case cohort with estimated weights, and
  # accounting for the estimation through the influences
  mod.estimated <- coxph(Surv(event.time, status) ~ X1 + X2 + X3,
                         data = casecohort, weight = weight.est, id = id,
                         robust = TRUE)

  est.estimated  <- influences.missingdata(mod.estimated,
                                           riskmat.phase2 = riskmat.phase2,
                                           dNt.phase2 = dNt.phase2,
                                           estimated.weights = TRUE,
                                           B.phase2 = B.phase2, Tau1 = Tau1,
                                           Tau2 = Tau2, x = x)

  beta.estimated <- est.estimated$beta.hat
  Lambda0.estimated <- est.estimated$Lambda0.Tau1Tau2.hat
  infl2.beta.estimated <- est.estimated$infl2.beta
  infl2.Lambda0.estimated <- est.estimated$infl2.Lambda0.Tau1Tau2
  infl3.beta.estimated <- est.estimated$infl3.beta
  infl3.Lambda0.estimated <- est.estimated$infl3.Lambda0.Tau1Tau2

  est.PR2.estimated <- influences.PR.missingdata(beta = beta.estimated,
                                                 Lambda0.Tau1Tau2 = Lambda0.estimated,
                                                 x = v,
                                                 infl2.beta = infl2.beta.estimated,
                                                 infl2.Lambda0.Tau1Tau2 = infl2.Lambda0.estimated,
                                                 infl3.beta = infl3.beta.estimated,
                                                 infl3.Lambda0.Tau1Tau2 = infl3.Lambda0.estimated)

  # print the influences on the pure risk estimate
  # est.PR2.estimated$infl.Pi.x.Tau1Tau2
  # print the phase-two influences on the pure risk estimate
  # est.PR2.estimated$infl2.Pi.x.Tau1Tau2
  # print the phase-three influences on the pure risk estimate
  # est.PR2.estimated$infl3.Pi.x.Tau1Tau2
}
}
