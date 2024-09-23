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
Gail (2024).
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{

\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, \code{\link{influences.CumBH.missingdata}},

\code{\link{influences.PR.missingdata}}, \code{\link{robustvariance}} and \code{\link{variance}}.
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

  # Estimation using the stratified case cohort with true known design weights

  mod.true <- coxph(Surv(event.time, status) ~ X1 + X2 + X3, data = casecohort,
                    weight = weight.true, id = id, robust = TRUE)

  est.true <- influences.missingdata(mod = mod.true,
                                     riskmat.phase2 = riskmat.phase2,
                                     dNt.phase2 = dNt.phase2, Tau1 = Tau1,
                                     Tau2 = Tau2, x = x)
  infl.beta.true <- est.true$infl.beta
  infl.Lambda0.true <- est.true$infl.Lambda0.Tau1Tau2
  infl.Pi.x.true <- est.true$infl.Pi.x.Tau1Tau2
  infl2.beta.true <- est.true$infl2.beta
  infl2.Lambda0.true <- est.true$infl2.Lambda0.Tau1Tau2
  infl2.Pi.x.true <- est.true$infl2.Pi.x.Tau1Tau2
  infl3.beta.true <- est.true$infl3.beta
  infl3.Lambda0.true <- est.true$infl3.Lambda0.Tau1Tau2
  infl3.Pi.x.true <- est.true$infl3.Pi.x.Tau1Tau2

  # variance estimate for the log-relative hazard estimate
  variance.missingdata(n = n, casecohort = casecohort,
                       casecohort.phase2 = phase2,
                       weights = casecohort$weight.true,
                       weights.phase2 = phase2$weight.true,
                       weights.p2.phase2 = phase2$weight.p2.true,
                       infl2 = infl2.beta.true, infl3 = infl3.beta.true,
                       stratified.p2 = TRUE)

  # variance estimate for the cumulative baseline hazard estimate
  variance.missingdata(n = n, casecohort = casecohort,
                       casecohort.phase2 = phase2,
                       weights = casecohort$weight.true,
                       weights.phase2 = phase2$weight.true,
                       weights.p2.phase2 = phase2$weight.p2.true,
                       infl2 = infl2.Lambda0.true, infl3 = infl3.Lambda0.true,
                       stratified.p2 = TRUE)

  # variance estimate for the pure risk estimate
  variance.missingdata(n = n, casecohort = casecohort,
                       casecohort.phase2 = phase2,
                       weights = casecohort$weight.true,
                       weights.phase2 = phase2$weight.true,
                       weights.p2.phase2 = phase2$weight.p2.true,
                       infl2 = infl2.Pi.x.true, infl3 = infl3.Pi.x.true,
                       stratified.p2 = TRUE)


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

  infl.beta.estimated <- est.estimated$infl.beta
  infl.Lambda0.estimated <- est.estimated$infl.Lambda0.Tau1Tau2
  infl.Pi.x.estimated <- est.estimated$infl.Pi.x.Tau1Tau2
  infl2.beta.estimated <- est.estimated$infl2.beta
  infl2.Lambda0.estimated <- est.estimated$infl2.Lambda0.Tau1Tau2
  infl2.Pi.x.estimated <- est.estimated$infl2.Pi.x.Tau1Tau2
  infl3.beta.estimated <- est.estimated$infl3.beta
  infl3.Lambda0.estimated <- est.estimated$infl3.Lambda0.Tau1Tau2
  infl3.Pi.x.estimated <- est.estimated$infl3.Pi.x.Tau1Tau2

  # variance estimate for the log-relative hazard
  variance.missingdata(n = n, casecohort = casecohort,
                       casecohort.phase2 = phase2,
                       weights = casecohort$weight.est,
                       weights.phase2 = phase2$weight.est,
                       weights.p2.phase2 = phase2$weight.p2.true,
                       infl2 = infl2.beta.estimated,
                       infl3 = infl3.beta.estimated,
                       stratified.p2 = TRUE, estimated.weights = TRUE)

  # variance estimate for the cumulative baseline hazard estimate
  variance.missingdata(n = n, casecohort = casecohort,
                       casecohort.phase2 = phase2,
                       weights = casecohort$weight.est,
                       weights.phase2 = phase2$weight.est,
                       weights.p2.phase2 = phase2$weight.p2.true,
                       infl2 = infl2.Lambda0.estimated,
                       infl3 = infl3.Lambda0.estimated,
                       stratified.p2 = TRUE, estimated.weights = TRUE)

  # variance estimate for the pure risk estimate
  variance.missingdata(n = n, casecohort = casecohort,
                       casecohort.phase2 = phase2,
                       weights = casecohort$weight.est,
                       weights.phase2 = phase2$weight.est,
                       weights.p2.phase2 = phase2$weight.p2.true,
                       infl2 = infl2.Pi.x.estimated,
                       infl3 = infl3.Pi.x.estimated,
                       stratified.p2 = TRUE, estimated.weights = TRUE)
}
}
