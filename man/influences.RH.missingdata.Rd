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
and Gail (2024).
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{
\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{estimation.PR}},
\code{\link{influences.missingdata}}, \code{\link{influences.CumBH.missingdata}},
\code{\link{influences.PR.missingdata}}, \code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}},
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

  # Estimation using the stratified case cohort with true known design weights
  mod.true <- coxph(Surv(event.time, status) ~ X1 + X2 + X3, data = casecohort,
                    weight = weight.true, id = id, robust = TRUE)

  est.true <- influences.missingdata(mod = mod.true, riskmat.phase2 = riskmat.phase2,
                                     dNt.phase2 = dNt.phase2, Tau1 = Tau1,
                                     Tau2 = Tau2, x = x)

  # print the influences on the log-relative hazard estimates
  # est.true$infl.beta
  # print the phase-two influences on the log-relative hazard estimates
  # est.true$infl2.beta
  # print the phase-three influences on the log-relative hazard estimates
  # est.true$infl3.beta

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

  # print the influences on the log-relative hazard estimates
  # est.estimated$infl.beta
  # print the phase-two influences on the log-relative hazard estimates
  # est.estimated$infl2.beta
  # print the phase-three influences on the log-relative hazard estimates
  # est.estimated$infl3.beta
}
}
