
\name{estimation.CumBH}
\alias{estimation.CumBH}
\title{estimation.CumBH}
\description{
Estimates the log-relative hazard, baseline hazards at each unique event time
and cumulative baseline hazard in a given time interval [Tau1, Tau2].
}
\usage{
estimation.CumBH(mod, Tau1 = NULL, Tau2 = NULL, missing.data = FALSE,
riskmat.phase2 = NULL, dNt.phase2 = NULL, status.phase2 = NULL)
}
\arguments{
  \item{mod}{a Cox model object, result of function coxph.}
  \item{Tau1}{left bound of the time interval considered for the cumulative
    baseline hazard. Default is the first event time.}
  \item{Tau2}{right bound of the time interval considered for the cumulative
    baseline hazard. Default is the last event time.}
  \item{missing.data}{was data on the \eqn{p} covariates missing for certain individuals
     in the phase-two data (i.e., was a third phase of sampling performed)? If
     \code{missing.data = TRUE}, the arguments below need to be provided.
      Default is FALSE.}
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

}
\details{
  \code{estimation.CumBH} returns the log-relative hazard estimates provided by
  \code{mod}, and estimates the baseline hazard point mass at any event time
non-parametrically.

  \code{estimation.CumBH} works for estimation from a case-cohort with design weights
or calibrated weights, when the case-cohort consists of the subcohort and cases
not in the subcohort (i.e., case-cohort obtained from two phases of sampling),
as well as with design weights when covariate data was missing for certain
individuals in the phase-two data (i.e., case-cohort obtained from three phases
of sampling).
}

\seealso{
 \code{\link{estimation}}, \code{\link{estimation.PR}}, \code{\link{influences}}, \code{\link{influences.RH}},
 \code{\link{influences.CumBH}}, \code{\link{influences.PR}},
 \code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, \code{\link{influences.CumBH.missingdata}},
and \code{\link{influences.PR.missingdata}}
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

  est.true <- estimation(mod.true, Tau1 = Tau1, Tau2 = Tau2, x = x,
                         missing.data = TRUE,
                         riskmat.phase2 = riskmat.phase2,
                         dNt.phase2 = dNt.phase2)

  est.true <- estimation.CumBH(mod.true, Tau1 = Tau1, Tau2 = Tau2,
                                            missing.data = TRUE,
                                            riskmat.phase2 = riskmat.phase2,
                                            dNt.phase2 = dNt.phase2)

  # print the cumulative baseline hazard estimate
  est.true$Lambda0.Tau1Tau2.hat
}}

\references{
Breslow, N. (1974). Covariance Analysis of Censored Survival Data. Biometrics, 30, 89-99.

Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}
