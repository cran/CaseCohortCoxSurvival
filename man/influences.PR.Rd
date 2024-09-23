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
(2024).

If \code{calibrated = FALSE}, the infuences are only provided for the individuals
in the case-cohort. If \code{calibrated = TRUE}, the influences are provided for
all the individuals in the cohort.
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{

\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{estimation.PR}},
\code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}},
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},
\code{\link{influences.CumBH.missingdata}}, \code{\link{influences.PR.missingdata}},
\code{\link{robustvariance}} and  \code{\link{variance}}.

}

\examples{
\donttest{
  data(dataexample.stratified, package="CaseCohortCoxSurvival")
  cohort <- dataexample.stratified$cohort
  casecohort <- cohort[which(cohort$status == 1 |
                       cohort$subcohort == 1),] # the stratified case-cohort
  casecohort$weights <- casecohort$strata.n / casecohort$strata.m
  casecohort$weights[which(casecohort$status == 1)] <- 1

  Tau1 <- 0
  Tau2 <- 8
  x <- c(-1, 1, -0.6) # given covariate profile for the pure risk

  # Estimation using the stratified case cohort with design weights
  mod <- coxph(Surv(event.time, status) ~ X1 + X2 + X3, data = casecohort,
               weight = weights, id = id, robust = TRUE)
  est <- influences(mod, Tau1 = Tau1, Tau2 = Tau2, x = x)

  # print the influences on the pure risk estimate
  # est$infl.Pi.x.Tau1Tau2
}}

