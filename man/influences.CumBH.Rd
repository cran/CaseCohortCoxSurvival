\name{influences.CumBH}
\alias{influences.CumBH}
\title{influences.CumBH}
\description{
Computes the influences on the log-relative hazard, baseline hazards at each
unique event time, and on the cumulative baseline hazard in a given time interval
[Tau1, Tau2]. Can take calibration of the design weights into account.
}

\usage{
influences.CumBH(mod, Tau1 = NULL, Tau2 = NULL, A=NULL, calibrated = NULL)
}

\arguments{
\item{mod}{a cox model object, result of function coxph.}
\item{Tau1}{left bound of the time interval considered for the cumulative
  baseline hazard and pure risk. Default is the first event time.}
\item{Tau2}{right bound of the time interval considered for the cumulative
  baseline hazard and pure risk. Default is the last event time.}
\item{A}{\eqn{n \times q} matrix with the values of the auxiliary variables used for
  the calibration of the weights in the whole cohort. Needs to be provided if
 \code{calibrated = TRUE}.}
\item{calibrated}{are calibrated weights used for the estimation of the
  parameters? If \code{calibrated = TRUE}, the argument below needs to be
  provided. Default is \code{FALSE}.}
}

\value{

\code{infl.beta}: matrix with the overall influences on the log-relative hazard
estimates.

\code{infl.lambda0.t}: matrix with the overall influences on the baseline hazards
estimates at each unique event time.

\code{infl.Lambda0.Tau1Tau2.hat}: vector with the overall influences on the
cumulative baseline hazard estimate in [Tau1, Tau2].

\code{infl2.beta}: matrix with the phase-two influences on the log-relative
hazard estimates. Returned if \code{calibrated = TRUE}.

\code{infl2.lambda0.t}: matrix with the phase-two influences on the baseline
hazards estimates at each unique event time. Returned if \code{calibrated = TRUE}.

\code{infl2.Lambda0.Tau1Tau2.hat}: vector with the phase-two influences on the
cumulative baseline hazard estimate in [Tau1, Tau2]. Returned if
\code{calibrated = TRUE}.

\code{beta.hat}: vector of length \eqn{p} with log-relative hazard estimates.

\code{lambda0.t.hat}: vector with baseline hazards estimates at each unique event
time.

\code{Lambda0.Tau1Tau2.hat}: cumulative baseline hazard estimate in [Tau1, Tau2].
}

\details{

\code{influences.CumBH} works for estimation from a case-cohort with design weights
or calibrated weights (case-cohort consisting of the subcohort and cases not in
the subcohort, i.e., case-cohort obtained from two phases of sampling).

If covariate information is missing for certain individuals in the phase-two data
(i.e., case-cohort obtained from three phases of sampling), use \code{\link{influences.CumBH.missingdata}}.

\code{influences.CumBH} uses the influence formulas provided in Etievant and Gail
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
\code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.PR}},
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},
\code{\link{influences.CumBH.missingdata}},\cr
\code{\link{influences.PR.missingdata}}, \code{\link{robustvariance}} and  \code{\link{variance}}.
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

  # print the influences on the cumulative baseline hazard estimate
  # est$infl.Lambda0.Tau1Tau2

}}
