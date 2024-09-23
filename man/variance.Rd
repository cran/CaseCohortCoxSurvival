\name{variance}
\alias{variance}
\title{variance}
\description{
Computes the variance estimate that follows the complete
variance decomposition, for a parameter such as log-relative hazard, cumulative
baseline hazard or covariate specific pure-risk.
}

\usage{
variance(n, casecohort, weights = NULL, infl, calibrated = NULL,
infl2 = NULL, cohort = NULL, stratified = NULL,
variance.phase2 = NULL)
}

\arguments{

\item{n}{number of individuals in the whole cohort.}

\item{casecohort}{If \code{stratified = TRUE}, data frame with \code{status}
(case status), \code{weights} (design, if they are not
provided in the argument below), \code{W} (the \eqn{J} strata), \code{strata.m}
(vector of length \eqn{J} with the numbers of sampled individuals in the strata) and
\code{strata.n} (vector of length \eqn{J} with the strata sizes in the cohort), for
each individual in the stratified case-cohort data. If \code{stratified = FALSE},
data frame with \code{weights} (design, if they are not
provided in the argument below), \code{m} (number of sampled individuals) and
\code{n} (cohort size), for each individual in the unstratified case-cohort data.}

\item{weights}{vector with design weights for the individuals
in the case-cohort data.}

\item{infl}{matrix with the overall influences on the parameter.}

\item{calibrated}{are calibrated weights used for the estimation of the
parameters? If \code{calibrated = TRUE}, the arguments below need to be
provided. Default is \code{FALSE}.}

\item{infl2}{matrix with the phase-two influences on the parameter. Needs to be
provided if \code{calibrated = TRUE}.}

\item{cohort}{If \code{stratified = TRUE}, data frame with \code{status} (case
status) and \code{subcohort} (subcohort sampling indicators) for each individual
in the stratified case-cohort data. If \code{stratified = FALSE}, data frame with
\code{status} (case status) and \code{unstrat.subcohort} (subcohort unstratified
sampling indicators) for each individual in the unstratified case-cohort data.
Needs to be provided if \code{calibrated = TRUE}.}

\item{stratified}{was the sampling of the case-cohort stratified on \code{W}?
Default is \code{FALSE}.}

\item{variance.phase2}{should the phase-two variance component also be returned?
Default is \code{FALSE}.}
}

\value{

\code{variance}: variance estimate.

\code{variance.phase2}: phase-two variance component.
}

\details{

\code{variance} works for estimation from a case-cohort with design weights
or calibrated weights (case-cohort consisting of the subcohort and cases not in
the subcohort, i.e., case-cohort obtained from two phases of sampling).

If covariate information is missing for certain individuals in the phase-two data
(i.e., case-cohort obtained from three phases of sampling), use \code{\link{variance.missingdata}}.

\code{variance} uses the variance formulas provided in Etievant and Gail
(2024).
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{
\code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}},
\code{\link{influences.PR}}, \code{\link{robustvariance}} and \code{\link{variance.missingdata}}.
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
  n <- nrow(cohort)

  # Estimation using the stratified case-cohort with design weights
  mod <- coxph(Surv(event.time, status) ~ X1 + X2 + X3, data = casecohort,
               weight = weights, id = id, robust = TRUE)

  # parameters and influences estimation
  est           <- influences(mod, Tau1 = Tau1, Tau2 = Tau2, x = x)
  beta.hat      <- est$beta.hat
  Lambda0.hat   <- est$Lambda0.Tau1Tau2.hat
  Pi.x.hat      <- est$Pi.x.Tau1Tau2.hat
  infl.beta     <- est$infl.beta
  infl.Lambda0  <- est$infl.Lambda0.Tau1Tau2
  infl.Pi.x     <- est$infl.Pi.x.Tau1Tau2

  # variance estimate for the log-relative hazard estimate
  variance(n = n, casecohort = casecohort, infl = infl.beta, stratified = TRUE)

  # variance estimate for the cumulative baseline hazard estimate
  variance(n = n, casecohort = casecohort, infl = infl.Lambda0,
           stratified = TRUE)

  # variance estimate for the pure risk estimate
  variance(n = n, casecohort = casecohort, infl = infl.Pi.x, stratified = TRUE)
}
}

