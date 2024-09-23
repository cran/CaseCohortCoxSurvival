
\name{auxiliary.construction}
\alias{auxiliary.construction}
\title{auxiliary.construction}
\description{
  Creates the auxiliary variables proposed by Breslow et al. (Stat. Biosci., 2009),
   Breslow and Lumley (IMS, 2013), and proposed by Shin et al. (Biometrics, 2020).
}
\usage{
auxiliary.construction(mod, Tau1 = NULL, Tau2 = NULL, method = "Breslow",
time.on.study = NULL, casecohort = NULL)
}
\arguments{
  \item{mod}{A cox model object, result of function coxph run on the cohort data with imputed covariate values.}
  \item{Tau1}{Left bound of the time interval considered for the cumulative baseline hazard.
           Default is the first event time.}
  \item{Tau2}{Right bound of the time interval considered for the cumulative baseline hazard.
             Default is the last event time.}
  \item{method}{"Breslow", "Breslow2013" or "Shin" to specify the algorithm to construct
                the auxiliary variables. The default is "Breslow".}
  \item{time.on.study}{Total folow-up time in \code{[Tau1, Tau2]}.
                       Required for \code{method} = "Shin".}
  \item{casecohort}{Data frame containing the casecohort data.
                   It must include columns "weights" containing
                   the design weights and "id" as an id variable.
                   Required for \code{method} = "Shin".}
}
\value{
  \code{A.RH.Breslow}: matrix with the influences on the log-relative hazard,
                    estimated from the cohort with imputed phase-two covariate
                    values for \code{method} = "Breslow" and
                    \code{method} = "Breslow2013".


  \code{A.CumBH.Breslow}: matrix with the influences on the cumulative baseline
                          hazard in \code{[Tau1, Tau2]},
                          estimated from the cohort with imputed phase-two
                          covariate values for \code{method} = "Breslow2013".


  \code{A.RH.Shin}: matrix with the influences on the log-relative hazard,
                    estimated from the cohort with imputed phase-two covariate
                    values for \code{method} = "Shin".


  \code{A.PR.Shin}: matrix with the products of total follow-up times in
                    \code{[Tau1, Tau2]} and estimated relative hazards,
                    estimated from the cohort with imputed phase-two
                    covariate values for \code{method} = "Shin".


}
\details{
  Construction of the auxiliary variables can follow Breslow et al. (2009),
Breslow and Lumley (2013), or Shin et al. (2020) (method).
It relies on predictions of the phase-two covariates for all members of the cohort.
The auxiliary variables are given by (i) the influences for the log-relative hazard parameters
estimated from the Cox model with imputed cohort data;
(ii) the influences for the cumulative baseline parameter estimated from the Cox model with
imputed cohort data;
(iii) the products of total follow-up time (on the time interval for which pure risk is to be estimated)
 with the estimated relative hazard for the imputed cohort data,
where the log-relative hazard parameters are estimated from the Cox model with case-cohort data and
weights calibrated with (i).
When method = Breslow, calibration of the design weights is against (i),
as proposed by Breslow et al. (2009) to improve efficiency of case-cohort estimates of relative hazard.
When method = Breslow2013, calibration of the design weights is against (i) and (ii),
as proposed by Breslow and Lumley (2013) to also improve efficiency of case-cohort estimates
of cumulative baseline hazard.
 When method = Shin, calibration is against (i) and (iii), as proposed by Shin et al. (2020)
to improve efficiency of relative hazard and pure risk estimates under the nested case-control design.
See Etievant and Gail (2024).


Following Etievant and Gail (2024), in function \code{\link{caseCohortCoxSurvival}}
 we only provide calibration
of the design weight as proposed by Breslow et al. (2009) or Shin et al. (2020).
}

\seealso{ \code{\link{calibration}}, \code{\link{influences}}, \code{\link{influences.RH}},
          \code{\link{influences.CumBH}} and \code{\link{influences.PR}}. }

\examples{
\donttest{
  data(dataexample.stratified, package="CaseCohortCoxSurvival")
  cohort <- dataexample.stratified$cohort
  Tau1   <- 0
  Tau2   <- 8

  # Running the coxph model on the imputed cohort data
  mod.imputedcohort <- coxph(Surv(event.time, status) ~ X1.pred + X2 + X3.pred,
                             data = cohort, robust = TRUE)

  # method = Breslow
  ret <- auxiliary.construction(mod.imputedcohort)
  # print auxiliary variables based on the log-relative hazard influences
  ret$A.RH.Breslow[1:5,]

  # Example for method = Shin, variables names must match with fitted model
  casecohort <- cohort[which(cohort$status == 1 |
                       cohort$subcohort == 1),] # the stratified case-cohort
  casecohort$weights <- casecohort$strata.n / casecohort$strata.m
  casecohort$weights[which(casecohort$status == 1)] <- 1
  casecohort[, "X1.pred"] <- casecohort[, "X1"]
  casecohort[, "X3.pred"] <- casecohort[, "X3"]

  time.on.study <- pmax(pmin(Tau2, cohort$event.time) - Tau1, 0)
  ret <- auxiliary.construction(mod.imputedcohort, method = "Shin",
                                time.on.study = time.on.study, casecohort = casecohort)
  ret$A.PR.Shin[1:5]

}}

\references{

Breslow, N.E. and Lumley, T. (2013). Semiparametric models and two-phase samples: Applications to Cox
regression. From Probability to Statistics and Back: High-Dimensional Models and Processes, 9, 65-78.

Breslow, N.E., Lumley, T., Ballantyne, C.M., Chambless, L.E. and Kulich, M. (2009). Improved Horvitz-Thompson Estimation of Model Parameters from Two-phase Stratified Samples: Applications in Epidemiology.
Statistics in Biosciences, 1, 32-49.

Shin Y.E., Pfeiffer R.M., Graubard B.I., Gail M.H. (2020) Weight calibration to improve
the efficiency of pure risk estimates from case-control samples nested in a cohort.
Biometrics, 76, 1087-1097.

Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.

}
