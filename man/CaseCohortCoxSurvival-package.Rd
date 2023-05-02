
\name{CaseCohortCoxSurvival-package}
\alias{CaseCohortCoxSurvival}
\docType{package}
\title{
  Case-Cohort Cox Survival Inference
}
\description{
This package uses case-cohort data to estimate log-relative hazard, 
baseline hazards at each unique event time, 
cumulative baseline hazard in a given time interval and pure risk on the time interval and for a given
covariate profile, under the Cox model. For the corresponding variance estimation, 
it relies on influence functions and follows the
complete variance decomposition, to enable correct analysis of case-cohort data 
with and without stratification,
weight calibration or missing phase-two covariate data.
}
\details{
  The package provides functions implementing the methods described in Etievant and Gail (2023). More
  precisely, it includes
  \itemize{
    \item a main driver function, \code{\link{caseCohortCoxSurvival}}.
    \item one function, \code{\link{estimatePureRisk}}, to estimate pure risks and the corresponding variances with additional covariate profiles.
    \item three functions, \code{\link{estimation}}, \code{\link{estimation.CumBH}} and 
          \code{\link{estimation.PR}}, for parameters estimation.
    \item four functions, \code{\link{influences}}, \code{\link{influences.RH}}, 
          \code{\link{influences.CumBH}} and  \code{\link{influences.PR}}, for influence functions
          derivation when estimation is with design or calibrated weights and from a 
          case-cohort consisting of the subcohort and cases not in the subcohort 
          (i.e., case-cohort obtained from two phases of sampling).
    \item four functions, \code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, 
          \code{\link{influences.CumBH.missingdata}} and \code{\link{influences.PR.missingdata}},
          for influence functions derivation when estimation is with design
          weights and from a case-cohort when covariate information was missing for 
          certain individuals in the phase-two data 
         (i.e., case-cohort obtained from three phases of sampling).
    \item two functions, \code{\link{variance}} and \code{\link{variance.missingdata}},
          for variance estimation following complete variance
          decomposition (with design or calibrated weights and without missing 
          phase-two data, or with design weights and missing phase-two covariate data).
    \item one function, \code{\link{robustvariance}}, for robust variance estimation.
    \item one function, \code{\link{auxiliary.construction}}, to compute the auxiliary
               variables proposed by Breslow et al. (Stat. Biosci., 2009), Breslow and Lumley (IMS, 2013), and Shin et al. (Biometrics, 2020),.
    \item one function, \code{\link{calibration}}, for weight calibration.
    \item one function, \code{\link{estimation.weights.phase3}}, for estimating the phase-three weights.
  }
}
\author{Lola Etievant, Mitchell H. Gail}
\keyword{package}
\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

