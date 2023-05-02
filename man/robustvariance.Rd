\name{robustvariance}
\alias{robustvariance}
\title{robustvariance}
\description{
Computes the robust variance estimate, i.e., the sum of the squared influence 
functions, for a parameter such as log-relative hazard, cumulative baseline 
hazard or covariate specific pure-risk.
}

\usage{
 robustvariance(infl) 
}

\arguments{
 \item{infl}{overall influences on a parameter such as log-relative hazard, 
cumulative baseline hazard or covariate specific pure-risk.}
}

\value{
 \code{robust.var}: robust variance estimate.
}

\details{
 \code{robustvariance} works for estimation with design or calibrated weights from
a case cohort obtained from two phases of sampling (i.e., case cohort consisting
of the subcohort and cases not in the subcohort), or when covariate information 
was missing for certain individuals in the phase-two data (i.e., case cohort 
obtained from three phases of sampling and consisting of individuals in the 
phase-two data without missing covariate information).
}

\references{
Barlow W. (1994). Robust Variance Estimation for the Case-Cohort Design. Biometrics, 50, 1064-1072.

Langholz B., Jiao J. (2007). Computational methods for case-cohort studies.
Computational Statistics & Data Analysis, 51, 3737-37.

Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{
\code{\link{influences.RH}}, \code{\link{influences.CumBH}}, \code{\link{influences.PR}}, 
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, 
\code{\link{influences.CumBH.missingdata}}, \code{\link{influences.PR.missingdata}} and \code{\link{variance}}.
}

\examples{

data(dataexample, package="CaseCohortCoxSurvival")

casecohort  <- dataexample$casecohort # a simulated stratified case cohort
A           <- dataexample$A # matrix with auxiliary variables values in the cohort

Tau1    <- 0 # given time interval for the pure risk
Tau2    <- 8
x       <- c(-1, 1, -0.6) # given covariate profile for the pure risk

# Estimation using the stratified case cohort with weights calibrated on A 

mod.calib <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                   weight = weights.calib, id = id, robust = TRUE)
estimation.calib    <- influences(mod.calib, A = A, calibrated = TRUE, 
                                  Tau1 = Tau1, Tau2 = Tau2, x = x)

# robust variance estimate for the log-relative hazard
robustvariance(estimation.calib$infl.beta)

# robust variance estimate for the cumulative baseline hazard estimate
robustvariance(estimation.calib$infl.Lambda0.Tau1Tau2)

# robust variance estimate for the pure risk estimate
robustvariance(estimation.calib$infl.Pi.x.Tau1Tau2)

}
