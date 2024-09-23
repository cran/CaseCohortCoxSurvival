
\name{estimatePureRisk}
\alias{estimatePureRisk}
\title{estimatePureRisk}
\description{
  Computes pure risk estimates and variances for new covariate values.
}
\usage{
estimatePureRisk(obj, x)
}
\arguments{
  \item{obj}{Return object from \code{\link{caseCohortCoxSurvival}}.}
  \item{x}{Data frame or a list containing values of the covariates that were used
           when \code{\link{caseCohortCoxSurvival}} was called,
           and for which the pure risk is to be estimated.}
}
\value{
  A list containing:
  \itemize{
   \item \code{var} Matrix of pure risk estimates in [Tau1, Tau2] and variance estimates
   \item \code{var.estimated} Matrix of pure risk estimates in [Tau1, Tau2] and variance estimates
                                 when the phase-three weights are estimated
   \item \code{var.design} Matrix of pure risk estimates in [Tau1, Tau2] and variance estimates
                                 when the phase-three weights are known
  }
  Depending on the analysis run, some of the above objects will be NULL.
}
%\details{}

\seealso{ \code{\link{caseCohortCoxSurvival}} }

\examples{
\donttest{
  data(dataexample.stratified, package="CaseCohortCoxSurvival")

  data <- dataexample.stratified$cohort
  cov1 <- "X2"
  cov2 <- c("X1", "X3")

  obj <- caseCohortCoxSurvival(data = data, status = "status",
                               time = "event.time", cox.phase1 = cov1,
                               cox.phase2 = cov2, strata = "W",
                               subcohort = "subcohort", Tau1 = 0, Tau2 = 8)

  # get pure risk estimate for every individual's profile in the cohort
  ret <- estimatePureRisk(obj, data)

  # get pure risk estimate for one given covariate profile
  ret <- estimatePureRisk(obj, list(X1 = 1, X2 = -1, X3 = 0.6))

  # get pure risk estimates for two given covariate profiles
  pr1 <- as.data.frame(cbind(X1 = -1, X2 = 1, X3 = -0.6))
  pr2 <- as.data.frame(cbind(X1 = 1, X2 = -1, X3 = 0.6))
  ret <- estimatePureRisk(obj, rbind(pr1, pr2))
  ret$var
}}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}
