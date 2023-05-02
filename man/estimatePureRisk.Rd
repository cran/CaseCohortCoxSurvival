
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
  data(dataexample, package="CaseCohortCoxSurvival")
  
  data <- dataexample$cohort
  
  obj <- caseCohortCoxSurvival(data, "status", "times", cox.phase2 = c("X1", "X2", "X3"), 
                              subcohort = "subcohort", Tau1 = 0, Tau2 = 8)

  # get pure risk estimate for every individual's profile in the cohort
  ret <- estimatePureRisk(obj, data) 

  # get pure risk estimate for one given covariate profile
  ret <- estimatePureRisk(obj, list(X1 = 1, X2 = 2, X3 = 3))
  
}
 
\references{ 
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}