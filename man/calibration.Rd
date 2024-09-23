
\name{calibration}
\alias{calibration}
\title{calibration}
\description{
  Calibrates the design weights using the raking procedure.
}
\usage{
calibration(A.phase2, design.weights, total, eta0 = NULL, niter.max = NULL,
epsilon.stop = NULL)
}
\arguments{
  \item{A.phase2}{matrix with the values of the q auxiliary variables to be used for
                  the calibration of the weights in the case-cohort (phase-two data).}
  \item{design.weights}{design weights to be calibrated.}
  \item{total}{vector of length q with un-weighted auxiliary variable totals in the whole cohort.}
  \item{eta0}{vector of length q with initial values for \code{eta} (the Lagrangian multipliers),
              to be used as seed in the iterative procedure. Default is (0, ... ,0).}
  \item{niter.max}{maximum number of iterations for the iterative optimization algorithm.
                   Default is 10^4 iterations.}
  \item{epsilon.stop}{threshold for the difference between the estimated weighted total
                 and the total in the whole cohort.
                 If this difference is less than the value of epsilon.stop,
                 no more iterations will be performed.
                 Default is 10^(-10).}
}

\value{
  \code{eta.hat}: vector of length q with final eta values.

  \code{calibrated.weights}: vector with the calibrated weights for the individuals
             in the case-cohort (phase-two data), computed from \code{design.weights},
             \code{A.phase2} and \code{eta.hat}.

  \code{estimated.total}: vector with the estimated totals, computed from the
                 \code{calibrated.weights} and \code{A.phase2}.
}
\details{
  Calibration matches the weighted total of the auxiliary variables in the case-cohort
  (with calibrated weights), to the un-weighted auxiliary variables total in the whole cohort.
  In other words, it solves in \eqn{\eta}
  \eqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} \lbrace \xi_{i,j} w_{i,j}
      \text{exp}(\eta' A_{i,j}) A_{i,j} - A_{i,j} \rbrace = 0},
  with \eqn{\xi_{i,j}} the sampling indicator and \eqn{w_{i,j}} the design weight of
  individual \eqn{i} in stratum \eqn{j}, and with
  \eqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} A_{i,j}}
 the total in the whole cohort. See Etievant and Gail (2024). The Newton Raphson method is used to solve the optimization problem.
In the end, the calibrated weights of the case-cohort individuals are given by
 \eqn{w_{i,j} \text{exp}(\hat \eta' A_{i,j})}, and
 \eqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} \lbrace \xi_{i,j} w_{i,j}
     \text{exp}(\hat \eta' A_{i,j}) A_{i,j} \rbrace}
 gives the estimated total.
}

\seealso{ \code{\link{auxiliary.construction}}, \code{\link{influences}}, \code{\link{influences.RH}}, \code{\link{influences.CumBH}} and \code{\link{influences.PR}}. }

\examples{
\donttest{
  data(dataexample.stratified, package="CaseCohortCoxSurvival")
  cohort <- dataexample.stratified$cohort
  casecohort <- cohort[which(cohort$status == 1 |
                       cohort$subcohort == 1),] # the stratified case-cohort
  casecohort$weight <- casecohort$strata.n / casecohort$strata.m
  casecohort$weight[which(casecohort$status == 1)] <- 1

  A <- dataexample.stratified$A # auxiliary variables values in the cohort
  indiv.phase2 <- casecohort$id
  q <- ncol(A)
  total <- colSums(A)
  A.phase2 <- A[indiv.phase2,]
  calib <- calibration(A.phase2 = A[indiv.phase2,], design.weights = casecohort$weight,
                       total = total, eta0 = rep(0, q), niter.max = 10^3, epsilon.stop = 10^(-10))
  #calib$calibrated.weights # print calibrated weights
}
}

\references{
Deville, J.C. and Sarndal, C.E. (1992). Calibration Estimators in Survey Sampling.
 Journal of the American Statistical Association, 87, 376-382.

Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}
