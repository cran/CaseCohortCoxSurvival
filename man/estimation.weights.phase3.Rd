
\name{estimation.weights.phase3}
\alias{estimation.weights.phase3}
\title{estimation.weights.phase3}
\description{
Estimates the weights for the third phase of sampling (due to missingness in phase two).
}

\usage{
estimation.weights.phase3(B.phase3, total.phase2, gamma0 = NULL, niter.max = NULL,
epsilon.stop = NULL)
}


\arguments{

\item{B.phase3}{matrix for the case-cohort (phase-three data), with phase-three
sampling strata indicators. It should have as many columns as phase-three strata
(\eqn{J^{(3)}}), with one 1 per row, to indicate the phase-three stratum position.}

\item{total.phase2}{vector of length \eqn{J^{(3)}} with un-weighted column totals
  for \code{B} in the phase-two data (i.e.,
using all the individuals, even the ones with missing covariate data).}

\item{gamma0}{vector of length \eqn{J^{(3)}} with initial values for \eqn{\gamma}
(Lagrangian multipliers), to be used as seed in the iterative procedure. Default
is (0,...,0).}

\item{niter.max}{maximum number of iterations for the iterative optimization
algorithm. Default is \code{10^4} iterations.}

\item{epsilon.stop}{threshold for the difference between the estimated weighted
total and the total in the whole cohort. If this difference is less than the
value of \code{epsilon.stop}, no more iterations will be performed. Default is \code{10^(-10)}.}
}

\value{

\code{gamma.hat}: vector of length \eqn{J^{(3)}} with final gamma values.

\code{estimated.weights}: vector with the estimated phase-three weights for the
individuals in the case-cohort (phase-three data), computed from \code{B.phase3}
and \code{gamma.hat}.

\code{estimated.total}: vector with the estimated totals, computed from the
 \code{estimated.weights} and \code{B.phase3}.
}

\details{

\code{estimation.weights.phase3} estimates the phase-three sampling weights by solving in
 \eqn{\gamma}

 \deqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} \lbrace \xi_{i,j} V_{i,j}
      \text{exp}( \gamma' B_{i,j}) B_{i,j} - \xi_{i,j} B_{i,j} \rbrace = 0,}
with \eqn{\xi_{i,j}} the phase-two sampling indicator and \eqn{V_{i,j}} the phase-three
sampling indicator of individual \eqn{i} in stratum \eqn{j}, and with
\eqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} \xi_{i,j} B_{i,j}} the total in the
phase-two data. See Etievant and Gail (2024).
The Newton Raphson method is used to solve the optimization problem.

In the end, the estimated weights are given by \eqn{\text{exp}(\hat \gamma' B_{i,j})},
and \eqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} \xi_{i,j} V_{i,j} \text{exp}(\hat \gamma' B_{i,j}) B_{i,j}}
gives the estimated total.
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},\code{\link{influences.CumBH.missingdata}} and \code{\link{influences.PR.missingdata}}.}

\examples{
\donttest{
  data(dataexample.missingdata.stratified, package="CaseCohortCoxSurvival")

  cohort <- dataexample.missingdata.stratified$cohort
  phase2 <- cohort[which(cohort$phase2 == 1),] # the phase-two sample
  casecohort <- cohort[which(cohort$phase3 == 1),] # the stratified case-cohort

  B.phase2 <- cbind(1 * (phase2$W3 == 0), 1 * (phase2$W3 == 1))
  rownames(B.phase2) <- cohort[cohort$phase2 == 1, "id"]
  B.phase3 <- cbind(1 * (casecohort$W3 == 0), 1 * (casecohort$W3 == 1))
  rownames(B.phase3) <- cohort[cohort$phase3 == 1, "id"]
  total.B.phase2 <- colSums(B.phase2)
  J3 <- ncol(B.phase3)

  estimation.weights.p3 <- estimation.weights.phase3(B.phase3 = B.phase3,
                                                  total.phase2 = total.B.phase2,
                                                  gamma0 = rep(0, J3),
                                                  niter.max = 10^(4),
                                                  epsilon.stop = 10^(-10))
}}

