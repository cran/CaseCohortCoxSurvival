
\name{estimation.weights.phase3}
\alias{estimation.weights.phase3}
\title{estimation.weights.phase3}
\description{
Estimates the weights for the third phase of sampling (due to missingness).
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
phase-two data. See Section 5.2 in Etievant and Gail (2023). 
The Newton Raphson method is used to solve the optimization problem.

In the end, the estimated weights are given by \eqn{\text{exp}(\hat \gamma' B_{i,j})},
and \eqn{\sum_{j=1}^J \sum_{i=1}^{n^{(j)}} \xi_{i,j} V_{i,j} \text{exp}(\hat \gamma' B_{i,j}) B_{i,j}}
gives the estimated total. 
}

\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},
\code{\link{influences.CumBH.missingdata}} and \code{\link{influences.PR.missingdata}}.
}

\examples{

data(dataexample.missingdata, package="CaseCohortCoxSurvival")

casecohort      <- dataexample.missingdata$casecohort # a simulated stratified case-cohort
# phase-two data: dataexample.missingdata$casecohort.phase2 
B.phase2        <- dataexample.missingdata$B.phase2
B.phase3        <- dataexample.missingdata$B.phase3

total.B.phase2  <- colSums(B.phase2)
J3              <- ncol(B.phase3)

estimation.weights.p3 <- estimation.weights.phase3(B.phase3 = B.phase3, 
                                        total.phase2 = total.B.phase2, 
                                        gamma0 = rep(0, J3), 
                                        niter.max = 10^(4),
                                        epsilon.stop = 10^(-10))

# print estimated phase-three weights
#estimation.weights.p3$estimated.weights
}

