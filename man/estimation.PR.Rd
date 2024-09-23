\name{estimation.PR}
\alias{estimation.PR}
\title{estimation.PR}
\description{
Estimates the pure risk in the time interval [Tau1, Tau2] and for a covariate
profile x, from the log-relative hazard and cumulative baseline hazard values.
}

\usage{
estimation.PR(beta, Lambda0.Tau1Tau2, x = NULL)
}

\arguments{
\item{beta}{vector of length \eqn{p} with log-relative hazard values.}
\item{Lambda0.Tau1Tau2}{cumulative baseline hazard in [Tau1, Tau2].}
\item{x}{vector of length \eqn{p}, specifying the covariate profile considered for the
pure risk. Default is (0,...,0).}
}

\value{
  \code{Pi.x.Tau1Tau2.hat}: pure risk estimate in [Tau1, Tau2] and for covariate profile \eqn{x}.
}

\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.
}

\seealso{
\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{influences}}, \code{\link{influences.RH}},
\code{\link{influences.CumBH}}, \code{\link{influences.PR}},
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}}, \code{\link{influences.CumBH.missingdata}},
and \code{\link{influences.PR.missingdata}}.
}

\examples{
\donttest{
  estimation.PR(beta = c(-0.2, 0.25, -0.3), Lambda0.Tau1Tau2 = 0.03,
                x = c(-1, 1, -0.6))
}
}
