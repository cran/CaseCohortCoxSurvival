\name{influences.RH}
\alias{influences.RH}
\title{influences.RH}
\description{
Computes the influences on the log-relative hazard. Can take calibration of the 
design weights into account.
}

\usage{
influences.RH(mod, calibrated = NULL, A = NULL)
}

\arguments{
\item{mod}{a cox model object, result of function coxph.}
\item{calibrated}{are calibrated weights used for the estimation of the
  parameters? If \code{calibrated = TRUE}, the argument below needs to be
  provided. Default is \code{FALSE}.}
\item{A}{\eqn{n \times q} matrix with the values of the auxiliary variables used for
  the calibration of the weights in the whole cohort. Needs to be provided if
  \code{calibrated = TRUE}.}
}

\value{

\code{infl.beta}: matrix with the overall influences on the log-relative hazard estimates.

\code{infl2.beta}: matrix with the phase-two influences on the log-relative hazard estimates. Returned if \code{calibrated = TRUE}.

\code{beta.hat}: vector of length \eqn{p} with log-relative hazard estimates.
}

\details{

\code{influences.RH} works for estimation from a case-cohort with design weights
or calibrated weights (case-cohort consisting of the subcohort and cases not in 
the subcohort, i.e., case-cohort obtained from two phases of sampling).

If covariate information is missing for certain individuals in the phase-two data 
(i.e., case-cohort obtained from three phases of sampling), use \code{\link{influences.RH.missingdata}}.

\code{influence.RH} uses the influence formulas provided in Etievant and Gail 
(2023). More precisely, as in Section 3.2 if \code{calibrated = FALSE}, and as in
Section 4.3 if \code{calibrated = TRUE}. 

If \code{calibrated = FALSE}, the infuences are only provided for the individuals
in the case-cohort. If \code{calibrated = TRUE}, the influences are provided for
all the individuals in the cohort.
}

\references{
Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.
}

\seealso{
\code{\link{estimation}}, \code{\link{estimation.CumBH}}, \code{\link{estimation.PR}},
\code{\link{influences}}, \code{\link{influences.CumBH}}, \code{\link{influences.PR}}, 
\code{\link{influences.missingdata}}, \code{\link{influences.RH.missingdata}},
\code{\link{influences.CumBH.missingdata}},\cr 
\code{\link{influences.PR.missingdata}}, \code{\link{robustvariance}} and  \code{\link{variance}}.
}

\examples{
data(dataexample, package="CaseCohortCoxSurvival")

cohort      <- dataexample$cohort # a simulated cohort
casecohort  <- dataexample$casecohort # a simulated stratified case-cohort
A           <- dataexample$A # matrix with auxiliary variables values in the cohort


# Estimation using the whole cohort

mod.cohort <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = cohort, 
                    robust = TRUE)
estimation.cohort <- influences.RH(mod.cohort)

# print the vector with log-relative hazard estimates
#estimation.cohort$beta.hat

# print the influences on the log-relative hazard estimates
#estimation.cohort$infl.beta

# Estimation using the stratified case-cohort with weights calibrated on A 

mod.calib <- coxph(Surv(times, status) ~ X1 + X2 + X3, data = casecohort, 
                   weight = weights.calib, id = id, robust = TRUE)
estimation.calib    <- influences.RH(mod.calib, A = A, calibrated = TRUE)

# print the influences on the log-relative hazard estimates
#estimation.calib$infl.beta

# print the phase-two influences on the log-relative hazard estimates
#estimation.calib$infl2.beta
}

