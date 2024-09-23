
\name{dataexample.stratified}
\alias{dataexample.stratified}

\title{
Example of case-cohort with stratified sampling of the subcohort, and set of auxiliary variables
}
\description{

List with \code{cohort} and \code{A}.

\code{cohort} is a simulated cohort with 20 000 subjects. It contains:

\code{id} is the subject identifier.

\code{X1} is a continuous baseline covariate. Its measurements are only available for subjects in the case-cohort, i.e., on subjects with \code{subcohort = 1} and/or \code{status = 1}.

\code{X2} is a categorical baseline covariate, with categories 0, 1, and 2. It is measured on all cohort subjects.

\code{X3} is a continuous baseline covariate. Its measurements are only available for subjects in the case-cohort.

\code{W} is a baseline categorical variable,  with categories 0, 1, 2, and 3. It depends on predictors of \code{X1} and \code{X2}. It is measured on all cohort subjects. The stratified sampling of the subcohort was based on the 4 strata defined by \code{W}.

\code{status} indicates case status.

\code{event.time} gives the event or censoring time. \code{status} indicates whether the subject experienced the event of interest or was censored.

97, 294, 300, and 380 subjects were sampled (independently of case status) from the 4 strata, respectively. \code{subcohort} indicates all these subjects included in the subcohort. The stratified case-cohort (phase-two sample) consists of the subcohort and any other cases not in the subcohort.

\code{strata.n} gives the number of subjects in the stratum in the cohort.

\code{strata.m} gives the number of subjects sampled from each of the 4 strata (i.e., 97, 294, 300, or 380).
\code{strata.m} and \code{strata.n} would be used to compute the stratum-specific design weights of non-cases. Because all the cases were included in the case-cohort, they would be assigned a design weight of 1.

\code{strata.n.cases} gives the number of cases in each of the 4 strata.

\code{n.cases} gives the number of cases in the entire cohort.

\code{X1.proxy} is a continuous baseline covariate. It is a proxy of \code{X1}, with 0.8 correlation. It is measured on all cohort subjects. It can be used for design weights calibration in the argument \code{predictors.cox.phase2} of function \code{caseCohortCoxSurvival}, as one would need to predict \code{X1} on the entire cohort.

\code{X3.proxy} is a continuous baseline covariate. It is a proxy of \code{X3}, with 0.8 correlation. It is measured on all cohort subjects. It can be used for design weights calibration in the argument \code{predictors.cox.phase2} of function \code{caseCohortCoxSurvival}, as one would need to predict \code{X3} on the entire cohort.

\code{X1.pred} is a prediction of \code{X1}, available for all cohort subjects. The predictions were obtained by weighted linear regression on \code{X1.proxy} and \code{W}, with the design weights.

\code{X3.pred} is a prediction of \code{X3}, available for all cohort subjects. The predictions were obtained by weighted linear regression on \code{X1.proxy},  \code{X2}, and \code{X3.proxy}, with the design weights.

\code{A} contains auxiliary variables, obtained as proposed by Breslow et al. (2009) and Shin et al. (2020). \code{A} can be used with argument \code{aux.var} of function \code{caseCohortCoxSurvival}.

Predictions of \code{X1} were obtained by weighted linear regression on \code{X1.proxy} and \code{W}, with the design weights. Predictions of \code{X3} were obtained by weighted linear regression on \code{X1.proxy},  \code{X2}, and \code{X3.proxy}, with the design weights. Then the Cox model with \code{X2} and the predicted values of \code{X1} and \code{X3} (available for all cohort subjects) was run. \code{A.X1}, \code{A.X2}, and \code{A.X3} contain the influences on the estimated log-RHs (available for all cohort subjects).

Second, design weights were then calibrated based on \code{A.1}, \code{A.X1}, \code{A.X2}, and \code{A.X3}, with \code{A.1} that is identically equal to 1. The log-RH parameter was then estimated from the case-cohort data with these calibrated weights. Finally, the log-RH estimate was used with \code{X2} and the predicted values of \code{X1} and \code{X3} (available for all cohort subjects), and exponentiated. \code{A.Shin} contains the product of this quantity with the total follow-up time on interval (0,8].

}
%\details{ }
\examples{
\donttest{
 data(dataexample.stratified, package="CaseCohortCoxSurvival")

 # Display some of the data
 dataexample.stratified$cohort[1:5, ]

 dataexample.stratified$A[1:5, ] # auxiliary variable values in the cohort
}}
\keyword{data}
%\seealso{\code{\link{caseCohortCoxSurvival}}}
\references{
Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.

Etievant, L., Gail, M. H. (2024). Software Application Profile: CaseCohortCoxSurvival: an R package for case-cohort inference for relative hazard and pure risk under the Cox model. Submitted.

Shin Y.E., Pfeiffer R.M., Graubard B.I., Gail M.H. (2020) Weight calibration to improve
the efficiency of pure risk estimates from case-control samples nested in a cohort.
Biometrics, 76, 1087-1097

Breslow, N.E., Lumley, T., Ballantyne, C.M., Chambless, L.E. and Kulich, M. (2009).
Improved Horvitz-Thompson Estimation of Model Parameters from Two-phase Stratified Samples:
Applications in Epidemiology. Statistics in Biosciences, 1, 32-49.
}



