
\name{dataexample.missingdata.stratified}
\alias{dataexample.missingdata.stratified}

\title{
Example of case-cohort with stratified sampling of the subcohort and missing covariate information in phase-two data
}
\description{

List with \code{cohort}.

\code{cohort} is a simulated cohort with 20 000 subjects. It contains:

\code{id} is the subject identifier.

\code{X1} is a continuous baseline covariate. Its measurements are only available for subjects in the case-cohort, i.e., with \code{phase3 = 1}.

\code{X2} is a categorical baseline covariate, with categories 0, 1, and 2. It is measured on all cohort subjects.

\code{X3} is a continuous baseline covariate. Its measurements are only available for subjects in the case-cohort.

\code{W} is a baseline categorical variable,  with categories 0, 1, 2, and 3. It depends on predictors of \code{X1} and \code{X2}. It is measured on all cohort subjects.

\code{status} indicates case status.

\code{event.time} gives the event or censoring time. \code{status} indicates whether the subject experienced the event of interest or was censored.

The stratified sampling of the subcohort was based on the 4 strata defined by \code{W}. 97, 294, 300, and 380 subjects were sampled (independently of case status) from the 4 strata, respectively. \code{subcohort} indicates all these subjects included in the subcohort.

The phase-two sample consisted of the subcohort and any other cases not in the subcohort. \code{phase2} indicates all these subjects included in the phase-two sample.

\code{W3} is a baseline binary variable, based on case status. It is measured on all cohort subjects.

The third phase of sampling was stratified based on the 2 strata defined by \code{W3}. Subjects were sampled from the 2 strata with sampling probabilities 0.9 and 0.8. \code{phase3} indicates all these subjects included in the case-cohort (phase-three sample).

\code{strata.n} gives the number of subjects in the stratum in the cohort.

\code{strata.m} gives the number of subjects sampled from each of the 4 phase-two strata to be included in the subcohort (i.e., 97, 294, 300, or 380).

\code{strata.m} and \code{strata.n} would be used to compute the phase-two design weights of non-cases. Because all the cases were included in the phase-two sample, they would be assigned a phase-two design weight of 1.

\code{strata.n.cases} gives the number of cases in each of the 4 phase-two strata in the cohort.

\code{n.cases} gives the number of cases in the entire cohort.

\code{strata.proba.missing} gives the the sampling probablity for the 2 phase-three strata based on \code{W3} and that were used for the third phase of sampling.

\code{weight.true} gives the true design weight (i.e., product of the phase-two and true phase-three design weight).

\code{weight.p2.true} gives true phase-two design weight. They are stratum-specific based on \code{W}.

\code{weight.p3.true} gives the true phase-three design weight. They are stratum-specific based on \code{W3}. \code{weight.p3.true} can be used with argument \code{weights.phase3} of function \code{caseCohortCoxSurvival}, along with argument \code{weights.phase3.type = "design"}.

\code{weight.p3.est} gives the estimated phase-three design weight. They were obtained from \code{W3}, \code{phase2} and \code{phase3}. \code{weight.p3.est} can be used with argument \code{weights.phase3} of function \code{caseCohortCoxSurvival}, along with argument \code{weights.phase3.type = "estimated"}. If in function \code{caseCohortCoxSurvival} \code{weights.phase3 = NULL} but \code{weights.phase3.type = "estimated"}, the phase-three design weights will be estimated from \code{W3}, \code{phase2} and \code{phase3} and should be identical.

\code{weight.est} gives the estimated design weight (i.e., product of the phase-two and estimated phase-three design weight).

}
%\details{ }
\examples{
\donttest{
 data(dataexample.missingdata.stratified, package="CaseCohortCoxSurvival")

 # Display some of the data
 dataexample.missingdata.stratified$cohort[1:5, ]
}}
\keyword{data}
%\seealso{\code{\link{caseCohortCoxSurvival}}}

\references{

Etievant, L., Gail, M. H. (2024). Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Analysis, 30, 572-599.

Etievant, L., Gail, M. H. (2024). Software Application Profile: CaseCohortCoxSurvival: an R package for case-cohort inference for relative hazard and pure risk under the Cox model. Submitted.

}



