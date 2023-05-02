
\name{caseCohortCoxSurvival}
\alias{caseCohortCoxSurvival}
\title{Parameter estimation and variance for case-cohort analyses under the Cox model}
\description{
  Function for estimating parameters (log-relative hazard, baseline hazards, 
  cumulative baseline hazard, pure risks) and their variance 
  (robust or the one accounting for sampling features) from cohort or case-cohort data,
   under the Cox model.
}
\usage{
 caseCohortCoxSurvival(data, status, time, 
                       cox.phase1=NULL, cox.phase2=NULL, other.covars=NULL, 
                       strata=NULL, weights.phase2=NULL, calibrated=FALSE, 
                       subcohort=NULL, subcohort.strata.counts=NULL,
                       predict=TRUE, predicted.cox.phase2=NULL,
                       predictors.cox.phase2=NULL,
                       aux.vars=NULL, aux.method="Shin",
                       phase3=NULL, strata.phase3=NULL,
                       weights.phase3=NULL, weights.phase3.type="both",
                       Tau1=NULL, Tau2=NULL, x=NULL,
                       weights.op=NULL, print=1)
}
\arguments{
  \item{data}{Data frame containing the cohort and all variables needed for the analysis.}
  \item{status}{Column name in \code{data} giving the case status for each
                individual in the cohort. This variable
                must be coded as 0 for non-cases and 1 for cases.}
  \item{time}{Column name(s) in \code{data} giving the time to event for each individual
              in the case-cohort. One variable
    is required for a time-on-study time scale, two variables for age-scale,
    with the first variable as the start age and second as the end age.}
  \item{cox.phase1}{Column name(s) in \code{data} giving the Cox model covariates 
                    measured on the entire cohort.
               See covariates and prediction in details.}
  \item{cox.phase2}{Column name(s) in \code{data} giving the Cox model covariates 
                    measured only on phase-two individuals.
               See covariates and prediction in details.}
  \item{other.covars}{Column name(s) in data giving other covariates 
                 measured on the entire cohort that might be useful, 
             alone or in combination with \code{cox.phase1}, 
          if predicted values of the phase-two covariates (\code{cox.phase2}) 
          need to be obtained on the whole cohort for the weight calibration.}
  \item{strata}{NULL or column name in data with the stratum value for each 
                individual in the cohort. The number of strata used for the 
                sampling of the subcohort equals the number of different stratum values. 
                For example, a stratum variable might take values 0,1,2,3 or 4. 
                The default is NULL.}
  \item{weights.phase2}{NULL or column name in data giving the phase-two design
        weights for each individual in the cohort. 
       For a whole cohort analysis (see \code{subcohort} below), 
      weights are not used in the \code{\link[survival]{coxph}} call. 
     If NULL but subcohort is not NULL, \code{subcohort.strata.counts}
    will be used to estimate \code{weights.phase2}.
                 The default is NULL.}
  \item{calibrated}{TRUE or FALSE to calibrate the \code{weights}.
                    Calibrated weights will be computed using the function 
                    \code{\link{calibration}}.   
                    If TRUE, then \code{phase3} (below) will be set to NULL.
                    See calibration in details.
                 The default is FALSE.}
  \item{subcohort}{NULL or column name in \code{data} giving the indicators 
           of membership in the subcohort. 
           The indicators are 1 if the individual belongs to the subcohort 
           and 0 otherwise. Some cases might be in the subcohort and others not.  
           If NULL, then a whole cohort analysis will be performed.
                 The default is NULL.}
   \item{subcohort.strata.counts}{NULL or a list of the number of individuals sampled 
            into the subcohort from each stratum of strata. 
            The names in the list must be the strata values and the length of the 
            list must be equal to the number of strata.
            If NULL, then the count for each stratum is estimated by the number of 
            subcohort individuals in each stratum. 
                 The default is NULL.}
  \item{predict}{TRUE or FALSE to predict the phase-two covariates using
                 \code{predictors.cox.phase2}. This option is not used if 
            \code{calibrated=FALSE}. If \code{calibrated=TRUE}, \code{aux.vars=NULL}
            and \code{predict=FALSE}, then \code{predicted.cox.phase2} must be specified. 
            See covariates and prediction in details. 
            This option is only used when \code{calibrated=TRUE}, \code{aux.vars=NULL}
            and \code{predicted.cox.phase2=NULL}. 
                 The default is TRUE.}
  \item{predicted.cox.phase2}{NULL or a named list giving the predicted values of the 
            phase-two covariates (\code{cox.phase2}) on the whole cohort. For example,
            if the phase-two covariates are \code{X1} and \code{X2}, then the list is
            of the form \code{list(X1=X1.pred, X2=X2.pred)}, where \code{X1.pred} and
            \code{X2.pred} are the predictions of \code{X1} and \code{X2} respectively.
           This option is only used when \code{calibrated=TRUE} and \code{aux.vars=NULL}. 
           If \code{calibrated=TRUE}, \code{aux.vars=NULL} and \code{predict=FALSE}, 
           then \code{predicted.cox.phase2} must be specified and must not contain missing values. 
           The default is NULL.}
  \item{predictors.cox.phase2}{NULL, a vector, or a list specifying the columns in data
        to use as predictor variables for obtaining the predicted values
        on the whole cohort for the phase-two covariates (\code{cox.phase2}). 
        A list allows for different proxy variables to be used for the different 
        phase-two covariates. The selected predictor variables must be from among 
        \code{cox.phase1} and \code{other.covars}.  
        See examples and covariates and prediction in details. 
        If NULL, then the phase-two covariates will be predicted using \code{cox.phase1}
        and \code{other.covars}. If NULL, \code{cox.phase1=NULL} and \code{other.covars=NULL}, 
        then the calibrated analysis will not be performed. 
        This option is only used when \code{calibrated=TRUE}, \code{aux.vars=NULL}, \cr
        \code{predicted.cox.phase2=NULL} and \code{predict=TRUE}. The default is NULL.}
  \item{aux.vars}{NULL or column name(s) in data giving the auxiliary variables for 
          each individual in the cohort. This option is only used when 
          \code{calibrated=TRUE}. If NULL, then auxiliary variables will be constructed 
          using method Breslow or Shin and predicted values on the whole cohort for the 
          phase-two covariates (see \code{aux.method}, \code{predict}, \code{predicted.cox.phase2}
          and \code{predictors.cox.phase2}). \code{aux.vars} must not contain missing values.
         The default is NULL.}
  \item{aux.method}{"Breslow", or "Shin" to specify the algorithm to construct the
                    auxiliary variables. This option is only used if \code{aux.vars=NULL} and
                    \code{calibrated=TRUE}.
                  The default is "Shin".}
  \item{phase3}{NULL or column name in data giving the indicators of membership in the in 
             the phase-three sample. The indicators are 1 if the individual belongs to the 
             phase-three sample and 0 otherwise. All individuals in the phase-three sample 
            must also belong to the phase-two sample. 
          This option is not used if \code{calibrated=TRUE}. 
          The default is NULL.}
  \item{strata.phase3}{NULL or column name in \code{data} giving the phase-three stratification
              for each individual in phase-two. The number of strata used for the third phase 
             of sampling equals the number of different stratum values. 
                 The default is NULL.}
  \item{weights.phase3}{NULL or column name in \code{data} giving the phase-three design weights
                 for each individual in phase-two.
                 If NULL but \code{phase3} is not NULL, then \code{phase3} and \code{subcohort}
                 will be used to estimate \code{weights.phase3} (see details in 
                 \code{\link{estimation.weights.phase3}}).
                 The default is NULL.}
  \item{weights.phase3.type}{One of NULL, "design", "estimated", or "both" to specify whether
                         the phase-three weights are design weights (known), or
                         to be estimated. The variance estimation differs for estimated
                         and design weights. If set to "both", then both variance estimates
                         will be computed. If not NULL, then only the first letter is matched for this option.
                 The default is "both".}
  \item{Tau1}{NULL or left bound of the time interval considered for the cumulative 
              baseline hazard and the pure risk. 
               If NULL, then the first event time is used.} 
  \item{Tau2}{NULL or right bound of the time interval considered for the cumulative 
              baseline hazard and the pure risk. 
               If NULL, then the last event time is used.} 
  \item{x}{Data frame containing \code{cox.phase1} and \code{cox.phase2} 
                   variables for which pure risk is estimated. 
                   The default is NULL so that no pure risk estimates will be computed.}
  \item{weights.op}{NULL or a list of options for calibration of phase-two design weights
                   or estimating phase-three design weights.
                 The available options are \code{niter.max}, and
                 \code{epsilon.stop} (see \code{\link{calibration}} or 
                 \code{\link{estimation.weights.phase3}}). 
                 The default is NULL.}
  \item{print}{0-3 to print information as the analysis is performed. 
               The larger the value, the more information will be printed. To not
               print any information, set \code{print = 0}.
              The default is 1.} 
}
\value{
  A list with class \code{casecohortcoxsurv} containing:
  \itemize{
   \item \code{beta} Estimated log-relative hazard estimates
   \item \code{Lambda0} Cumulative baseline hazard estimate in [Tau1, Tau2]
   \item \code{beta.var} Influence-based variance estimate for \code{beta} 
   \item \code{Lambda0.var} Influence-based variance estimate for \code{Lambda0}
   \item \code{beta.var.estimated} Influence-based variance estimate for \code{beta} with estimated
                                   phase-three weights 
   \item \code{Lambda0.var.estimated} Influence-based variance estimate for \code{Lambda0} with estimated
                                   phase-three weights 
   \item \code{beta.var.design} Influence-based variance estimate for \code{beta} with design
                                   phase-three weights
   \item \code{Lambda0.var.design} Influence-based variance estimate for \code{Lambda0} with design
                                   phase-three weights 
   \item \code{beta.robustvar} Robust variance estimate for \code{beta} 
   \item \code{Lambda0.robustvar} Robust variance estimate for \code{Lambda0}
   \item \code{beta.robustvar.estimated} Robust variance estimate for \code{beta} with estimated
                                   phase-three weights
   \item \code{Lambda0.robustvar.estimated} Robust variance estimate for \code{Lambda0} with estimated
                                   phase-three weights
   \item \code{beta.robustvar.design} Robust variance estimate for \code{beta} with design
                                   phase-three weights
   \item \code{Lambda0.robustvar.design} Robust variance estimate for \code{Lambda0} with design
                                   phase-three weights
   \item \code{Pi.var} Matrix of pure risk estimates in [Tau1, Tau2] and variance estimates
   \item \code{Pi.var.estimated} Matrix of pure risk estimates in [Tau1, Tau2] and variance estimates
                                 with estimated phase-three weights
   \item \code{Pi.var.design} Matrix of pure risk estimates in [Tau1, Tau2] and variance estimates
                                 with design phase-three weights 
   \item \code{coxph.fit} Return object from \code{\link[survival]{coxph}} of the model fit
   \item \code{changed.times} Matrix of original and new event times for individuals who had their event times
                              changed due to ties. Will be NULL if event times were not changed.
   \item \code{args} List containing the values of the input arguments (except \code{data})
   \item \code{risk.obj} List containing objects needed to compute pure risk estimates and variances
                         for a different set of data
 }

}
\details{
  The different scenarios covered by the function are:\cr
  1) Whole cohort (\code{subcohort = NULL})\cr

  2) (stratified) case-cohort (= stratified phase-two sample with no missing covariate data)\cr
  a. With design weights (\code{subcohort}, \code{strata}, \code{calibrated = FALSE})\cr
  b. With calibrated weights and proxies to predict phase-two covariates and the
     auxiliary variables (\code{subcohort}, \code{strata}, \code{calibrated=TRUE},
     \code{predict=TRUE}, \code{predictors.cox.phase2}, \code{aux.method}) \cr
  c. With calibrated weights and externally supplied predicted values of phase-two covariates
     (\code{calibrated=TRUE}, \code{strata}, \code{predict=FALSE}, \code{predicted.cox.phase2}) \cr

  3) (unstratified) case-cohort (= unstratified phase-two sample with no missing covariate data)\cr
  a. With design weights (\code{subcohort}, \code{strata=NULL}, \code{calibrated=FALSE})\cr
  b. With calibrated weights and proxies to predict phase-two covariates and obtain the
     auxiliary variables (\code{subcohort}, \code{strata=NULL}, \code{calibrated=TRUE},
     \code{predict=TRUE}, \code{predictors.cox.phase2}, \code{aux.method}) \cr
  c. With calibrated weights and externally supplied predicted values of phase-two covariates
     (\code{calibrated=TRUE}, \code{strata=NULL}, \code{predict=FALSE}, \code{predicted.cox.phase2}) \cr

  4) Case-cohort (= phase-three sample, because of missing covariate information in phase-two 
             data, with stratified or unstratified phase-two sampling)\cr
  a. With known phase-three design weights (\code{subcohort}, \code{strata}, \code{phase3},
       \code{strata.phase3}, \cr
       \code{weights.phase3.type="design"}) \cr
  b. With estimated phase-three design weights (\code{subcohort}, \code{strata}, \code{phase3},
       \code{strata.phase3}, \cr
       \code{weights.phase3.type="estimated"})\cr

  covariates and prediction \cr
  Prediction of phase-two covariates is performed when \code{calibrated = TRUE}, \code{predict = TRUE},
  \code{aux.vars = NULL} and \code{predicted.cox.phase2 = NULL}. If \code{predictors.cox.phase2 = NULL},
  all the covariates measured on the entire cohort will be used for the prediction 
  (see \code{cox.phase1} and \code{other.covars}). 
  Prediction of phase-two covariates is performed by linear regression for a continuous variable, 
     logistic regression for a binary variable and the function \code{\link[nnet]{multinom}} for a 
     categorical variable. Dummy variables should not be used for categorical covariates,
     because independent logistic (or linear) regressions will be performed using the dummy variables. \cr
  Alternatively, predicted values of phase-two covariates on the whole cohort can be specified with 
    \code{predicted.cox.phase2}.
  
  calibration \cr
  Calibrating the design weights against some informative auxiliary variables, 
  measured on all cohort members, can increase efficiency. 
  When \code{calibrated = TRUE}, the user can either provide the auxiliary variables
  (\code{aux.vars}), or let the driver function build the auxiliary variables (\code{aux.method}).
  Construction of the auxiliary variables follows Breslow et al. (2009) or Shin et al. (2020)
  (see \code{aux.method}), and relies on predictions of the phase-two covariates for all members 
  of the cohort (see covariates and prediction above). 
  The auxiliary variables are given by (i) the influences for the log-relative hazard parameters 
  estimated from the Cox model with imputed cohort data; and (ii) the products of total 
  follow-up time (on the time interval for which pure risk is to be estimated) with the estimated 
  relative hazard for the imputed cohort data, where the log-relative hazard parameters are 
  estimated from the Cox model with case-cohort data and weights calibrated with (i). 
  When \code{aux.method = Breslow}, calibration of the design weights is against (i), 
  as proposed by Breslow et al. (2009) to improve efficiency of case-cohort estimates 
  of relative hazard. When \code{aux.method = Shin}, calibration is against (i) and (ii), 
  as proposed by Shin et al. (2020) to improve efficiency of relative hazard and pure risk 
  estimates under the nested case-control design.


  Note \cr
  If \code{subcohort = NULL}, then a whole cohort analysis will be run and only robust variance estimates
           will be computed.\cr
  
}

%\seealso{ \code{\link{calibration}} }

\examples{
  data(dataexample.missingdata, package="CaseCohortCoxSurvival")

  data <- dataexample.missingdata$cohort
  cov1 <- "X1"
  cov2 <- c("X2", "X3")

  # Whole cohort, get pure risk estimate for every individual's profile in the cohort
  # Only robust variance estimates are computed for a whole cohort analysis.
  caseCohortCoxSurvival(data, "status", "times", cox.phase1 = cov1, x = data)

  # Stratified case-cohort analysis with missing covariate information in
  #   phase-two data, with phase-three strata
  caseCohortCoxSurvival(data, "status", "times", cox.phase1 = cov1, 
                      cox.phase2 = cov2, strata = "W", subcohort = "subcohort",  
                      phase3 = "phase3", strata.phase3 = "W3")
           
           
  data(dataexample, package="CaseCohortCoxSurvival")

  data <- dataexample$cohort
  cov2 <- c("X1", "X2", "X3")

  # Stratified case-cohort (phase-two) analysis with weight calibration and default 
  #   proxies to predict the phase-two covariates.
  caseCohortCoxSurvival(data, "status", "times", cox.phase2 = cov2, strata = "W", 
                        subcohort = "subcohort", calibrated = TRUE)

  \donttest{

  # Stratified case-cohort (phase-two) analysis with weight calibration specifying 
  #  a different set of proxy variables to predict each phase-two covariate.
  caseCohortCoxSurvival(data, "status", "times", cox.phase2 = cov2, 
                      strata = "W", subcohort = "subcohort", calibrated = TRUE, 
                      predictors.cox.phase2 = list(X1 = c("X1.proxy"), 
                      X2 = c("X1.proxy", "X2.proxy", "W"), X3 = c("X1.proxy", "X3.proxy")))

  # Stratified case-cohort (phase-two) analysis with weight calibration, get pure
  #   risk estimate for one given covariate profile.
  caseCohortCoxSurvival(data, "status", "times", cox.phase2 = cov2, 
                      strata = "W", subcohort = "subcohort", calibrated = TRUE, 
                      predictors.cox.phase2=list(X1 = c("X1.proxy"), 
                      X2 = c("X1.proxy", "X2.proxy", "W"), X3 = c("X1.proxy", "X3.proxy")),
                      x = list(X1 = 1, X2 = 1, X3 = 0.6), Tau1 = 0, Tau2 = 8)

  # Set the correct sampling counts in phase-two for each level of strata.
  # The strata variable W has levels 0-3.
  caseCohortCoxSurvival(data, "status", "times", cox.phase2 = cov2, 
                      subcohort = "subcohort", calibrated = TRUE, strata = "W",
                      subcohort.strata.counts = list("0"=129, "1"=313, "2"=308, "3"=311))
  }
}
 
\references{ 

Etievant, L., Gail, M.H. (2023). Cox model inference for relative hazard and 
pure risk from stratified weight-calibrated case-cohort data. Submitted.

Shin Y.E., Pfeiffer R.M., Graubard B.I., Gail M.H. (2020) Weight calibration to improve
the efficiency of pure risk estimates from case-control samples nested in a cohort.
Biometrics, 76, 1087-1097

Breslow, N.E., Lumley, T., Ballantyne, C.M., Chambless, L.E. and Kulich, M. (2009). 
Improved Horvitz-Thompson Estimation of Model Parameters from Two-phase Stratified Samples: 
Applications in Epidemiology. Statistics in Biosciences, 1, 32-49.

}