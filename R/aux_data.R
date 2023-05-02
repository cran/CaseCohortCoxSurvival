getAuxData <- function(data, obj) {

  if (!obj$phase2Flag) return(NULL)
  if (!obj$calibrated) return(NULL)
  print <- obj$print

  if (!obj$predict) {
    aux.vars <- obj[["aux.vars", exact=TRUE]]
    if (!is.null(aux.vars)) {
      if (print) cat("NOTE: aux.vars will be used as auxiliary variables for weight calibration.\n")
      ret <- getModelMatrix(data, aux.vars, rem.intercept=TRUE)
    } else {
     if (print) cat("NOTE: predicted.cox.phase2 variables will be used to build auxiliary variables for weight calibration.\n") 
     ret <- cc_getAuxData(data, obj) 
    }
  } else {
    ret <- cc_getAuxData(data, obj)
  }
  ret
}

cc_getAuxData <- function(data, obj) {

    if (obj$print) {
      if (obj$shin.method) {
        cat("Building the auxiliary variables proposed by Shin.\n")
      } else {
        cat("Building the auxiliary variables proposed by Breslow.\n")
      }
    }

    # Get predicted data 
    data.pred <- getPredictedData(data, obj) 

    # Fit cox Model on cox.phase1 and imputed cox.phase2
    weights  <- data.pred[, obj$weights.phase2, drop=TRUE]
    ids      <- data.pred[, obj$id.var, drop=TRUE]
    covars   <- obj[["covars.orig", exact=TRUE]]
    cntl     <- obj$coxph.control
    if (!length(covars)) stop("INTERNAL CODING ERROR 1")
    form     <- as.formula(getCoxFormula(obj$status, obj$time, covars))
    fit      <- coxph(form, data=data.pred, robust=TRUE, id=ids, control=cntl)
    if (obj$print > 2) {
      cat("***************************************************\n")
      cat("Fitted Cox model with predicted phase-2 covariates:\n")
      cat("***************************************************\n")
      print(fit)
      cat("\n") 
    }
 
    # Check the fitted coefficients, if any missing, then drop and refit
    beta <- fit$coefficients
    tmp  <- is.na(beta)
    if (any(tmp)) {
      nms    <- names(beta)[tmp]
      tmp    <- !(covars %in% nms)
      covars <- covars[tmp] 

      # If no parms could be estimated, stop
      if (!length(covars)) stop("ERROR fitting model with phase-two covariates")

      # Issue a warning if some parms were estimated
      str <- getQuotedVarStr(nms)
      msg <- paste0("In the auxiliary variable construction, the parameters ", str, 
                    " were not estimated in the Cox model adjusting",
                    " for cox.phase1 and predicted cox.phase2 covariates.",
                    " To view the fit, set print=3\n")
      warning(msg)
      form   <- as.formula(getCoxFormula(obj$status, obj$time, covars))
      fit    <- coxph(form, data=data.pred, robust=TRUE, weights=weights, id=ids, control=cntl)
      beta   <- fit$coefficients
    } 
    if (any(is.na(beta))) stop("ERROR fitting model with cox.phase1 and predicted cox.phase2 covariates")
    resid      <- residuals(fit, type="dfbeta", weighted=TRUE)
    fit.detail <- coxph.detail(fit, riskmat=FALSE)

    tmp <- cc_aux.constr(fit.detail, beta, weights, resid, obj$tau1, obj$tau2)
    rm(fit.detail, fit, resid, weights)
    gc()

    # Return A.RH or A.CumBH
    ret <- cbind(1, tmp$A.RH)  # Breslow

    # Shin's method
    if (obj$shin.method) {
      tmp <- cc_aux.shin(data, data.pred, obj, ret) 
      ret <- cbind(ret, tmp) 
    }

    ret
}

cc_aux.shin <- function(data, data.pred, obj, A) {

  # Use covars.orig now because of change in prediction
  # Building the auxiliary variables proposed by Shin et al. (2020) --------------

  Tau1   <- obj$tau1
  Tau2   <- obj$tau2
  timev  <- obj$time
  ntimev <- length(timev)
  if (ntimev == 1) {
    time.on.study <- pmax(pmin(Tau2, data[, timev, drop=TRUE]) - Tau1, 0)
  } else {
    vec1          <- data[, timev[1], drop=TRUE] 
    vec2          <- data[, timev[2], drop=TRUE] 
    time.on.study <- pmax(pmin(Tau2, vec2) - pmax(Tau1, vec1),  0)
  }

  total        <- colSums(A)
  indiv.phase2 <- data[, obj$phase2, drop=TRUE] %in% 1
  wgtop        <- obj$weights.op
  wgtv         <- obj$weights.phase2
  tmp          <- calibration(A.phase2 = A[indiv.phase2,], 
                             design.weights = data[indiv.phase2, wgtv, drop=TRUE], 
                             total = total, eta0 = rep(0, ncol(A)), 
                             niter.max=wgtop$niter.max, epsilon.stop=wgtop$epsilon.stop)
  weights.calib1 <- as.vector(tmp$calibrated.weights)
  ids            <- data[indiv.phase2, obj$id.var, drop=TRUE]
  covars         <- obj$covars.orig
  form           <- as.formula(getCoxFormula(obj$status, obj$time, covars))
  mod.calib1     <- coxph(form, data=data[indiv.phase2, , drop=FALSE], 
                          weights=weights.calib1, robust=TRUE, id=ids, control=obj$coxph.control)
  coef           <- mod.calib1$coefficients
  mat            <- getModelMatrix(data.pred, covars, rem.intercept=TRUE)
  cnames         <- names(coef)
  mnames         <- colnames(mat)
  nms            <- intersect(cnames, mnames)
  m              <- length(nms)
  if (!m) stop("ERROR: names do not agree")
  if ((m != length(cnames)) || (m != length(mnames))) {
    warning("Not all parameter names agree between coxph model and predicted data")
  }
  coef    <- coef[nms]
  mat     <- mat[, nms, drop=FALSE]
  coef    <- matrix(coef, ncol=1)
  A.Shin  <- time.on.study*exp(mat %*% coef)

  A.Shin
}

cc_aux.constr <- function(mod.detail, beta.hat, weights, resid, Tau1, Tau2) {

  
  X                 <- getFitDetailX(mod.detail)
  n                 <- nrow(X)
  p                 <- ncol(X) 
  mod.detail.time   <- mod.detail$time
  number.times      <- length(mod.detail.time)
  tmp               <- getStatusTimesCoxDetail(mod.detail$y)
  y                 <- tmp$status
  time              <- tmp$time
  time2             <- tmp[["time2", exact=TRUE]]
  if (is.null(time2)) {
    obstimes           <- cc_getTimeId(time, mod.detail.time)
    obstimes0          <- NULL
  } else {
    obstimes           <- cc_getTimeId(time2, mod.detail.time)
    obstimes0          <- cc_getTime0Id(time, mod.detail.time)
    if (any(obstimes < obstimes0)) stop("INTERNAL CODING ERROR")
  }
  if (is.null(Tau1)) Tau1 <- floor(min(mod.detail.time))
  if (is.null(Tau2)) Tau2 <- floor(max(mod.detail.time))
  Tau1Tau2.times    <- which((Tau1 < mod.detail.time) & (mod.detail.time <= Tau2)) 
  if (!length(Tau1Tau2.times)) stop("ERROR: no event times in (tau1, tau2], adjust tau1 and/or tau2")
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  S0t               <- cc_getS0t(obstimes, obstimes0, number.times, exp.X.weighted)
  if (any(S0t == 0)) stop("ERROR: S0t = 0")
  S1t               <- cc_getS1t(obstimes, obstimes0, exp.X.weighted, X, number.times)
  rm(X, time, time2)
  gc()

  lambda0.t.hat         <- cc_getdNtColSums(obstimes, number.times, y)/S0t
  infl.Lambda0.Tau1Tau2 <- cc_infl_lambda0_tau12_noCalib(resid, S1t, 
                               lambda0.t.hat, S0t, obstimes, obstimes0, y, exp.X.weighted,                      
                               number.times, Tau1Tau2.times) 
  rm(S0t, S1t, obstimes, obstimes0, y, exp.X.weighted, Tau1Tau2.times)
  gc()

  list(A.RH=resid, A.CumBH=infl.Lambda0.Tau1Tau2)
}

## -----------------------------------------------------------------------------
## Function: auxiliary.construction()
## -----------------------------------------------------------------------------
## Description: This function returns the auxiliary variables proposed by 
##              Breslow et al. (Stat. Biosci., 2009) and Breslow et al. (IMS,
##              2013). They are based on the relative hazard and cumulative 
##              baseline hazard influences, when using the cohort data with 
##              imputed covariates values
## -----------------------------------------------------------------------------
auxiliary.construction <- function (mod, Tau1=NULL, Tau2=NULL, method="Breslow",
                                    time.on.study=NULL, casecohort=NULL) {
  
  if (!(method %in% c("Breslow", "Breslow2013", "Shin"))) stop("ERROR: method must be 'Breslow', 'Breslow2013' or 'Shin'")
  if (method == "Shin") { 
    if (!is.data.frame(casecohort)) stop("ERROR: casecohort must be a data frame for method = 'Shin'")
    if (is.null(time.on.study)) stop("ERROR: time.on.study must be specified for method = 'Shin'")
  } 

  infl.Lambda0.Tau1Tau2 <- NULL
  A.RH.Shin             <- NULL
  A.Shin                <- NULL

  # ----------------------------------------------------------------------------
  # Computation of the influences for log-relative hazard, beta ----------------
  infl.beta <- residuals(mod, type = "dfbeta", weighted = T) 

  if (method == "Breslow2013") {
    # ----------------------------------------------------------------------------
    # Quantities needed for the influences ---------------------------------------
    mod.detail        <- coxph.detail(mod, riskmat = T)
    riskmat           <- mod.detail$riskmat 
    number.times      <- ncol(riskmat) 
    n                 <- nrow(riskmat) 
    X                 <- model.matrix(mod) 
    p                 <- ncol(X) 
    weights           <- mod$weights # weights used for the fit, omega_i
    if (is.null(weights)) {
      weights         <- 1
    } # when using the whole cohort, omega_i = 1. Should be the case here
  
    beta.hat          <- mod$coefficients
  
    exp.X.weighted    <- weights * exp(X %*% beta.hat)
    Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                          ncol = number.times, byrow = FALSE) 
    S0t               <- t(riskmat) %*% (exp.X.weighted) 
    S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                  ncol = p, byrow = FALSE)) 
    observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
    dNt               <- matrix(0, n, number.times) 
    dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
    dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
    lambda0.t.hat     <- t(colSums(dNt) / S0t)
    if (is.null(Tau1)) {
      Tau1            <- floor(min(mod.detail$time))
    }
    if (is.null(Tau2)) {
      Tau2            <- floor(max(mod.detail$time))
    }
    Tau1Tau2.times  <- which((Tau1 < mod.detail$time) & (mod.detail$time <= Tau2)) 
  
    # ----------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] ---------------------------------------------------------------
  
    infl.lambda0.t    <- (dNt - (Y.exp.X.weighted + infl.beta %*% t(S1t)) * 
                          matrix(lambda0.t.hat, nrow = n, ncol = number.times,
                                byrow = TRUE)) /
                        matrix(S0t, nrow = n, ncol = number.times,
                                byrow = TRUE) 
    infl.Lambda0.Tau1Tau2 <- rowSums(infl.lambda0.t[, Tau1Tau2.times])
  }

  if (tolower(method) == "shin") {
    indiv.phase2 <- rownames(casecohort)
    A            <- cbind(1, infl.beta)
    total        <- colSums(A)
    tmp          <- calibration(A.phase2 = A[indiv.phase2, , drop=FALSE], 
                         design.weights = casecohort[, "weights", drop=TRUE], 
                         total = total, eta0 = rep(0, ncol(A)))
    form         <- mod$formula
    weights.calib1 <- as.vector(tmp$calibrated.weights)
    id <- NULL
    casecohort$weights.calib1 <- weights.calib1
    mod.calib1 <- coxph(form, data=casecohort, 
                        weights=weights.calib1, robust=TRUE, id=id)
    coef <- mod.calib1$coefficients
    cnms <- names(coef)
    mat  <- model.matrix(mod)    
    nms  <- intersect(cnms, colnames(mat))
    if (!length(nms)) stop("ERROR: no common names with coefficients and model.matrix") 
    coef      <- coef[nms]
    mat       <- mat[, nms, drop=FALSE]  
    coef      <- matrix(mod.calib1$coefficients, ncol=1)
    A.Shin    <- time.on.study*exp(mat %*% coef)
    A.RH.Shin <- infl.beta
  }

  list(A.RH.Breslow = infl.beta, A.CumBH.Breslow = infl.Lambda0.Tau1Tau2,
       A.RH.Shin = A.RH.Shin, A.PR.Shin = A.Shin)
}

