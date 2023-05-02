caseCohortCoxSurvival <- function(data, status, time, 
                               cox.phase1=NULL, cox.phase2=NULL, other.covars=NULL, 
                               strata=NULL, weights.phase2=NULL, 
                               calibrated=FALSE, 
                               subcohort=NULL, 
                               subcohort.strata.counts=NULL,
                               predict=TRUE, predicted.cox.phase2=NULL,
                               predictors.cox.phase2=NULL,
                               aux.vars=NULL, aux.method="Shin",
                               phase3=NULL, strata.phase3=NULL,
                               weights.phase3=NULL, weights.phase3.type="both",
                               Tau1=NULL, Tau2=NULL, x=NULL,
                               weights.op=NULL, print=1) {

  check_data(data)
  check_status(status, data)
  check_time_vars(data, time)
  tmp <- colnames(data)
  check_covars(cox.phase1, data, name="cox.phase1")
  check_covars(cox.phase2, data, name="cox.phase2")
  if (!length(cox.phase1) && !length(cox.phase2)) stop("ERROR: cox.phase1 and/or cox.phase2 must be specified")
  check_covars(other.covars, data, name="other.covars")
  check_strata(strata, data)
  check_binary(calibrated, "calibrated")
  check_binary(predict, "predict")
  check_predicted.covars(predicted.cox.phase2, data)
  check_proxy.vars(predictors.cox.phase2, data) 
  check_auxvars(aux.vars, data)
  aux.method <- check_aux.method(aux.method)
  check_phase23(subcohort, phase3, data, status) 
  check_weights(weights.phase2, data, name="weights.phase2")
  check_weights(weights.phase3, data, name="weights.phase3")
  if (!is.null(weights.phase3.type)) {
    weights.phase3.type <- check_opString(weights.phase3.type, "weights.phase3.type", 
                                          c("both", "design", "estimated"))
  }
  check_strata(strata.phase3, data)
  if (!is.null(Tau1)) check_numeric(Tau1, "Tau1", pos=0)
  if (!is.null(Tau2)) check_numeric(Tau2, "Tau2")
  x <- check_risk.data(x, c(cox.phase1, cox.phase2), data) 
  check_calibrated.op(weights.op) 
  check_integer(print, "print", valid=0:3) 
  check_phase2.strata.counts(subcohort.strata.counts, data, strata)
 
  objlist <- list(status=status, time=time, 
                  cox.phase1=cox.phase1, cox.phase2=cox.phase2, other.covars=other.covars,
                  strata=strata, subcohort=subcohort, weights.phase2=weights.phase2,
                  predict=predict, predicted.cox.phase2=predicted.cox.phase2,
                  predictors.cox.phase2=predictors.cox.phase2,
                  calibrated=calibrated, 
                  aux.vars=aux.vars, aux.method=aux.method, 
                  phase3=phase3, weights.phase3.type=weights.phase3.type,
                  strata.phase3=strata.phase3, weights.phase3=weights.phase3,
                  tau1=Tau1, tau2=Tau2, phase2.strata.counts=subcohort.strata.counts,
                  weights.op=weights.op, print=print, min.nrow=5, DEBUG=0)
  ret <- cc.main(data, x, objlist)
  ret
}

cc.main <- function(data, risk.data, obj) {

  # Set up the data
  tmp  <- setUpData(data, obj)
  data <- tmp$data
  rn0  <- tmp$rownames0
  obj  <- tmp$obj.list
  rm(tmp)
  gc()

  # Get auxiliary data if needed
  A <- getAuxData(data, obj)

  # Get calibrated weights if needed, will be stored in variable obj$weights.calibrated
  if (obj$calibrated) {
    tmp  <- getCalibratedWgts(data, obj, A)
    data <- tmp$data
    obj  <- tmp$obj.list
    rm(tmp); gc()
  }

  # Get the weights for coxph
  weights <- getCoxph_weights(data, obj, A)
  subset  <- getSubsetVec(data, obj) # All TRUE for whole cohort analysis

  # Fit the model (case-cohort or whole cohort)
  tmp        <- fitCoxModel1(data, subset, weights, obj) 
  fit        <- tmp$fit
  obj$covars <- tmp$covars
  rm(tmp, subset); gc()

  # Influences
  infl <- getInfluences(data, obj, A, fit)
  rm(weights, A); gc()

  if (obj$print) cat("Computing robust variance estimates for parameters\n")
  # Robust variance for estimates
  rv.est <- getRobustVarForEst(infl) 

  # Variance for estimates
  var.est <- getVarianceForEst(data, obj, infl, fit)

  # Variance for risk data
  var.data <- try(getVarianceForData(risk.data, data, obj, infl), silent=TRUE)
  if ("try-error" %in% class(var.data)) {
    if (obj$print) {
      print(var.data)
      cat("Error computing risk.data estimates, see above error message\n") 
    }
    var.data <- NULL
  }

  # Set return object
  ret <- setReturnObject(data, obj, infl, var.est, rv.est, var.data, fit)
  if (obj$print) cat("Finished\n")

  ret
}

getInputArgNames <- function() {

  ret <- c("status", "time", 
           "cox.phase1", "cox.phase2", "other.covars", 
           "strata", "weights.phase2", "calibrated", 
           "subcohort", "subcohort.strata.counts",
           "predict", "predicted.cox.phase2", "predictors.cox.phase2",
           "aux.vars", "aux.method",
           "phase3", "strata.phase3", "weights.phase3", "weights.phase3.type",
           "weights.op", "print")
  ret
}

getArgsFromObj <- function(obj) {

  args <- getInputArgNames()
  ret  <- list(Tau1=obj$tau1, Tau2=obj$tau2)
  for (arg in args) {
    val <- obj[[arg, exact=TRUE]]
    if (!length(val)) {
      ret[arg]   <- list(NULL)
    } else {
      ret[[arg]] <- val 
    }
  }
  ret
}

getObjFromInfl <- function(infl, name) {

  ret      <- NULL
  infl.nms <- names(infl)
  if (!length(infl.nms)) return(ret)

  for (nm in infl.nms) {
    lst <- infl[[nm, exact=TRUE]]
    if (!length(lst)) next
    ret <- lst[[name, exact=TRUE]]
    if (length(ret)) break
  }
  ret
}

setReturnObject <- function(data, obj, infl, var.est, rv.est, var.data, fit) {

  ret      <- list()
  dat      <- NULL
  vars     <- c(obj$status, obj$weights, obj$covars.orig, obj[["phase2", exact=TRUE]],
                obj[["strata", exact=TRUE]], obj[["weights.phase2", exact=TRUE]],
                obj[["phase3", exact=TRUE]], obj[["strata.phase3", exact=TRUE]])
  tmp      <- vars %in% colnames(data)
  vars     <- vars[tmp]
  if (length(vars)) dat <- data[, vars, drop=FALSE]
  chtimes <- obj[["changed.times", exact=TRUE]]
  obj$changed.times <- NULL
  risk.obj <- list(data=dat, infl=infl, obj.list=obj)
  rm(dat); gc()

  beta   <- getObjFromInfl(infl, "beta.hat")
  lambda <- getObjFromInfl(infl, "Lambda0.Tau1Tau2.hat")
  ret[["beta"]]   <- beta
  ret[["Lambda0"]] <- lambda

  tmp <- var.est[["var", exact=TRUE]]
  if (!is.null(tmp)) {
    ret[["beta.var"]]   <- tmp$var.beta
    ret[["Lambda0.var"]] <- tmp$var.lam
  }
  tmp <- var.est[["var.T", exact=TRUE]]
  if (!is.null(tmp)) {
    ret[["beta.var.estimated"]]   <- tmp$var.beta
    ret[["Lambda0.var.estimated"]] <- tmp$var.lam
  }
  tmp <- var.est[["var.F", exact=TRUE]]
  if (!is.null(tmp)) {
    ret[["beta.var.design"]]   <- tmp$var.beta
    ret[["Lambda0.var.design"]] <- tmp$var.lam
  }

  tmp <- rv.est[["infl", exact=TRUE]]
  if (!is.null(tmp)) {
    ret[["beta.robustvar"]]   <- tmp$infl.beta
    ret[["Lambda0.robustvar"]] <- tmp$infl.Lambda0.Tau1Tau2
  }
  tmp <- rv.est[["infl.T", exact=TRUE]]
  if (!is.null(tmp)) {
    ret[["beta.robustvar.estimated"]]   <- tmp$infl.beta
    ret[["Lambda0.robustvar.estimated"]] <- tmp$infl.Lambda0.Tau1Tau2
  }
  tmp <- rv.est[["infl.F", exact=TRUE]]
  if (!is.null(tmp)) {
    ret[["beta.robustvar.design"]]   <- tmp$infl.beta
    ret[["Lambda0.robustvar.design"]] <- tmp$infl.Lambda0.Tau1Tau2
  }

  tmp <- var.data[["rv", exact=TRUE]]
  if (!is.null(tmp)) ret[["Pi.var"]] <- tmp
  tmp <- var.data[["rv.T", exact=TRUE]]
  if (!is.null(tmp)) ret[["Pi.var.estimated"]] <- tmp
  tmp <- var.data[["rv.F", exact=TRUE]]
  if (!is.null(tmp)) ret[["Pi.var.design"]] <- tmp

  ret[["coxph.fit"]] <- fit
  if (length(chtimes)) ret[["changed.times"]] <- chtimes
  ret[["args"]]      <- getArgsFromObj(obj)
  ret[["risk.obj"]]  <- risk.obj
  
  class(ret) <- "casecohortcoxsurv"
  ret
}

fitCoxModel1 <- function(data, subset, weights, obj) {

  status <- obj$status
  time   <- obj$time
  covars <- obj$covars
  id.var <- obj$id.var
  print  <- obj$print
  cntl   <- obj$coxph.control

  n    <- nrow(data)
  if (is.null(subset)) subset <- rep(TRUE, n)
  form   <- as.formula(getCoxFormula(status, time, covars))
  wcflag <- all(subset) # whole cohort analysis

  # For whole-cohort analysis, do not use weights
  # Weights could be the calibrated weights from calibration(). In this case,
  #  they are for the phase2 sample
  if (!wcflag) { 
    if (length(weights) == nrow(data)) weights <- weights[subset]
    weights <- as.vector(weights)
  } else {
    weights <- NULL
  }
  ids <- data[subset, id.var, drop=TRUE]
  
  # Fit the model, check for warning
  op0 <- options()
  on.exit(options(op0))
  options(warn=2)
  fit <- try(coxph(form, data=data[subset, , drop=FALSE], robust=TRUE,
               weights=weights, id=ids, control=cntl))
  options(op0)
  if ("try-error" %in% class(fit)) stop("ERROR fitting Cox model")
  if (print > 1) {
      cat("*********************************************\n") 
    if (!wcflag) {
      cat("Fitted Cox model using the case-cohort data:\n")
    } else {
      cat("Fitted Cox model using the whole cohort data:\n")
    }
      cat("*********************************************\n") 
    print(fit)
    cat("\n")
  }
  beta <- fit$coefficients
  tmp  <- is.na(beta)
  if (any(tmp)) {
    nms    <- names(beta)[tmp]
    tmp    <- !(covars %in% nms)
    covars <- covars[tmp] 
    if (!length(covars)) stop("ERROR fitting model, check covars")
    if (print) {
      cat("NOTE: not all parameters were estimated in the Cox model.\n")
      if (print > 1) cat("NOTE: parameters not estimated in the above model will be dropped.\n")
    }
    form   <- as.formula(getCoxFormula(status, time, covars))
    fit    <- coxph(form, data=data[subset, , drop=FALSE], robust=TRUE,
               weights=weights, id=ids, control=cntl)
    beta   <- fit$coefficients
  } 
  if (any(is.na(beta))) stop("ERROR fitting Cox model, check covars")
  list(fit=fit, covars=covars)
}

getSubsetVec <- function(data, obj) {

  subsetv <- obj[["phase3", exact=TRUE]]
  if (is.null(subsetv)) subsetv <- obj[["phase2", exact=TRUE]]
  if (is.null(subsetv)) {
    ret <- rep(TRUE, nrow(data))
  } else {
    ret <- data[, subsetv, drop=TRUE]
  }
  ret
}

getStrataVec <- function(data, obj) {

  v <- obj[["strata", exact=TRUE]]
  if (is.null(v)) {
    ret <- NULL
  } else {
    ret <- data[, v, drop=TRUE]
  }
  ret
}

oneMinusPiExpx <- function(x, beta.hat, Lambda0.Tau1Tau2.hat) {

  if (is.null(dim(x))) dim(x) <- c(1, length(x))
  if (is.null(dim(beta.hat))) dim(beta.hat) <- c(length(beta.hat), 1)

  exp.x    <- as.vector(exp(x %*% beta.hat))
  Pi.x.hat <- 1 - exp(-exp.x*Lambda0.Tau1Tau2.hat)
  ret      <- (1 - Pi.x.hat)*exp.x 
  list(ret=ret, Pi.x.hat=Pi.x.hat)
}

getStratCounts <- function(data, obj) {

  stratVec <- NULL
  if (obj$stratFlag) stratVec <- data[, obj$strata, drop=TRUE] 
  casecohort <- data[, obj$subcohort, drop=TRUE]

  if (is.null(stratVec)) {
    ret <- list(cnts.cohort=length(casecohort), cnts.casecohort=sum(casecohort))
    return(ret)
  }
  p2counts        <- obj[["phase2.strata.counts", exact=TRUE]]
  p2Flag          <- !is.null(p2counts)
  ustrat          <- obj$strata.new
  ustrat.orig     <- obj$strata.orig 
  N               <- length(ustrat)
  cnts.cohort     <- rep(0, N)
  cnts.casecohort <- rep(0, N)

  for (i in 1:N) {
    val            <- ustrat[i]
    val.orig       <- ustrat.orig[i]
    tmp            <- stratVec %in% val
    cnts.cohort[i] <- sum(tmp)
    if (p2Flag) {
      cnt <- p2counts[[as.character(val.orig), exact=TRUE]]
      if (is.null(cnt)) stop("INTERNAL CODING ERROR")
      cnts.casecohort[i] <- cnt
    } else {
      tmp2 <- tmp & casecohort 
      cnt  <- sum(tmp2)
      if (cnt < 2) {
        msg <- paste0("There is ", cnt, " subcohort individual(s) in stratum ", val.orig)
        warning(msg)
      }
      cnts.casecohort[i] <- cnt
    }
  }

  list(cnts.cohort=cnts.cohort, cnts.casecohort=cnts.casecohort)
}

getPhase2Weights <- function(data, obj) {

  # strata.counts have been set at this point, and the strata are 0, 1, ...
  scounts <- obj[["strata.counts", exact=TRUE]]
  if (!length(scounts)) stop("INTERNAL CODING ERROR")
  cohort.counts    <- scounts$cnts.cohort
  subcohort.counts <- scounts$cnts.casecohort

  n       <- nrow(data)
  ret     <- rep(NA, n)
  strata  <- getStrataVec(data, obj)
  if (is.null(strata)) strata <- rep(0, n)
  N       <- obj$n.strata
  for (strat in 0:(N-1)) {
    tmp      <- strata %in% strat
    if (!any(tmp)) stop("INTERNAL CODING ERROR 1")
    n        <- cohort.counts[strat+1]
    m        <- subcohort.counts[strat+1]
    ret[tmp] <- n/m
  }
  # Cases are set to 1
  tmp <- data[, obj$status, drop=TRUE] %in% 1
  ret[tmp] <- 1

  if (any(is.na(ret))) stop("INTERNAL CODING ERROR 2")
  ret
}

cc_addDummyVars <- function(data, var) {

  form <- as.formula(paste0("~ ", var))
  mat  <- model.matrix(form, data=data)
  mat  <- mat[, -1, drop=FALSE]
  cx   <- colnames(data)
  cm   <- colnames(mat)
  rm   <- rownames(mat)
  if (!any(cm %in% cx)) {
    data[, cm]   <- NA
    data[rm, cm] <- mat
  } else {
    while (1) {
      cm <- paste0(".", cm)
      if (!any(cm %in% cx)) break
    }
    data[, cm]   <- NA
    data[rm, cm] <- mat
  }
  list(data=data, dvars=cm)
}

cc_addAllDummyVars <- function(data, vars, obj) {

  vars    <- unique(vars)
  nvars   <- length(vars)
  if (!nvars) return(list(data=data, obj.list=obj))
  vartype   <- obj$vartype
  catvars   <- NULL
  dummyvars <- list() 
  varlevels <- list()
  for (i in 1:nvars) {
    var  <- vars[i]
    type <- vartype[[var, exact=TRUE]]
    if (type %in% "cat") {
      data[, var]      <- as.factor(data[, var, drop=TRUE])
      tmp              <- cc_addDummyVars(data, var)
      data             <- tmp$data
      dvars            <- tmp$dvars
      catvars          <- c(catvars, var)
      dummyvars[[var]] <- dvars
      varlevels[[var]] <- levels(data[, var]) 
    }
  } 
  obj$catvars   <- catvars
  obj$dummyvars <- dummyvars
  obj$varlevels <- varlevels

  list(data=data, obj.list=obj)

}

cc_getAllVars <- function(obj, cox.phase2=1, time=1, weights.phase3=1) {

  ret <- c(obj$status, 
           obj[["cox.phase1", exact=TRUE]],
           obj[["other.covars", exact=TRUE]],
           obj[["strata", exact=TRUE]],
           obj[["weights.phase2", exact=TRUE]],
           obj[["subcohort", exact=TRUE]],
           unlist(obj[["predicted.cox.phase2", exact=TRUE]]),
           unlist(obj[["predictors.cox.phase2", exact=TRUE]]),
           obj[["aux.vars", exact=TRUE]],
           obj[["phase2", exact=TRUE]],
           obj[["phase3", exact=TRUE]],
           obj[["strata.phase3", exact=TRUE]])
  if (cox.phase2) ret <- c(ret, obj[["cox.phase2", exact=TRUE]])
  if (time) ret <- c(ret, obj$time)
  if (weights.phase3) ret <- c(ret, obj[["weights.phase3", exact=TRUE]])

  ret <- unique(ret)
  ret

}

cc_checkForMissing <- function(data, obj) {

  miss.list <- obj$varmissing
  vars      <- cc_getAllVars(obj, cox.phase2=0, time=0, weights.phase3=0)
  for (v in vars) {
    miss <- miss.list[[v, exact=TRUE]]
    if (is.null(miss)) next # Some variables were defined later (phase2, phase3, etc)
    if (miss) {
      msg <- paste0("ERROR: variable ", v, " cannot contain missing/non-finite values")
      stop(msg)
    } 
  }

  # time must be defined on the case-cohort
  subset <- getSubsetVec(data, obj)
  vars   <- obj$time
  for (v in vars) {
    vec  <- data[subset, v, drop=TRUE]
    miss <- any(!is.finite(vec))
    if (miss) {
      msg <- paste0("ERROR: variable ", v, " cannot contain missing/non-finite values",
                    " for individuals in the case-cohort")
      stop(msg)
    } 
  }

  # weights.phase3 must be defined on phase2
  v    <- obj[["weights.phase3", exact=TRUE]]
  subv <- obj[["phase2", exact=TRUE]]
  if (length(v) && length(subv)) {
    subset <- data[, subv, drop=TRUE] %in% 1
    vec    <- data[subset, v, drop=TRUE]
    miss   <- any(!is.finite(vec))
    if (miss) {
      msg <- paste0("ERROR: variable ", v, " cannot contain missing/non-finite values",
                    " for individuals in phase-two")
      stop(msg)
    } 
  }

  NULL
}

cc_updateVarsWithDummy <- function(vars, obj) {

  if (!length(vars)) return(NULL)
  catvars <- obj[["catvars", exact=TRUE]]
  if (!length(catvars)) return(vars)

  tmp   <- vars %in% catvars 
  if (!any(tmp)) return(vars)

  cvars <- vars[tmp]
  dlist <- obj$dummyvars
  if (!length(dlist)) stop("INTERNAL CODING ERROR")
  for (cvar in cvars) {
    id    <- match(cvar, vars)
    dvars <- dlist[[cvar, exact=TRUE]]
    if (!length(dvars)) stop("INTERNAL CODING ERROR 2")
    n <- length(vars)
    if (id == 1) {
      vars <- c(dvars, vars[-1])
    } else if (id == n) {
      vars <- c(vars[-n], dvars)
    } else {
      vars <- c(vars[1:(id-1)], dvars, vars[(id+1):n])
    }
  }

  vars
}

checkRiskDataType <- function(risk.data, vars, data.vartype) {

  if (!length(risk.data)) return(NULL)
  vartype <- getVarTypes(risk.data, vars, only.finite=1)
  for (var in vars) {
    risk.type <- vartype[[var, exact=TRUE]]
    data.type <- data.vartype[[var, exact=TRUE]]
    if (!length(risk.type) || !length(data.type)) stop("INTERNAL CODING ERROR")
    tmp1 <- (risk.type == "cat") && (data.type != "cat")
    tmp2 <- (risk.type != "cat") && (data.type == "cat")
    if (tmp1 || tmp2) {
      msg <- paste0("ERROR: variable ", var, " has a different type in data and in risk.data")
      stop(msg)
    }
  }
  vartype
}

setStrataInts <- function(vec) {

  vec  <- unfactor(vec)
  ret  <- vec
  uvec <- sort(unique(vec)) 
  nvec <- length(uvec)
  for (i in 1:nvec) {
    tmp      <- vec %in% uvec[i]
    ret[tmp] <- i - 1
  }
  list(vec=ret, strata.orig=uvec, strata.new=0:(nvec-1))
}

setStrataVars <- function(data, obj) {

  var <- obj[["strata", exact=TRUE]]
  if (!is.null(var)) {
    tmp         <- setStrataInts(data[, var, drop=TRUE]) 
    data[, var] <- tmp$vec
    # Need orig strata
    obj$strata.orig <- tmp$strata.orig
    obj$strata.new  <- tmp$strata.new
  }
  var <- obj[["strata.phase3", exact=TRUE]]
  if (!is.null(var)) data[, var] <- setStrataInts(data[, var, drop=TRUE])$vec 
  
  list(data=data, obj=obj)
}

predictedCovarsCheck <- function(data, obj) {

  # A predicted column must be specified for each phase 2 covariate, when predict=FALSE 
  #   and weight calibration
 
  p2vars <- obj[["cox.phase2", exact=TRUE]]
  if (is.null(p2vars)) return(NULL)
    
  predvars <- obj[["predicted.cox.phase2", exact=TRUE]]
  if (is.null(predvars)) return(NULL)
  
  # If a list then check that it must contain all phase 2 covars
  listFlag <- isList(predvars)
  if (!listFlag) stop("INTERNAL CODING ERROR 1")
  
  nms <- names(predvars)
  tmp <- !(p2vars %in% nms)
  if (any(tmp)) {
    str <- getQuotedVarStr(p2vars[tmp])
    msg <- paste0("ERROR: predicted.cox.phase2 is missing predictions for ", str)
    stop(msg) 
  } 

  # Predicted columns cannot have missing values
  vartype <- obj$vartype
  nvals   <- c("cont", "binary")
  for (i in 1:length(p2vars)) {
    p2v   <- p2vars[i]
    predv <- predvars[[p2v, exact=TRUE]]
    if (length(predv) != 1) stop("INTERNAL CODING ERROR 2")

    # Variable must be of the same type as the phase 2 covariate
    p2type   <- vartype[[p2v, exact=TRUE]]
    predtype <- vartype[[p2v, exact=TRUE]]
    tmp      <- (p2type == "cat") && (predtype != "cat")
    if (tmp) {
      msg <- paste0("ERROR: the phase-two covariate ", p2v, " is categorical, but \n",
                    "predicted.cox.phase2 = ", predv, " is not categorical\n")
      stop(msg)
    }
    tmp      <- (p2type != "cat") && (predtype == "cat")
    if (tmp) {
      msg <- paste0("ERROR: the phase-two covariate ", p2v, " is not categorical, but \n",
                    "predicted.cox.phase2 = ", predv, " is categorical\n")
      stop(msg)
    }
  }
  
  NULL
}

removeMiss_phase2 <- function(data, obj) {

  if (!obj$calibrated) return(data)
  if (!obj$phase2Flag) return(data)

  keep   <- rep(TRUE, nrow(data))
  p2     <- data[, obj$phase2, drop=TRUE]
  wgtv   <- obj[["calibrated.weights", exact=TRUE]]
  avars  <- obj[["aux.vars", exact=TRUE]]
  covars <- obj[["covars", exact=TRUE]]
  vars   <- unique(c(wgtv, avars, covars))
  p2vec  <- p2 %in% 1
  if (length(vars)) {
    # covars and aux.vars could contain categorical vars, so loop over each one
    for (v in vars) {
      vec <- data[, v, drop=TRUE]
      if (is.numeric(vec)) {
        keep <- keep & is.finite(vec)
      } else {
        keep <- keep & !is.na(vec)
      }
    }
  }
  if (!all(keep)) data <- data[keep, , drop=FALSE]
  data
}

tau12.check <- function(data, obj) {

  print <- obj$print
  tmp   <- data[, obj$status, drop=TRUE] %in% 1
  timev <- obj$time
  if (length(timev) > 1) timev <- timev[2]
  time  <- data[tmp, timev, drop=TRUE]
  if (length(time) < 2) stop("ERROR: too few event times")
  Tau1  <- obj[["tau1", exact=TRUE]]
  Tau2  <- obj[["tau2", exact=TRUE]]
  if (is.null(Tau1)) {
    Tau1 <- floor(min(time, na.rm=TRUE))
    if (print) cat(paste0("NOTE: setting tau1 to ", Tau1, " (minimum event time)\n"))
  }
  if (is.null(Tau2)) {
    Tau2 <- floor(max(time, na.rm=TRUE))
    if (print) cat(paste0("NOTE: setting tau2 to ", Tau2, " (maximum event time)\n"))
  }
  if (Tau1 >= Tau2) stop("ERROR: tau1 >= tau2")
  tmp <- (Tau1 < time) & (time <= Tau2)
  if (!any(tmp)) stop("ERROR: no event times in (tau1, tau2], adjust tau1 and/or tau2")
 
  obj$tau1 <- Tau1
  obj$tau2 <- Tau2

  obj

}

checkB.phase3Cols <- function(B) {

  ret <- NULL
  for (i in 1:ncol(B)) {
    tmp <- B[, i, drop=TRUE] == 0
    tmp[is.na(tmp)] <- TRUE
    if (all(tmp)) ret <- c(ret, i)
  }
  ret
}

computePhase3Weights <- function(data, obj) {

  if (!obj$phase3Flag) return(NULL)
  est.op   <- obj$weights.op
  B.phase2 <- getB.phase2(data, obj) 
  tmp      <- data[, obj$phase3, drop=TRUE]  
  p3subs   <- rownames(data[tmp, , drop=FALSE])
  B.phase3 <- B.phase2[p3subs, , drop=FALSE]
  rem      <- checkB.phase3Cols(B.phase3)
  if (length(rem)) {
    warning("Removing columns in B.phase3 to compute phase-three weights")
    B.phase2 <- B.phase2[, -rem, drop=FALSE]
    B.phase3 <- B.phase3[, -rem, drop=FALSE]
  }
  ret      <- calibration(B.phase3, 1, colSums(B.phase2), eta0=NULL, 
                         niter.max=est.op[["niter.max", exact=TRUE]], 
                         epsilon.stop=est.op[["epsilon.stop", exact=TRUE]])$calibrated.weights
  ret
}

getCoxph_weights <- function(data, obj, A) {

  if (obj$phase3Flag) {
    wgt2v <- obj[["weights.phase2", exact=TRUE]]
    if (is.null(wgt2v)) stop("INTERNAL CODING ERROR 1")
    wgt3v <- obj[["weights.phase3", exact=TRUE]]
    if (is.null(wgt3v)) stop("INTERNAL CODING ERROR 2")
    ret <- data[, wgt3v, drop=TRUE]*data[, wgt2v, drop=TRUE] 
  } else if (obj$phase2Flag) {
    if (obj$calibrated) {
      wgtv <- obj[["weights.calibrated", exact=TRUE]]
      if (is.null(wgtv)) stop("INTERNAL CODING ERROR 3")
      ret  <- data[, wgtv, drop=TRUE]
    } else {
      ret <- data[, obj$weights.phase2, drop=TRUE]
    }
  } else {
    ret <- rep(1, nrow(data))
  }
  
  ret 
}

cc_removeMiss <- function(data, obj) {

  # Remove missing values
  n0   <- nrow(data)
  vars <- c(obj$status, obj$time, obj$strata, obj[["weights", exact=TRUE]])
  if (obj$phase2Flag) vars <- c(vars, obj[["phase2", exact=TRUE]])
  if (obj$phase3Flag) vars <- c(vars, obj[["phase3", exact=TRUE]],
                                  obj[["strata.phase3", exact=TRUE]])
  if (!obj$phase2Flag) vars <- c(vars, obj[["covars", exact=TRUE]])       
  if (obj$calibrated) {
    auxvars <- obj[["aux.vars", exact=TRUE]]
    if (length(auxvars)) {
      vars <- c(vars, auxvars)
    } else {
      pvars <- unlist(obj[["proxy.vars", exact=TRUE]])
      if (is.null(pvars)) pvars <- obj$covars
      vars <- c(vars, pvars)
    }
  } 
  if (obj$phase3Flag) vars <- c(vars, obj$covars)
  if (!obj$phase3Flag && !obj$calibrated) vars <- c(vars, obj$covars)
  vars <- unique(vars)   
  data <- removeMiss(data, NULL, vars, NULL, min.nrow=obj$min.nrow, print=obj$print)$data

  # Remove missing from phase2
  data <- removeMiss_phase2(data, obj)
  n1   <- nrow(data)
  if (obj$print) cat(paste0(n0-n1, " rows have been removed due to missing values\n"))

  data
}

setup_checkForNoMiss <- function(data, obj, obj.name) {

  obj.vars <- obj[[obj.name, exact=TRUE]]
  vars     <- unlist(obj.vars)
  n        <- length(vars)
  if (!n) return(obj)

  miss    <- obj$varmissing
  ok      <- rep(TRUE, n)
  prnt    <- obj$print
  for (i in 1:n) {
    var  <- vars[i]
    flag <- miss[[var, exact=TRUE]]
    if (!length(flag)) stop("INTERNAL CODING ERROR 0")
    if (flag) {
      ok[i] <- FALSE
      if (prnt) {
        v   <- getQuotedVarStr(var)
        msg <- paste0("NOTE: ", obj.name, " variable ", v, 
                      " contains missing values in the cohort, ",
                      "removing ", v, " from ", obj.name, ".\n")
        cat(msg) 
      }
    } 
  }
  remFlag <- !all(ok)
  vars    <- vars[ok]
  if (!length(vars)) {
    obj[[obj.name]] <- NULL
    if (prnt) cat(paste0("NOTE: all variables have been removed from ", obj.name, " other.covars.\n"))
  } else if (remFlag) {
    obj[[obj.name]] <- keepValsInVecList(obj.vars, vars)
  }

  obj
}

setup_check_covars_subcohort <- function(data, obj) {

  # check the 3 sets of covariates
  print <- obj$print
  #if (print) cat("Checking variables ...\n")
  varmissing <- obj$varmissing

  # Other covariates - no missing in cohort
  obj <- setup_checkForNoMiss(data, obj, "other.covars")

  # predicted variables - no missing in cohort
  obj <- setup_checkForNoMiss(data, obj, "predicted.cox.phase2")

  # predictor variables - no missing in cohort
  obj <- setup_checkForNoMiss(data, obj, "predictors.cox.phase2")

  # auxiliary variables - no missing in cohort
  obj <- setup_checkForNoMiss(data, obj, "aux.vars")

  covars1 <- obj[["cox.phase1", exact=TRUE]]
  covars2 <- obj[["cox.phase2", exact=TRUE]]
  n1      <- length(covars1)
  n2      <- length(covars2)
  if (!n1 && !n2) stop("ERROR: specify cox.phase1 and/or cox.phase2")
  phase2Flag <- !is.null(obj[["subcohort", exact=TRUE]])
  #if (n2 && !phase2Flag) stop("ERROR: subcohort=NULL but cox.phase2 != NULL")

  both <- intersect(covars1, covars2)
  if (length(both)) {
    str <- getQuotedVarStr(both)
    msg <- paste0("NOTE: cox.phase1 and cox.phase2 contain overlapping variable(s) ", str,
                  ", they will be moved to either cox.phase1 or cox.phase2, but not both.\n")
    cat(msg) 
  }

  # Phase1 - no missing values
  ret1 <- NULL
  ret2 <- NULL
  if (n1) {
    for (i in 1:n1) {
      var  <- covars1[i]
      flag <- varmissing[[var, exact=TRUE]]
      if (!length(flag)) stop("INTERNAL CODING ERROR 1")
      if (flag) {
        ret2 <- c(ret2, var)
        if (print) {
          v   <- getQuotedVarStr(var)
          msg <- paste0("NOTE: cox.phase1 variable ", v, " contains missing values in the cohort, ",
                        "moving ", v, " to cox.phase2.\n")
          cat(msg) 
        }
      } else {
        # OK
        ret1 <- c(ret1, var)
      }
    }
  }

  obj$cox.phase1 <- ret1
  obj$cox.phase2 <- ret2
  #if (!phase2Flag) return(obj)

  # Phase 2 - check for missing in subcohort
  if (n2) {
    for (i in 1:n2) {
      var  <- covars2[i]
      flag <- varmissing[[var, exact=TRUE]]
      if (!length(flag)) stop("INTERNAL CODING ERROR 2")
      if (!flag) {
        ret1 <- c(ret1, var)
        if (print) {
          v   <- getQuotedVarStr(var)
          msg <- paste0("NOTE: cox.phase2 variable ", v, " contains no missing values in the cohort, ",
                        "moving ", v, " to cox.phase1.\n")
          cat(msg) 
        }
      } else {
        # OK
        ret2 <- c(ret2, var)
      }
    }
  }
  obj$cox.phase1 <- ret1
  obj$cox.phase2 <- ret2 
  if (!length(obj[["cox.phase2", exact=TRUE]])) {
    if (n2 && print) cat("NOTE: no phase-two covariates left.\n")
    obj <- setup_setArgsMsg(obj, "predict", FALSE, str1="since there are no phase-two covariates.")
    obj <- setup_setArgsMsg(obj, "predicted.cox.phase2", NULL, str1="since there are no phase-two covariates.")
    obj <- setup_setArgsMsg(obj, "predictors.cox.phase2", NULL, str1="since there are no phase-two covariates.")
  }

  # Check that predictors are in other.covars or cox.phase1
  pvars <- unlist(obj[["predictors.cox.phase2", exact=TRUE]])
  if (length(pvars)) {
    ovars <- obj[["other.covars", exact=TRUE]]
    tmp   <- !(pvars %in% c(obj[["cox.phase1", exact=TRUE]], ovars))
    if (any(tmp)) {
      miss <- getQuotedVarStr(pvars[tmp])
      msg  <- paste0("NOTE: the predictors.cox.phase2 variables ", miss, 
                     " should also be in cox.phase1 or other.covars, ",
                     " putting them in other.covars.\n")
      if (print) cat(msg)
      obj[["other.covars"]] <- c(ovars, pvars)  
    }
  }

  # Check that predicted.cox.phase variables have the correct type
  predictedCovarsCheck(data, obj)

  obj
}

setup_getAllVarTypesAndMissing <- function(data, obj) {

  vars  <- cc_getAllVars(obj) 
  lst   <- list()
  miss  <- list()
  nvars <- length(vars)
  nvals <- c("cont", "binary") # for numeric variables
  for (i in 1:nvars) {
    var   <- vars[i]
    vec   <- data[, var, drop=TRUE]
    type  <- getVarType(vec, only.finite=0)
    if (type %in% nvals) {
      tmp <- !is.finite(vec)
    } else {
      tmp <- is.na(vec)
    }
    lst[[var]]  <- type
    miss[[var]] <- any(tmp)
    
    # Check for all missing
    if (all(tmp)) {
      v <- getQuotedVarStr(var)
      stop(paste0("ERROR: variable ", v, " contains all missing/non-finite values"))
    }
  }
  obj$vartype    <- lst
  obj$varmissing <- miss
  obj
}

setup_getMissingVec <- function(data, obj, var) {

  vec  <- data[, var, drop=TRUE]
  type <- obj$vartype[[var, exact=TRUE]]
  if (!length(type)) stop("INTERNAL CODING ERROR 1")
  if ((type == "cont") || (type == "binary")) {
    ret <- !is.finite(vec)
  } else {
    ret <- is.na(vec)
  } 
  ret
}

setup_phase2And3 <- function(data, obj) {

  print <- obj$print

  # Define phase2
  obj$phase2     <- NULL
  subcohortv     <- obj[["subcohort", exact=TRUE]]
  if (!is.null(subcohortv)) {
    tmp         <- (data[, subcohortv, drop=TRUE] %in% 1) |
                   (data[, obj$status, drop=TRUE] %in% 1)
    new         <- getRandomVariable("phase2_")
    data[, new] <- as.numeric(tmp)
    obj$phase2  <- new 
  } else {
    obj$phase2 <- NULL
    obj$phase3 <- NULL
    obj        <- setup_setArgsMsg(obj, "phase3", NULL, str1="since subcohort = NULL")
  }

  p2v <- obj[["phase2", exact=TRUE]]
  if (length(p2v)) {
    p2vec  <- data[, p2v, drop=TRUE] %in% 1
    p3v    <- obj[["phase3", exact=TRUE]]
    p3Flag <- !is.null(p3v)
    if (p3Flag) {
      # Check the user input phase3 variable
      p3vec <- data[, p3v, drop=TRUE] %in% 1
      tmp   <- p3vec & !p2vec
      if (any(tmp)) stop("ERROR: there are phase-three subjects not belonging to phase-two")
    } else {
      p3vec <- p2vec  
    }
 
    # At this point, we know the cox.phase2 vars contain missing values.
    # Check phase2 covs for missing values within phase2 subs
    p2covs <- obj[["cox.phase2", exact=TRUE]]
    n2     <- length(p2covs)
    if (n2) { 
      missInP2Flag <- FALSE
      for (i in 1:n2) {
        var      <- p2covs[i]
        missvec  <- setup_getMissingVec(data, obj, var)
        missInP2 <- p2vec & missvec
        if (any(missInP2)) {
          missInP2Flag <- TRUE
          # Check that these subs are not in phase3 (if user input phase 3), or
          #  assign phase3 vec for the non-missing subs
          if (p3Flag) {
            tmp <- p3vec & missInP2
            if (any(tmp)) stop("ERROR: phase-three subjects with missing data for cox.phase2 variables")
          } else {
            p3vec <- p3vec & !missInP2
          }
        } 
      }
      if (missInP2Flag) {
        if (!p3Flag) {
          # Add a phase 3 variable
          msg <- paste0("NOTE: some cox.phase2 covariates contain missing values within the phase-two sample,",
                        " a phase-three sample (phase3) will be defined assuming the missingness is at random.\n")
          if (print) cat(msg)
          new         <- getRandomVariable("phase3_")
          data[, new] <- p3vec
          obj$phase3  <- new
        }

        # Set calibrated to FALSE if it is TRUE
        if (obj$calibrated) {
          obj$calibrated <- FALSE
          if (print) cat("NOTE: setting calibrated to FALSE, since phase3 is not NULL.\n")
        }
      }
    }
  }

  list(data=data, obj=obj)
}

setup_setArgsMsg <- function(obj, arg.name, correctValue, str1="") {

  arg.val <- obj[[arg.name, exact=TRUE]]
  flag1 <- length(correctValue)
  flag2 <- length(arg.val)
  if (flag1 && flag2) {
    ok <- arg.val == correctValue
  } else if (!flag1 && !flag2) {
    ok <- TRUE
  } else {
    ok <- FALSE
  }
  if (!ok) {
    cval0 <- correctValue
    if (!flag1) correctValue <- "NULL"
    correctValue <- as.character(correctValue)
    msg <- paste0("NOTE: setting ", arg.name, " to ", correctValue, " ", str1, "\n")
    if (obj$print) cat(msg)
    obj[[arg.name]] <- cval0
  }
  obj
}

setup_check_args1 <- function(obj) {

  print <- obj$print
  phase2Flag <- !is.null(obj[["subcohort", exact=TRUE]])
  if (!phase2Flag) {
    # Check and set arguments to NULL
    null.args <- c("other.covars", "strata", "weights.phase2", "subcohort.strata.counts",
                   "predicted.cox.phase2", "predictors.cox.phase2", "aux.vars",
                   "phase3", "strata.phase3", "weights.phase3", "weights.op")
    for (arg in null.args) obj <- setup_setArgsMsg(obj, arg, NULL, str1="since subcohort = NULL")
    false.args <- c("calibrated", "predict")
    for (arg in false.args) obj <- setup_setArgsMsg(obj, arg, FALSE, str1="since subcohort = NULL")
  }
  if (!obj$calibrated) {
    null.args <- c("other.covars", "predicted.cox.phase2", "predictors.cox.phase2", "aux.vars")
    for (arg in null.args) obj <- setup_setArgsMsg(obj, arg, NULL, str1="since calibrated = FALSE")
    false.args <- c("predict")
    for (arg in false.args) obj <- setup_setArgsMsg(obj, arg, FALSE, str1="since calibrated = FALSE")
  }  
  if (!obj$predict) {
    null.args <- c("predictors.cox.phase2")
    for (arg in null.args) obj <- setup_setArgsMsg(obj, arg, NULL, str1="since predict = FALSE")

    if (obj$calibrated) {
      # Must have aux.vars or predicted.cox.phase2
      if (is.null(obj[["aux.vars", exact=TRUE]]) && is.null(obj[["predicted.cox.phase2", exact=TRUE]])) {
        obj$calibrated <- FALSE
        if (print) cat("NOTE: calibrated is set to FALSE since predict=FALSE, aux.vars = NULL, and predicted.cox.phase2 = NULL\n")
      }
    }
  }
  # If predict=TRUE and predicted.cox.phase2 != NULL, set predict to FALSE
  if (obj$predict) {
    if (length(obj[["predicted.cox.phase2", exact=TRUE]])) {
      obj <- setup_setArgsMsg(obj, "predict", FALSE, str1="since predicted.cox.phase2 is specified")
    }
  }  

  obj
}

setup_predictors <- function(obj) {

  print <- obj$print

  p2vars <- obj[["cox.phase2", exact=TRUE]]
  if (!length(p2vars)) stop("INTERNAL CODING ERROR 1")

  # Default variables
  defvars <- c(obj[["cox.phase1", exact=TRUE]], obj[["other.covars", exact=TRUE]])
  defFlag <- length(defvars)

  # Check predictors.cox.phase2 
  predvars <- obj[["predictors.cox.phase2", exact=TRUE]]
  if (!length(predvars)) {
    if (defFlag) {
      # Set to default vars
      predvars <- defvars
      msg <- "NOTE: since predictors.cox.phase2 = NULL, it will be set to c(cox.phase1, other.covars).\n"
    } else {
      # calibration cannot be done
      obj$calibrated <- FALSE
      predvars       <- NULL
      msg            <- "NOTE: calibrated analysis will not performed as cox.phase1=NULL and other.covars=NULL.\n"
    }
    if (print) cat(msg)
  } else if (isList(predvars)) {
   # Check each phase 2 covariate
   nms <- names(predvars)
   if (is.null(nms)) stop("INTERNAL CODING ERROR 2")
   tmp <- !(p2vars %in% nms)
   if (any(tmp)) {
     miss <- p2vars[tmp]
     msg  <- paste0("NOTE: the cox.phase2 variables ", getQuotedVarStr(miss), 
                    " are not in predictors.cox.phase2.\n")
     if (print) cat(msg)
     if (defFlag) {
       msg <- "Predictors for the above variables will be set to the cox.phase1 and other.covars variables.\n"
     } else {
       obj$calibrated <- FALSE
       predvars       <- NULL
       msg            <- "NOTE: calibrated analysis will not be performed as cox.phase1=NULL and other.covars=NULL.\n"
     }
     if (print) cat(msg)
     for (v in miss) predvars[[v]] <- defvars
   }
  }
  obj$predictors.cox.phase2 <- predvars
  obj
}

setup_check_constant_vars <- function(data, obj, vars, vars.name) {

  # aux.vars is a special case, do not allow constant categorical variable
  auxFlag <- vars.name == "aux.vars"

  nvars <- length(vars)
  if (!nvars) return(obj)
  vartype <- obj$vartype
  print   <- obj$print
  ntype   <- c("cont", "binary")
  rem     <- NULL
  for (i in 1:nvars) {
    v    <- vars[i]
    type <- vartype[[v]]
    vec  <- data[, v, drop=TRUE]
    if (type %in% ntype) {
      tmp   <- is.finite(vec)
      nflag <- 1
    } else {
      tmp   <- !is.na(vec)
      nflag <- 0
    }
    vec <- unique(vec[tmp])
    if (length(vec) < 2) {
      if ( !auxFlag || (auxFlag && !nflag) ) {
        rem <- c(rem, v)
        if (print) {
          if (!auxFlag) {
            msg <- paste0("NOTE: variable ", getQuotedVarStr(v), " is constant and will be removed from ", vars.name, "\n")
          } else {
            msg <- paste0("NOTE: categorical variable ", getQuotedVarStr(v), " is constant and will be removed from ", vars.name, "\n")
          }
          cat(msg)
        }
      }
    }
  }
  if (length(rem)) {
    tmp  <- !(vars %in% rem)
    vars <- vars[tmp]
    if (!length(vars)) {
      vars <- NULL
      if (print) cat(paste0("NOTE: all variables have been removed from ", vars.name, "\n"))
    }
    obj[[vars.name]] <- vars
  }
  obj
}

setup_check_constant <- function(data, obj) {
  
  nm  <- "cox.phase1"
  obj <- setup_check_constant_vars(data, obj, obj[[nm, exact=TRUE]], nm)
  nm  <- "cox.phase2"
  obj <- setup_check_constant_vars(data, obj, obj[[nm, exact=TRUE]], nm)
  nm  <- "other.covars"
  obj <- setup_check_constant_vars(data, obj, obj[[nm, exact=TRUE]], nm)
  nm  <- "aux.vars"
  obj <- setup_check_constant_vars(data, obj, obj[[nm, exact=TRUE]], nm)

  obj

}

setup_check_args2 <- function(obj) {

  print          <- obj$print
  obj            <- setup_check_args1(obj)
  obj$phase2Flag <- !is.null(obj[["subcohort", exact=TRUE]])
  obj$phase3Flag <- !is.null(obj[["phase3", exact=TRUE]])
  if (obj$phase3Flag && !obj$phase2Flag) stop("INTERNAL CODING ERROR 1")
  if (obj$calibrated && obj$predict) obj <- setup_predictors(obj)
  if (obj$calibrated && !obj$predict) {
    if (is.null(obj[["aux.vars", exact=TRUE]]) && is.null(obj[["predicted.cox.phase2", exact=TRUE]])) {
      msg <- paste0("NOTE: calibrated analysis will not be performed as predict = FALSE,\n",
                    " predicted.cox.phase2 = NULL and aux.vars = NULL.\n")
      obj$calibrated <- FALSE
      if (print) print(msg)
    }
  }

  obj <- setup_check_args1(obj)

  stratv        <- obj[["strata", exact=TRUE]]
  obj$stratFlag <- !is.null(stratv)
  if (obj$calibrated) {
    null.args <- c("phase3", "strata.phase3")
    for (arg in null.args) obj <- setup_setArgsMsg(obj, arg, NULL, str1="since calibrated = TRUE")

    # Check args 
    aflag <- length(obj[["aux.vars", exact=TRUE]])
    if (aflag) {
      msg1 <- "NOTE: since aux.vars != NULL, aux.vars will be used as auxiliary variables.\n"
      msg2 <- "      To override this default, set aux.vars = NULL.\n"
      if (print) {
        cat(msg1)
        cat(msg2)
      }
      null.args <- c("predicted.cox.phase2", "predictors.cox.phase2")
      for (arg in null.args) obj <- setup_setArgsMsg(obj, arg, NULL, str1="since aux.vars != NULL")
      obj <- setup_setArgsMsg(obj, "predict", FALSE, str1="since aux.vars != NULL")
    }
    if (obj$predict) {
      pflag  <- length(obj[["predicted.cox.phase2", exact=TRUE]])
      pflag2 <- length(obj[["predictors.cox.phase2", exact=TRUE]])
      if (pflag && pflag2) {
        msg1 <- "NOTE: since predicted.cox.phase2 != NULL, predicted.cox.phase2 will be used as the predicted values for the phase-two covariates.\n"
        msg2 <- "      To override this default, set predicted.cox.phase2 = NULL.\n"
        if (print) {
          cat(msg1)
          cat(msg2)
        }
        obj <- setup_setArgsMsg(obj, "predictors.cox.phase2", NULL, str1="since predicted.cox.phase2 != NULL")
        obj <- setup_setArgsMsg(obj, "predict", FALSE, str1="since predicted.cox.phase2 != NULL")
      }
    }
  }
  if (!obj$phase3Flag) {
    obj <- setup_setArgsMsg(obj, "strata.phase3", NULL, str1="since phase3 = NULL")
  }
  if (obj$phase3Flag) {
    if (is.null(obj[["weights.phase3.type", exact=TRUE]])) {
      if (!is.null(obj[["weights.phase3", exact=TRUE]])) {
        obj$weights.phase3.type <- "design"
      } else {
        obj$weights.phase3.type <- "estimated"
      }
    }
  }

  obj
}

getEventTimeVar <- function(data, obj) {

  time  <- obj$time
  if (length(time) == 1) {
    ret  <- time
  } else {
    ret  <- time[2]
  }
  ret
}


getEventTime <- function(data, obj) {

  timev <- getEventTimeVar(data, obj)
  ret   <- data[, timev, drop=TRUE]

  ret
}

checkForEqualTimes <- function(data, obj) {

  # Call this function before it is sorted by time

  time <- getEventTime(data, obj) 
  tmp  <- breakTies(time, data[, obj$status, drop=TRUE])
  if (!length(tmp)) return(list(data=data))

  ret  <- tmp$time
  rows <- tmp[["rows.changed", exact=TRUE]]
  if (length(ret) != nrow(data)) stop("INTERNAL CODING ERROR")
  warning("Data contains non-unique times that will be broken randomly") 

  # Save the data that changed
  timev <- getEventTimeVar(data, obj)
  time0 <- data[, timev, drop=TRUE]
  data[, timev] <- ret

  mat <- NULL
  if (length(rows)) {
    ids <- rownames(data)[rows]
    old <- time0[rows]
    new <- data[rows, timev, drop=TRUE]
    mat <- cbind(old, new)
    colnames(mat) <- paste0(timev, ".", c("orig", "new"))
    rownames(mat) <- ids
  }

  list(data=data, changed=mat)  
}


breakTies <- function(time, status) {

  # Get all duplicated times
  tmp <- duplicated(time)
  if (!any(tmp)) return(NULL)
  dup.times <- unique(time[tmp])

  # Get event times
  status1 <- status %in% 1
  etimes  <- unique(time[status1])

  # Get the duplicated times that are also event times
  tmp <- etimes %in% dup.times
  if (!any(tmp)) return(NULL)
  dups <- etimes[tmp]

  # Get the maximum number of duplicates for a given time (for events only)
  tmp  <- (time %in% dups) & status1
  maxn <- max(table(time[tmp])) + 1

  # save original time
  ret   <- time
  time  <- unique(sort(time))
  ntime <- length(time)

  # Loop over each set of duplicates and add a small amount to the events only.
  # Return row numbers with changed data.
  all.rows <- NULL
  rowvec   <- 1:length(ret)
  for (dup in dups) {
    tmp             <- (ret == dup) & status1
    tmp[is.na(tmp)] <- FALSE
    m               <- sum(tmp)
    if (!m) next

    # Get the next time after dup 
    ii <- which(time > dup)
    if (length(ii)) {
      ii   <- ii[1]  # First time after dup
      mind <- abs(time[ii] - dup) 
    } else {
      mind <- 0.5 # dup was largest time
    }
    mind     <- mind/maxn
    vec      <- sample(1:m, m)
    add      <- mind*vec
    ret[tmp] <- ret[tmp] + add
    all.rows <- c(all.rows, rowvec[tmp])     
  }  

  list(time=ret, rows.changed=all.rows)
}

setUpData.checkStrata <- function(data, obj) {

  stratav <- obj[["strata", exact=TRUE]]
  if (!length(stratav)) return(data)

  subcohortv <- obj[["subcohort", exact=TRUE]]
  if (!length(subcohortv)) return(data)

  strata    <- data[, stratav, drop=TRUE]
  ustrata   <- unique(strata)
  nstrata   <- length(ustrata)
  subcohort <- data[, subcohortv, drop=TRUE] %in% 1
  rem       <- NULL
  for (i in 1:nstrata) {
    lev <- ustrata[i] 
    tmp <- (strata %in% lev) & subcohort
    m   <- sum(tmp)
    if (m < 1) {
      rem <- c(rem, lev)
      msg <- paste0("Stratum '", lev, "' contains too few individuals in the subcohort and will be removed.")
      warning(msg)
    }
  }
  if (length(rem)) {
    tmp  <- !(strata %in% rem)
    data <- data[tmp, , drop=FALSE]
  }
  data
}

setUpData <- function(data, obj) {

  if (nrow(data) < obj$min.nrow) stop("ERROR: too few rows in data for analysis")
  status <- data[, obj$status, drop=TRUE]
  if (sum(status %in% 0) < 1) stop("ERROR: no subjects with status = 0")
  if (sum(status %in% 1) < 1) stop("ERROR: no subjects with status = 1")

  # Check stratification
  data <- setUpData.checkStrata(data, obj)

  obj <- setup_check_args1(obj)

  # Get the variable type for each variable and determine which have missing values
  obj <- setup_getAllVarTypesAndMissing(data, obj)

  # Check for constant vars, must be called after setup_getAllVarTypesAndMissing
  obj <- setup_check_constant(data, obj)

  # Check covars and determine the phase 2 covariates
  obj <- setup_check_covars_subcohort(data, obj)

  # Check phase2 and phase3 samples 
  tmp  <- setup_phase2And3(data, obj)
  data <- tmp$data
  obj  <- tmp$obj
  rm(tmp); gc()

  obj <- setup_check_args2(obj)
  obj <- setup_check_args2(obj) # Call a second time to make sure objects are correctly set

  rn <- rownames(data)
  if (is.null(rn)) {
    rn <- as.character(1:nrow(data))
    rownames(data) <- rn
  }

  # Check for missing values
  cc_checkForMissing(data, obj) 

  # Set this after cc_checkForMissing call
  tmp               <- obj[["cox.phase2", exact=TRUE]]
  obj$covars        <- c(obj[["cox.phase1", exact=TRUE]], tmp)
  obj$phase2.covars <- tmp
  
  # Check tau
  obj <- tau12.check(data, obj)

  # Keep only columns we need
  vars <- cc_getAllVars(obj) 
  data <- data[, vars, drop=FALSE]

  # Set the strata variables
  tmp  <- setStrataVars(data, obj)
  data <- tmp$data
  obj  <- tmp$obj

  # Add dummy variables for cat vars, adds factors
  tmp  <- cc_addAllDummyVars(data, obj$covars, obj)
  data <- tmp$data
  obj  <- tmp$obj.list
  rm(tmp)
  gc()

  # Replace cat vars with dummy variables
  obj$covars.orig        <- obj$covars
  obj$covars             <- cc_updateVarsWithDummy(obj$covars, obj)
  obj$phase2.covars.orig <- obj$phase2.covars
  obj$phase2.covars      <- cc_updateVarsWithDummy(obj$phase2.covars, obj)
  
  # For times that are equal, break ties
  tmp               <- checkForEqualTimes(data, obj)
  data              <- tmp$data
  obj$changed.times <- tmp[["changed", exact=TRUE]]

  # Order by event time, or time[2] for age-scale
  obj$ntime <- length(obj$time)
  timevec   <- getEventTime(data, obj)
  ord       <- order(timevec)
  data      <- data[ord, , drop=FALSE]

  # Set phase2, phase3 to logical
  if (obj$phase2Flag) {
    data[, obj$subcohort] <- as.logical(data[, obj$subcohort, drop=TRUE])
    data[, obj$phase2]    <- as.logical(data[, obj$phase2, drop=TRUE])
  }
  if (obj$phase3Flag) data[, obj$phase3] <- as.logical(data[, obj$phase3, drop=TRUE])

  # Get stratification counts. Must be called after new strata are defined
  if (obj$phase2Flag) {
    obj$strata.counts <- getStratCounts(data, obj) 
    obj$n.strata      <- length(obj$strata.counts$cnts.cohort)
    obj$n.phase2      <- sum(data[, obj$phase2, drop=TRUE])

    # Get phase 2 weights
    if (is.null(obj[["weights.phase2", exact=TRUE]])) {
      wgtv               <- getRandomVariable("wgtp2_")
      data[, wgtv]       <- getPhase2Weights(data, obj) 
      obj$weights.phase2 <- wgtv
    }
  }
  if (obj$phase3Flag) {
    obj$strat3Flag <- !is.null(obj[["strata.phase3", exact=TRUE]])

    # Compute weights.phase3 if needed
    if (is.null(obj[["weights.phase3", exact=TRUE]])) {
      wgtv               <- getRandomVariable("wgtp3_")
      data[, wgtv]       <- 1
      tmp                <- data[, obj$phase3, drop=TRUE]
      data[tmp, wgtv]    <- computePhase3Weights(data, obj)  
      obj$weights.phase3 <- wgtv 
    }
  }

  # Add id variable for coxph
  idv         <- getRandomVariable("id_")
  data[, idv] <- 1:nrow(data) 
  obj$id.var  <- idv

  obj$nrow.cohort <- nrow(data)
  obj$shin.method <- obj$aux.method == "shin"
  if (obj$phase3Flag && !obj$phase2Flag) stop("INTERNAL CODING ERROR 9")

  weights.op <- check_calibrated.op(obj[["weights.op", exact=TRUE]]) 

  if (obj$print) {
    str1 <- ifelse(obj$stratFlag, "A stratified ", "An unstratified ")
    str2 <- ifelse(obj$phase3Flag, "phase-three ", "phase-two ")
    if (!obj$phase2Flag) {
      str <- "A whole cohort analysis will be run\n"
    } else {
      str <- paste0(str1, str2, "analysis will be run\n")
    }
    cat(str)
  }

  # control list for coxph
  obj$coxph.control <- coxph.control(timefix=FALSE)

  list(data=data, rownames0=rn, obj.list=obj)
}

