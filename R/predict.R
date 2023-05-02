getPredictedData <- function(data, obj) {

  if (obj$predict) {
    data <- predictData(data, obj)
  } else {
    # First use aux.vars, then predicted.cox.phase2,  
    predvars <- obj[["predicted.cox.phase2", exact=TRUE]]
    if (is.null(predvars)) stop("INTERNAL CODING ERROR 1")
    if (!isList(predvars)) stop("INTERNAL CODING ERROR 2")
    covars.p2 <- obj[["phase2.covars.orig", exact=TRUE]]
    if (is.null(covars.p2)) stop("INTERNAL CODING ERROR 3")
    catvars <- obj[["catvars", exact=TRUE]]
    catFlag <- length(catvars)
    for (v in covars.p2) {
      predv <- predvars[[v, exact=TRUE]]
      if (length(predv) != 1) stop("INTERNAL CODING ERROR 4")
      data[, v] <- data[, predv, drop=TRUE]
      if (catFlag && (v %in% catvars)) data[, v] <- factor(data[, v, drop=TRUE])
    }
  }

  data
}

getVarTypes <- function(data, vars, only.finite=0) {

  vars  <- unique(vars)
  nvars <- length(vars)
  if (!length(vars)) return(NULL)
  ret <- list()
  for (i in 1:nvars) {
    var        <- vars[i]
    vec        <- data[, var, drop=TRUE]
    type       <- getVarType(vec, only.finite=only.finite)
    ret[[var]] <- type
  }
  ret
}

predictData <- function(data, obj) {

  # Predict phase2 covariates from linear, logistic, or multinomial
  # Each model is fit in the phase2 sample

  print <- obj$print
  if (print) cat("Computing predicted values for phase-two covariates\n")

  # Check which of aux.vars and predictors.cox.phase2 are defined
  pclist <- obj[["aux.vars", exact=TRUE]]
  if (!is.null(pclist)) {
    if (print) cat("NOTE: aux.vars will be used as predictors for the phase-two covariates.\n")
  } 
  if (is.null(pclist)) {
    # Try predictors.cox.phase2, must not be NULL
    pclist <- obj[["predictors.cox.phase2", exact=TRUE]]
    if (is.null(pclist)) stop("INTERNAL CODING ERROR 1")
    if (print) cat("NOTE: predictors.cox.phase2 will be used as predictors for the phase-two covariates.\n")
  }
  covars.p2 <- obj[["phase2.covars.orig", exact=TRUE]]
  if (is.null(covars.p2)) stop("INTERNAL CODING ERROR 2")
 
  if (is.null(pclist)) pclist <- list()
  if (isList(pclist)) {
    pclistFlag <- 1
  } else {
    pclistFlag <- 0
    pclist     <- unlist(pclist)
  }
  subset     <- getSubsetVec(data, obj) 
  vartype    <- obj$vartype
  wgts       <- data[, obj$weights.phase2, drop=TRUE]
  ret        <- data
  dvars.list <- obj[["dummyvars", exact=TRUE]] 
  for (var in covars.p2) {
    if (pclistFlag) {
      vv <- pclist[[var, exact=TRUE]]
    } else {
      vv <- pclist # could be a vector of proxy vars
    }
    if (!length(vv)) {
      msg <- paste0("ERROR: no covariates available to predict phase-2 covariate ",  var, ",\n",
                    " specify predictors.cox.phase2")
      stop(msg)
    }
    fstr <- getGlmFormula(var, vv)
    if (print > 1) {
      msg <- paste0("Model formula for predicting phase-two covariate ", var, " is:\n",
                    fstr, "\n")
      cat(msg)
    }
    form <- as.formula(fstr)
    type <- vartype[[var, exact=TRUE]]
    if (is.null(type)) stop("INTERNAL CODING ERROR 3")

    catFlag <- FALSE
    if (type == "binary") {
      tmp <- predictVar_binary(data, subset, form, wgts) 
      fit <- tmp$fit
      obj <- tmp$predicted
    } else if (type == "cont") {
      tmp <- predictVar_cont(data, subset, form, wgts) 
      fit <- tmp$fit
      obj <- tmp$predicted
    } else {
      # Categorical
      tmp     <- predictVar_cat(data, subset, form, wgts) 
      fit     <- tmp$fit
      obj     <- tmp$predicted
      catFlag <- TRUE
    }
    ret[, var] <- obj
    if (catFlag) ret[, var] <- factor(ret[, var, drop=TRUE])
    if (print > 2) {
      cat("*********************************************\n")
      cat(paste0("Fitted model for covariate ", var, ":\n"))
      cat("*********************************************\n")
      print(fit)
    }
  }

  ret
}

predictVar_cont <- function(data, subset, form, wgts) {

  fit <- suppressWarnings(glm(form, family=gaussian(), data=data[subset, , drop=FALSE], weights=wgts[subset]))
  vec <- predict(fit, type="response", newdata=data)
  if (length(vec) != nrow(data)) stop("ERROR 1")

  list(fit=fit, predicted=vec)
}

predictVar_binary <- function(data, subset, form, wgts) {

  mu   <- mean(wgts[subset], na.rm=TRUE)
  wgts <- wgts/mu
  fit  <- suppressWarnings(glm(form, family=binomial(), data=data[subset, , drop=FALSE], weights=wgts[subset])) 
  vec  <- predict(fit, type="response", newdata=data)
  if (length(vec) != nrow(data)) stop("ERROR 1")
  list(fit=fit, predicted=vec)
}

predictVar_cat <- function(data, subset, form, wgts) {

  fit <- suppressWarnings(nnet::multinom(form, data[subset, , drop=FALSE], wgts[subset], trace=FALSE)) 
  vec <- predict(fit, newdata=data)
  if (length(vec) != nrow(data)) stop("ERROR 1")
  list(fit=fit, predicted=vec)
}

getVarType <- function(vec, only.finite=0) {

  ret <- NULL
  if (is.character(vec) || is.factor(vec)) {
    ret <- "cat"
  } else if (is.numeric(vec) || is.logical(vec)) {
    if (only.finite) vec <- vec[is.finite(vec)]  
    if (all(vec %in% 0:1)) {
      ret <- "binary"
    } else {
      ret <- "cont"
    }
  } else {
    stop("ERROR 1")
  }  
  ret
}

