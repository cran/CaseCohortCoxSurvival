getVarianceForEst <- function(data, obj, infl, fit) {

  if (obj$print) cat("Computing variance estimates for parameters\n")

  var.lam <- NULL
  varlst  <- NULL
  varlstT <- NULL
  varlstF <- NULL

  if (!obj$phase2Flag) {
    var.beta <- getCoxphVar(fit)
    # Set to NA, perhaps change later
    var.beta[] <- NA
    varlst     <- list(var.beta=var.beta)
  } else if (obj$phase3Flag) {
    estT     <- !is.null(infl[["infl.T", exact=TRUE]])
    estF     <- !is.null(infl[["infl.F", exact=TRUE]])
    data.p2  <- data[data[, obj$phase2, drop=TRUE], , drop=FALSE]
    if (estT) {
      x        <- infl[["infl.T"]]
      var.beta <- cc_var_miss(data.p2, obj, x$infl2.beta, x$infl3.beta, estimated.weights=TRUE)
      var.lam  <- cc_var_miss(data.p2, obj, x$infl2.Lambda0.Tau1Tau2, 
                              x$infl3.Lambda0.Tau1Tau2, estimated.weights=TRUE)
      varlstT  <- list(var.beta=var.beta, var.lam=var.lam)
    }
    if (estF) {
      x        <- infl[["infl.F"]]
      var.beta <- cc_var_miss(data.p2, obj, x$infl2.beta, x$infl3.beta, estimated.weights=FALSE)
      var.lam  <- cc_var_miss(data.p2, obj, x$infl2.Lambda0.Tau1Tau2, 
                              x$infl3.Lambda0.Tau1Tau2, estimated.weights=FALSE)
      varlstF  <- list(var.beta=var.beta, var.lam=var.lam)
    }
  } else {
    phase2  <- data[, obj$phase2, drop=TRUE]
    status  <- data[, obj$status, drop=TRUE]
    wgtv <- obj$weights.phase2
    if (!length(wgtv)) stop("INTERNAL CODING ERROR 1")
    weights <- data[, wgtv, drop=TRUE]
    strata  <- getStrataVec(data, obj)
    x       <- infl[["infl", exact=TRUE]]
    if (is.null(x)) stop("INTERNAL CODING ERROR 2")

    var.beta <- cc_variance(x$infl.beta, status, phase2, weights=weights, 
                            infl2=x[["infl2.beta", exact=TRUE]], 
                            strata=strata, calibrated=obj$calibrated, 
                            variance.phase2=FALSE, strata.counts=obj$strata.counts)
    var.lam  <- cc_variance(x$infl.Lambda0.Tau1Tau2, status, phase2, 
                            weights=weights, infl2=x[["infl2.Lambda0.Tau1Tau2", exact=TRUE]], 
                            strata=strata, calibrated=obj$calibrated, 
                            variance.phase2=FALSE, strata.counts=obj$strata.counts)
    varlst  <- list(var.beta=var.beta, var.lam=var.lam)
  }
    
  
  list(var=varlst, var.T=varlstT, var.F=varlstF)
}

getVarianceForData <- function(risk.data, data, obj, infl) {

  if (is.null(risk.data)) return(NULL)
  
  covars   <- obj$covars
  infl.nms <- names(infl)
  if (!length(infl.nms)) stop("INTERNAL CODING ERROR 1")
  beta <- getObjFromInfl(infl, "beta.hat")
  if (!length(beta)) stop("INTERNAL CODING ERROR 2")

  rv   <- NULL
  rv.T <- NULL
  rv.F <- NULL

  if (obj$print) cat("Computing variance estimates for pure risk\n")

  # Get matrix of risk data. Original covars must be used, because some
  #  may be categorical
  rmat <- setUpRiskData(risk.data, obj$covars.orig, beta, obj) 
  if (!length(rmat)) return(NULL)

  # For phase3, same call as with calibrated=TRUE. Different infl lists for
  #   estimate=TRUE or FALSE
  if (obj$phase3Flag) {
    obj$calibrated <- TRUE 
    x <- infl[["infl.T", exact=TRUE]] 
    if (!is.null(x)) {
      obj$estimated.weights <- TRUE
      rv.T <- cc.var.risk(rmat, data, obj, x)
    }
    x <- infl[["infl.F", exact=TRUE]] 
    if (!is.null(x)) {
      obj$estimated.weights <- FALSE
      rv.F <- cc.var.risk(rmat, data, obj, x)
    }
  } else {
    x <- infl[["infl", exact=TRUE]]
    if (is.null(x)) stop("INTERNAL CODING ERROR 3")    
    rv <- cc.var.risk(rmat, data, obj, x)
  }
  list(rv=rv, rv.T=rv.T, rv.F=rv.F)

}

cc.var.risk <- function(riskmat, data, obj, infl) {

  # Get unique covariate patterns
  xids  <- myPasteCols(riskmat, sep=":")
  uids  <- unique(xids)
  rows  <- match(xids, uids)
  tmp   <- !duplicated(xids)
  x     <- riskmat[tmp, , drop=FALSE] 
  gc()

  tmp      <- oneMinusPiExpx(x, infl$beta.hat, infl$Lambda0.Tau1Tau2.hat)
  obj1     <- tmp$ret
  Pi.x.hat <- tmp$Pi.x.hat
  obj2     <- obj1*infl$Lambda0.Tau1Tau2.hat
  # obj1 and obj2 are vectors 

  if (obj$calibrated) {
    ret <- cc.var.risk.Calib(x, data, obj, infl, obj1, obj2) 
  } else {
    ret <- cc.var.risk.noCalib(x, data, obj, infl, obj1, obj2)
  }
 
  # ret contains variance and robust variance estimates
  ret           <- cbind(Pi.x.hat, ret$var, ret$robust.var)
  colnames(ret) <- c("risk", "variance", "robust.variance")
  ret           <- ret[rows, , drop=FALSE]
  nms           <- rownames(riskmat)
  if (!is.null(nms)) rownames(ret) <- nms
  ret

}

cc.var.risk.noCalib <- function(x, data, obj, infl, obj1, obj2) {

  phase2Flag            <- obj$phase2Flag
  missing               <- rowSums(!is.finite(x))
  nx                    <- nrow(x)
  ret                   <- rep(NA, nx)
  tx                    <- t(x)
  infl.beta             <- infl[["infl.beta", exact=TRUE]]
  if (!length(infl.beta)) stop("INTERNAL CODING ERROR 1")
  infl.Lambda0.Tau1Tau2 <- infl[["infl.Lambda0.Tau1Tau2", exact=TRUE]]
  if (!length(infl.Lambda0.Tau1Tau2)) stop("INTERNAL CODING ERROR 2")
  if (phase2Flag) subcohort <- data[, obj$phase2, drop=TRUE]
  status                <- data[, obj$status, drop=TRUE]
  weights               <- data[, obj$weights.phase2, drop=TRUE]
  strata                <- getStrataVec(data, obj)
  strata.counts         <- obj$strata.counts
  ret.rvar              <- rep(NA, nx)
  ret.var               <- rep(NA, nx)

  for (i in 1:nx) {
    if (!missing[i]) {
      infl.x      <- obj2[i]*infl.beta %*% tx[, i, drop=FALSE] + obj1[i]*infl.Lambda0.Tau1Tau2
      ret.rvar[i] <- robustvariance(infl.x)

      if (phase2Flag) { 
        ret.var[i]  <- cc_variance(infl.x, status, subcohort, weights=weights, strata=strata, 
                         calibrated=FALSE, variance.phase2=FALSE, strata.counts=strata.counts)
      }
    }
  }
  list(var=ret.var, robust.var=ret.rvar)
}

cc.var.risk.Calib <- function(x, data, obj, infl, obj1, obj2) {

  phase3Flag <- obj$phase3Flag
  subcohort  <- data[, obj$phase2, drop=TRUE]
  if (!phase3Flag) {
    infl1.beta             <- infl$infl1.beta 
    infl1.Lambda0.Tau1Tau2 <- infl$infl1.Lambda0.Tau1Tau2 
    infl2.beta             <- infl$infl2.beta 
    infl2.Lambda0.Tau1Tau2 <- infl$infl2.Lambda0.Tau1Tau2 
    status                 <- data[, obj$status, drop=TRUE]
    weights                <- data[, obj$weights.phase2, drop=TRUE]
    strata                 <- getStrataVec(data, obj)
    strata.counts          <- obj$strata.counts 
  } else {
    estWgt                 <- obj$estimated.weights
    data.p2                <- data[subcohort, , drop=FALSE]
    infl1.beta             <- infl$infl2.beta 
    infl1.Lambda0.Tau1Tau2 <- infl$infl2.Lambda0.Tau1Tau2 
    infl2.beta             <- infl$infl3.beta 
    infl2.Lambda0.Tau1Tau2 <- infl$infl3.Lambda0.Tau1Tau2 
  }
  missing  <- rowSums(!is.finite(x))
  nx       <- nrow(x)
  ret.rvar <- rep(NA, nx)
  ret.var  <- rep(NA, nx)
  tx       <- t(x)
  
  
  for (i in 1:nx) {
    if (!missing[i]) {
      # Get the influences for each row of risk data x 
      obj2i       <- obj2[i]
      obj1i       <- obj1[i]   
      xvec        <- tx[, i, drop=FALSE]
      infl1.Pi.x  <- obj2i*infl1.beta %*% xvec + obj1i*infl1.Lambda0.Tau1Tau2
      infl2.Pi.x  <- obj2i*infl2.beta %*% xvec + obj1i*infl2.Lambda0.Tau1Tau2
      infl.Pi.x   <- infl1.Pi.x + infl2.Pi.x
      ret.rvar[i] <- robustvariance(infl.Pi.x)
      if (!phase3Flag) {
        ret.var[i] <- cc_variance(infl.Pi.x, status, subcohort, weights=weights, infl2=infl2.Pi.x, strata=strata, 
                                  calibrated=TRUE, variance.phase2=FALSE, strata.counts=strata.counts)
      } else {
        # infl2 and infl3 becomes infl1 and infl2
        ret.var[i] <- cc_var_miss(data.p2, obj, infl1.Pi.x, infl2.Pi.x, estimated.weights=estWgt)
      }
    }
  }
  list(var=ret.var, robust.var=ret.rvar)
}

cc_variance <- function(infl, status, subcohort, weights=NULL, infl2=NULL, strata=NULL, 
                     calibrated=FALSE, variance.phase2=FALSE, strata.counts=NULL) {
  
  stratified <- !is.null(strata)
  if (is.null(weights)) stop("ERROR: weights must be specified")
  if (is.null(strata.counts)) stop("INTERNAL CODING ERROR")

  # Weights must have same length as infl rows
  nri <- nrow(infl)

  # Whole cohort uses infl for phase-2 var, phase-2 uses infl2
  if (!length(infl2)) {
    phase2.var <- cc_getPhase2Var(subcohort, status, infl, strata,
                                strata.counts$cnts.casecohort, strata.counts$cnts.cohort)
  } else {
    phase2.var <- cc_getPhase2Var(subcohort, status, infl2, strata,
                                strata.counts$cnts.casecohort, strata.counts$cnts.cohort)
  }

  n <- length(status)

  if (!calibrated) {
    if (length(weights) > nri) weights <- weights[subcohort]
    if (length(weights) != nri) stop("INTERNAL CODING ERROR")
    superpop.var  <- n / (n - 1) * crossprod((infl / sqrt(weights)))
    var           <- diag(superpop.var) + phase2.var
  } else { 
    omega            <- rep(1, n)
    omega[subcohort] <- weights[subcohort]
    superpop.var     <- n / (n - 1) * colSums((infl - infl2) ^ 2 + 2 * 
                                               (infl - infl2) * infl2 + 
                                               (infl2 / sqrt(omega)) ^ 2)
    var              <- superpop.var + phase2.var
  }
  
  if (!variance.phase2) {
    return(var)
  } else {
    return(list(variance = var, variance.phase2 = diag(phase2.var)))
  }
}

getCoxphVar <- function(fit) {

  mat <- summary(fit)$coefficients
  ret <- mat[,  "se(coef)", drop=TRUE]
  ret <- ret*ret
  names(ret) <- rownames(mat)
  ret
}

