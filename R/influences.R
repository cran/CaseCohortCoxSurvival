getInfluences <- function(data, obj, A, fit) {

  if (obj$print) cat("Computing influences\n")
  beta       <- fit$coefficients
  fit.detail <- coxph.detail(fit, riskmat=FALSE)
  Tau1       <- obj[["tau1", exact=TRUE]]
  Tau2       <- obj[["tau2", exact=TRUE]]
  weights    <- fit[["weights", exact=TRUE]]
  if (is.null(weights)) weights <- rep(1, nrow(getFitDetailX(fit.detail)))
  #names(weights) <- rownames(fit$y)

  infl   <- NULL
  infl.F <- NULL
  infl.T <- NULL

  if (!obj$phase3Flag) {
    resid <- residuals(fit, type="dfbeta", weighted=TRUE)
    infl  <- cc_influences(fit.detail, beta, weights=weights, Tau1=Tau1, 
                  Tau2=Tau2, calibrated=obj$calibrated, resid=resid, A=A, DEBUG=obj$DEBUG)
  } else {
    ymat.p2      <- getCoxDetailYmat.phase2(data, obj)
    type         <- obj$weights.phase3.type
    est.T        <- type %in% c("both", "estimated")
    est.F        <- type %in% c("both", "design")
    B.phase2     <- NULL

    if (est.T) {
      # Get B.phase2
      B.phase2 <- getB.phase2(data, obj)
      infl.T   <- cc_influences_miss(fit.detail, beta, weights, ymat.p2, 
                       Tau1, Tau2, estimated.weights=TRUE, B.phase2=B.phase2, est.op=obj$weights.op)
    }
    B.phase2 <- NULL
    if (est.F) {
      # Get B.phase2
      B.phase2 <- getB.phase2(data, obj)
      infl.F   <- cc_influences_miss(fit.detail, beta, weights, ymat.p2, 
                       Tau1, Tau2, estimated.weights=FALSE, B.phase2=B.phase2, est.op=obj$weights.op)
    }
  }

  list(infl=infl, infl.T=infl.T, infl.F=infl.F)
}

getFitDetailX <- function(fit.detail) {

  x <- fit.detail$x
  d <- dim(x)
  if (is.null(d)) {
    nms         <- names(x)
    dim(x)      <- c(length(x), 1)
    rownames(x) <- nms
  }
  x
}

getStatusTimesCoxDetail <- function(y) {

  # y is coxph.detail$y matrix

  status <- as.numeric(y[, "status", drop=TRUE])
  if ("time" %in% colnames(y)) {
    time  <- as.numeric(y[, "time", drop=TRUE])
    time2 <- NULL
  } else {
    time  <- as.numeric(y[, "start", drop=TRUE])
    time2 <- as.numeric(y[, "stop", drop=TRUE])
  }
  list(status=status, time=time, time2=time2)
}

getCoxDetailYmat.phase2 <- function(data, obj) {

  # Data was sorted in the beginning 

  tmp    <- data[, obj$phase2, drop=TRUE]
  time   <- obj$time
  ret    <- data[tmp, c(obj$status, time), drop=FALSE]
  ret    <- as.matrix(ret)
  if (length(time) == 1) {
    cx <- c("status", "time")
  } else {
    cx <- c("status", "start", "stop")
  }
  colnames(ret) <- cx

  ret
}

cc_getTimeId <- function(alltimes, eventTimes) {

  n      <- length(alltimes)
  ret    <- rep(0, n)
  ntimes <- length(eventTimes)
  if (!ntimes) return(ret)
  for (i in 1:ntimes) {
    tmp <- alltimes >= eventTimes[i]
    ret <- ret + as.numeric(tmp)
  }
  ret  

}

cc_getTime0Id <- function(alltimes, eventTimes) {

  n      <- length(alltimes)
  ntimes <- length(eventTimes)
  ret    <- rep(ntimes, n)
  for (i in 1:ntimes) {
    tmp <- alltimes < eventTimes[i]
    ret <- ret - as.numeric(tmp)
  }
  ret  

}

cc_influences  <- function (mod.detail, beta.hat, weights=NULL, Tau1=NULL, Tau2=NULL, 
                            calibrated=FALSE, A=NULL, resid=NULL, DEBUG=0) {
  
  # Sort original data by event time

  if (is.null(calibrated)) calibrated <- FALSE
  if (calibrated && is.null(A)) stop("A matrix with auxiliary variables values need to be provided")

  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  X                 <- getFitDetailX(mod.detail)
  n                 <- nrow(X)
  p                 <- ncol(X) 
  if (is.null(rownames(X))) stop("ERROR: design matrix does not have rownames")
  #if (is.null(names(weights))) stop("ERROR: weights do not have names")
  if (!is.null(A)) {
    tmp <- !(rownames(X) %in% rownames(A))
    if (any(tmp)) stop("ERROR with rownames(A)") 
  }
  #tmp <- !(rownames(X) %in% names(weights))
  #if (any(tmp)) stop("ERROR with names(weights)") 
  #weights <- weights[rownames(X)]
  #if (length(resid)) {
  #  if (is.null(rownames(resid))) stop("ERROR: resid does not have rownames")
  #  tmp <- !(rownames(X) %in% rownames(resid))
  #  if (any(tmp)) stop("ERROR with rownames(resid)") 
  #  resid <- resid[rownames(X), , drop=FALSE]
  #}


  mod.detail.time   <- mod.detail$time 
  number.times      <- length(mod.detail.time)
  tmp               <- getStatusTimesCoxDetail(mod.detail$y)
  y                 <- tmp$status
  time              <- tmp$time
  time2             <- tmp[["time2", exact=TRUE]]
  imat              <- mod.detail$imat
  if (!is.null(dim(imat))) {
    infomat         <- apply(X=imat, MARGIN=c(1,2), FUN=sum)
  } else  {
    infomat         <- sum(imat)
  }
  infomat.inv       <- solve(infomat)
  rm(tmp, imat, infomat)
  gc()
  
  if (is.null(Tau1)) Tau1 <- floor(min(mod.detail.time))
  if (is.null(Tau2)) Tau2 <- floor(max(mod.detail.time))
  Tau1Tau2.times <- which((Tau1 < mod.detail.time) & (mod.detail.time <= Tau2)) 
  if (!length(Tau1Tau2.times)) stop("ERROR: no event times in (tau1, tau2], adjust tau1 and/or tau2")
  if (is.null(weights)) weights <- rep(1, n)
  
  exp.X.weighted       <- weights*exp(X %*% beta.hat)
  if (is.null(time2)) {
    obstimes           <- cc_getTimeId(time, mod.detail.time)
    obstimes0          <- NULL
  } else {
    obstimes           <- cc_getTimeId(time2, mod.detail.time)
    obstimes0          <- cc_getTime0Id(time, mod.detail.time)
    if (any(obstimes < obstimes0)) stop("INTERNAL CODING ERROR")
  }

  S0t                  <- cc_getS0t(obstimes, obstimes0, number.times, exp.X.weighted, DEBUG=DEBUG)
  if (any(S0t == 0)) stop("ERROR: S0t = 0")
  S1t                  <- cc_getS1t(obstimes, obstimes0, exp.X.weighted, X, number.times, DEBUG=DEBUG)
  lambda0.t.hat        <- cc_getdNtColSums(obstimes, number.times, y, DEBUG=DEBUG)/S0t
  Lambda0.Tau1Tau2.hat <- sum(lambda0.t.hat[Tau1Tau2.times])
  rm(time, time2)
  gc()

  # ----------------------------------------------------------------------------
  # If calibrated weights ------------------------------------------------------
  
  if (calibrated) {

    # --------------------------------------------------------------------------
    # Computation of the influences for the Lagrangian multipliers, eta --------
    A.phase2            <- A[rownames(X), , drop=FALSE]
    q                   <- ncol(A.phase2) 
    A.weighted          <- sweep(A.phase2, 1, weights, "*")
    sum.AA.weighted     <- cc_getSumAAwgt(A.phase2, A.weighted, DEBUG=DEBUG)
    sum.AA.weighted.inv <- solve(sum.AA.weighted)
    
    infl1.eta               <- A %*% sum.AA.weighted.inv
    infl2.eta               <- 0 * A
    infl2.eta[rownames(X),] <- - A.weighted %*% sum.AA.weighted.inv
    infl.eta                <- infl1.eta + infl2.eta 
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    tmp <- cc_get_drond_U_eta(X, A.phase2, obstimes, obstimes0, number.times, y, 
                              S1t, S0t, exp.X.weighted, weights, DEBUG=DEBUG)
    rm(A.phase2)
    gc()
    drond.S0t.eta    <- tmp$drond.S0t.eta
    drond.U.eta      <- tmp$drond.U.eta
    infl1.beta       <- infl1.eta %*% t(drond.U.eta) %*% infomat.inv
    score.beta       <- cc_getBetaScore(X, obstimes, obstimes0, number.times, y, 
                                        weights, exp.X.weighted, S1t, S0t, DEBUG=DEBUG)
    infl2.beta       <- infl2.eta %*% t(drond.U.eta) %*% infomat.inv
    tmp              <- rownames(X)
    infl2.beta[tmp,] <- infl2.beta[tmp,] + score.beta %*% infomat.inv
    infl.beta        <- infl1.beta + infl2.beta

    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    rn  <- rownames(X)
    tmp <- rownames(A) %in% rownames(X)
    rm(X)
    gc()
    tmp <- cc_infl_lambda0_tau12(infl1.beta, S1t, infl1.eta, drond.S0t.eta, 
                                 infl2.beta[rn, , drop=FALSE], infl2.eta[rn, , drop=FALSE], 
                                 lambda0.t.hat, S0t, obstimes, obstimes0, y, exp.X.weighted, 
                                 number.times, Tau1Tau2.times, tmp, DEBUG=DEBUG)
    rm(obstimes, y, exp.X.weighted, Tau1Tau2.times, drond.S0t.eta, S0t, S1t)
    gc()
    infl1.Lambda0.Tau1Tau2 <- tmp$infl1.Lambda0.Tau1Tau2
    infl2.Lambda0.Tau1Tau2 <- tmp$infl2.Lambda0.Tau1Tau2
    rm(tmp, rn)
    gc()
    infl.Lambda0.Tau1Tau2  <- infl1.Lambda0.Tau1Tau2 + infl2.Lambda0.Tau1Tau2 

    return(list(infl.beta = as.matrix(infl.beta), infl1.beta=as.matrix(infl1.beta),
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl2.beta = as.matrix(infl2.beta), 
                infl1.Lambda0.Tau1Tau2 = as.matrix(infl1.Lambda0.Tau1Tau2),
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }else{
    # --------------------------------------------------------------------------
    # If design weights --------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    # infl.beta = resid

    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    infl.Lambda0.Tau1Tau2 <- cc_infl_lambda0_tau12_noCalib(resid, S1t, 
                          lambda0.t.hat, S0t, obstimes, obstimes0, y, exp.X.weighted,                      
                          number.times, Tau1Tau2.times, DEBUG=DEBUG) 

    rm(X, S0t, S1t, obstimes, y, exp.X.weighted, Tau1Tau2.times)
    gc()
    
    return(list(infl.beta = as.matrix(resid), 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }
}

