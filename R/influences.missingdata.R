getB.phase2 <- function(data, obj) {

  phase2 <- data[, obj$phase2, drop=TRUE]
  if (!is.logical(phase2)) stop("INTERNAL CODING ERROR")
  n <- sum(phase2)
  strat3v <- obj[["strata.phase3", exact=TRUE]]
  if (length(strat3v)) {
    strat  <- data[phase2, strat3v, drop=TRUE]
    ustrat <- sort(unique(strat)) 
    nstrat <- length(ustrat)
    ret    <- matrix(data=NA, nrow=n, ncol=nstrat)
    for (i in 1:nstrat) ret[, i] <- as.numeric(strat %in% ustrat[i]) 
  } else {
    ret <- matrix(data=1, nrow=n, ncol=1)
  }
  rownames(ret) <- rownames(data[phase2, , drop=FALSE])
  ret
}

cc_influences_miss <- function(mod.detail, beta.hat, weights, ymat.p2,
                               Tau1, Tau2, estimated.weights=FALSE, B.phase2=NULL, est.op=NULL) {

  # mod.detail           from phase 3 data
  # ymat.p2              cox.detail$y matrix for phase 2 subjects
  # indiv.phase2         rownames for phase 2 subs

  if (estimated.weights && is.null(B.phase2)) {
    stop("B.phase2: Strata for the third phase of sampling (missingness) must be provided")
  }
  if (!length(est.op)) est.op <- list()

  X            <- getFitDetailX(mod.detail)
  n            <- nrow(X)
  p            <- ncol(X) 
  indiv.phase3 <- rownames(X)
  if (is.null(indiv.phase3)) stop("ERROR: design matrix does not have rownames")
  indiv.phase2 <- rownames(ymat.p2)

  # Observed times for phase2
  mat                    <- ymat.p2
  tmp                    <- getStatusTimesCoxDetail(mat)
  y.phase2               <- tmp$status
  time.phase2            <- tmp$time
  time2.phase2           <- tmp[["time2", exact=TRUE]]
  tmp                    <- y.phase2 %in% 1 
  if (is.null(time2.phase2)) {
    mod.detail.time.phase2 <- sort(unique(time.phase2[tmp]))
  } else {
    mod.detail.time.phase2 <- sort(unique(time2.phase2[tmp]))
  }
  number.times.phase2    <- length(mod.detail.time.phase2)
  if (is.null(time2.phase2)) {
    obstimes.phase2      <- cc_getTimeId(time.phase2, mod.detail.time.phase2)
    obstimes0.phase2     <- NULL
  } else {
    obstimes.phase2      <- cc_getTimeId(time2.phase2, mod.detail.time.phase2)
    obstimes0.phase2     <- cc_getTime0Id(time.phase2, mod.detail.time.phase2)
  } 
  rm(time.phase2, tmp, time2.phase2)
  gc()

  subset.phase3          <- match(indiv.phase3, rownames(mat))
  if (any(is.na(subset.phase3))) stop("INTERNAL CODING ERROR 1")
  mat                    <- mat[indiv.phase3, , drop=FALSE]
  tmp                    <- getStatusTimesCoxDetail(mat)
  y.phase3               <- tmp$status
  time.phase3            <- tmp$time
  time2.phase3           <- tmp[["time2", exact=TRUE]]
  tmp                    <- y.phase3 %in% 1
  if (is.null(time2.phase3)) {
    mod.detail.time.phase3 <- sort(unique(time.phase3[tmp])) 
  } else {
    mod.detail.time.phase3 <- sort(unique(time2.phase3[tmp]))
  }
  number.times.phase3    <- length(mod.detail.time.phase3)
  obstimes.phase3        <- obstimes.phase2[subset.phase3]  
  if (!is.null(obstimes0.phase2)) {
    obstimes0.phase3     <- obstimes0.phase2[subset.phase3] 
  } else {
    obstimes0.phase3     <- NULL
  }
  rm(time.phase3, tmp, mod.detail.time.phase3, mat, time2.phase3)
  gc()
  
  # obstimes.phase3 must be used to get a column of the matrix riskmat.phase3

  exp.X.weighted    <- weights*exp(X %*% beta.hat)

  mod.detail.time   <- mod.detail$time  # Event times
  number.times      <- length(mod.detail.time)
  tmp               <- getStatusTimesCoxDetail(mod.detail$y)
  y                 <- tmp$status
  time              <- tmp$time
  time2             <- tmp[["time2", exact=TRUE]]
  imat              <- mod.detail$imat
  if (!is.null(dim(imat))) {
    infomat         <- apply(X=imat, MARGIN=c(1,2), FUN=sum)
  } else {
    infomat         <- sum(imat)
  }
  infomat.inv       <- solve(infomat) 

  if (is.null(time2)) {
    obstimes        <- cc_getTimeId(time, mod.detail.time)
    obstimes0       <- NULL
  } else {
    obstimes        <- cc_getTimeId(time2, mod.detail.time)
    obstimes0       <- cc_getTime0Id(time, mod.detail.time)
  }
  rm(imat, infomat, tmp, time, time2)
  gc()

  if (is.null(Tau1)) Tau1 <- floor(min(mod.detail.time))
  if (is.null(Tau2)) Tau2 <- floor(max(mod.detail.time))
  eventtimes.phase2    <- mod.detail.time.phase2
  Tau1Tau2.times       <- which((Tau1 < eventtimes.phase2) & (eventtimes.phase2 <= Tau2)) 

  if (!length(Tau1Tau2.times)) stop("ERROR: no event times in (tau1, tau2], adjust tau1 and/or tau2")
  rm(eventtimes.phase2, mod.detail.time.phase2)
  gc()

  S0t             <- cc_getS0t(obstimes, obstimes0, number.times, exp.X.weighted)
  if (any(S0t == 0)) stop("ERROR: S0t = 0")

  S1t             <- cc_getS1t(obstimes, obstimes0, exp.X.weighted, X, number.times)
  S0t.casestimes  <- cc_getS0t(obstimes.phase3, obstimes0.phase3, number.times.phase2, 
                               exp.X.weighted)
  tmp             <- S0t.casestimes == 0
  if (any(tmp)) {
    # Check if any of them are included in the Tau1Tau2.times. If yes, then remove them.
    tmp0 <- S0t.casestimes[Tau1Tau2.times] == 0
    if (any(tmp0)) {
      warning("Restricting event times in (tau1, tau2] to exclude event times with no risk set")
      Tau1Tau2.times <- Tau1Tau2.times[!tmp0]
      if (!length(Tau1Tau2.times)) stop("ERROR: no event times in (tau1, tau2], adjust tau1 and/or tau2")
    }
    
    # Set to a small value to prevent errors later. NOTE: these small values will not be used
    #   in the program, just to prevent errors in a .C call
    S0t.casestimes[tmp] <- 1.0e-6
  }
  S1t.casestimes  <- cc_getS1t(obstimes.phase3, obstimes0.phase3, exp.X.weighted, X, number.times.phase2)

  score.beta           <- cc_getBetaScore(X, obstimes, obstimes0, number.times, y, 
                                          weights, exp.X.weighted, S1t, S0t)
  lambda0.t.hat        <- cc_getdNtColSums(obstimes.phase2, number.times.phase2, y.phase2)/S0t.casestimes
  Lambda0.Tau1Tau2.hat <- sum(lambda0.t.hat[Tau1Tau2.times])


  if (estimated.weights) {
    # --------------------------------------------------------------------------
    # Computation of the influences for gamma ----------------------------------
    J3              <- ncol(B.phase2)
    B.phase3        <- B.phase2[indiv.phase3, , drop=FALSE]
    weights.p3.est  <- calibration(B.phase3, 1, colSums(B.phase2), eta0=est.op[["eta0", exact=TRUE]], 
                         niter.max=est.op[["niter.max", exact=TRUE]], 
                         epsilon.stop=est.op[["epsilon.stop", exact=TRUE]])$calibrated.weights
    B.p3weights          <- sweep(B.phase3, 1, weights.p3.est, "*")
    sum.BB.p3weights     <- cc_getSumAAwgt(B.phase3, B.p3weights)
    sum.BB.p3weights.inv <- solve(sum.BB.p3weights)
    
    infl2.gamma                <- B.phase2 %*% sum.BB.p3weights.inv
    infl3.gamma                <- 0 * B.phase2
    infl3.gamma[indiv.phase3,] <- - B.p3weights %*% sum.BB.p3weights.inv
    infl.gamma                 <- infl2.gamma + infl3.gamma 

    rm(weights.p3.est, B.p3weights, sum.BB.p3weights, sum.BB.p3weights.inv)
    gc()
   
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    tmp <- cc_get_drond_U_eta(X, B.phase3, obstimes, obstimes0, number.times, y, 
                              S1t, S0t, exp.X.weighted, weights)
    drond.U.gamma              <- tmp$drond.U.eta

    tmp                        <- t(drond.U.gamma) %*% infomat.inv
    infl2.beta                 <- infl2.gamma %*% tmp
    infl3.beta                 <- infl3.gamma %*% tmp
    infl3.beta[indiv.phase3, ] <- infl3.beta[indiv.phase3,] + score.beta %*% infomat.inv
    infl.beta                  <- infl2.beta + infl3.beta
    rm(tmp, score.beta, drond.U.gamma)
    gc()    

    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    drond.S0t.gamma.casestimes <- cc_getS0GammaCasetimes(B.phase3, obstimes.phase3, obstimes0.phase3,
                                    number.times.phase2, exp.X.weighted)
    infl2.Lambda0.Tau1Tau2 <- cc_infl2_lambda0t(obstimes.phase2, y.phase2, infl2.beta, S1t.casestimes,
                      infl2.gamma, drond.S0t.gamma.casestimes, S0t.casestimes,
                      lambda0.t.hat, Tau1Tau2.times)
    infl3.Lambda0.Tau1Tau2 <- cc_infl3_lambda0t(obstimes.phase3, obstimes0.phase3, exp.X.weighted, 
                      infl3.beta, S1t.casestimes,
                      infl3.gamma, drond.S0t.gamma.casestimes, S0t.casestimes,
                      lambda0.t.hat, Tau1Tau2.times, subset.phase3)

    infl.Lambda0.Tau1Tau2 <- infl2.Lambda0.Tau1Tau2 + infl3.Lambda0.Tau1Tau2
    
    return(list(infl.beta = infl.beta, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl2.beta = infl2.beta, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl3.beta = infl3.beta, 
                infl3.Lambda0.Tau1Tau2 = as.matrix(infl3.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  } else {
    # --------------------------------------------------------------------------
    # If missingness weights known ---------------------------------------------
    
    # --------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta --------------
    
    infl2.beta                <- matrix(0, nrow=length(indiv.phase2), ncol=p)
    infl3.beta                <- infl2.beta
    rownames(infl3.beta)      <- indiv.phase2
    infl3.beta[indiv.phase3,] <- score.beta %*% infomat.inv
    infl.beta                 <- infl2.beta + infl3.beta
    rm(score.beta)
    gc()
    # --------------------------------------------------------------------------
    # Computation of the influences for the baseline hazards at each unique 
    # event time and for the cumulative baseline hazard in time interval 
    # [Tau1, Tau2] -------------------------------------------------------------
    
    infl2.Lambda0.Tau1Tau2 <- cc_infl2_lambda0t_noEst(obstimes.phase2, y.phase2, 
                                            S0t.casestimes, Tau1Tau2.times) 
    infl3.Lambda0.Tau1Tau2 <- cc_infl3_lambda0t_noEst(obstimes.phase3, obstimes0.phase3, exp.X.weighted, 
                            infl3.beta, S1t.casestimes, S0t.casestimes, lambda0.t.hat, 
                            Tau1Tau2.times, subset.phase3)
    infl.Lambda0.Tau1Tau2  <- infl2.Lambda0.Tau1Tau2 + infl3.Lambda0.Tau1Tau2
     
    return(list(infl.beta = infl.beta, 
                infl.Lambda0.Tau1Tau2 = as.matrix(infl.Lambda0.Tau1Tau2), 
                infl2.beta = infl2.beta, 
                infl2.Lambda0.Tau1Tau2 = as.matrix(infl2.Lambda0.Tau1Tau2), 
                infl3.beta = infl3.beta, 
                infl3.Lambda0.Tau1Tau2 = as.matrix(infl3.Lambda0.Tau1Tau2), 
                beta.hat = beta.hat, 
                lambda0.t.hat = lambda0.t.hat, 
                Lambda0.Tau1Tau2.hat = Lambda0.Tau1Tau2.hat))
  }

}



