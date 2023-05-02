
cc_getRiskMatCol <- function(col, timeInd, time0Ind) {

  if (is.null(time0Ind)) time0Ind <- -9
  n   <- length(timeInd)
  ret <- rep(-9999.9, n)
  tmp <- .C("C_getRiskMatCol", as.integer(col-1), as.integer(timeInd-1), as.integer(time0Ind-1),
            as.integer(n), ret=as.numeric(ret))
  tmp$ret

}

cc_getRiskMatCol2 <- function(col, timeInd, time0Ind) {

  if (any(timeInd < time0Ind)) stop("INTERNAL CODING ERROR")
  n   <- length(timeInd)
  ret <- rep(-9999.9, n)
  tmp <- .C("C_getRiskMatCol2", as.integer(col-1), as.integer(timeInd-1), as.integer(time0Ind-1), 
             as.integer(n), ret=as.numeric(ret))
  tmp$ret

}


cc_getYexpXWgtCol <- function(col, timeInd, time0Ind, expXwgt) {

  if (!length(time0Ind)) time0Ind <- -9
  nsub <- length(timeInd)
  if (nsub != length(expXwgt)) stop("INTERNAL CODING ERROR 1")
  ret  <- rep(-9999.9, nsub) 
  tmp  <- .C("C_getYexpXWgtCol", as.integer(col), as.integer(timeInd-1), as.integer(time0Ind-1),
             as.integer(nsub), as.numeric(expXwgt), ret=as.numeric(ret))
  tmp$ret

}

cc_getS0t <- function(timeInd, time0Ind, ntimes, exp.X.weighted, DEBUG=0) {

  if (!length(time0Ind)) time0Ind <- -9
  nsub <- length(timeInd)
  ret  <- rep(0, ntimes)
  if (DEBUG) cat("Begin: C_getS0t\n")
  tmp  <- .C("C_getS0t", as.integer(timeInd-1), as.integer(time0Ind-1), as.numeric(exp.X.weighted), 
             as.integer(nsub), as.integer(ntimes), ret=as.numeric(ret),
             PACKAGE="CaseCohortCoxSurvival")
  if (DEBUG) cat("End: C_getS0t\n")
  ret  <- tmp$ret
  ret
}

cc_getS1t <- function(timeInd, time0Ind, expXwgt, X, ntimes, DEBUG=0) {

  if (!length(time0Ind)) time0Ind <- -9
  nsub <- length(timeInd)
  p    <- ncol(X)
  ret  <- rep(0, ntimes*p)
  if (DEBUG) cat("Begin: C_getS1t\n")
  tmp <- .C("C_getS1t", as.integer(timeInd-1), as.integer(time0Ind-1), as.numeric(expXwgt), 
            as.numeric(X), as.integer(nsub), as.integer(ntimes), as.integer(p),
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival")
  if (DEBUG) cat("End: C_getS1t\n")
  ret <- matrix(tmp$ret, nrow=ntimes, ncol=p, byrow=TRUE)
  ret
}

cc_getdNtColSums <- function(timeInd, ntimes, y, DEBUG=0) {

  nsub <- length(timeInd)
  ret  <- rep(0, ntimes)
  if (DEBUG) cat("Begin: C_getdNtColSums\n")
  tmp  <- .C("C_getdNtColSums", as.integer(timeInd-1), as.integer(nsub), 
              as.integer(ntimes), as.numeric(y), 
              ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival")
  if (DEBUG) cat("End: C_getdNtColSums\n")
  ret  <- tmp$ret
  ret

}

cc_getdNtWgtColSums <- function(timeInd, ntimes, y, wgt, DEBUG=0) {

  nsub <- length(timeInd)
  ret  <- rep(0, ntimes)
  if (DEBUG) cat("Begin: C_getdNtWgtColSums\n")
  tmp  <- .C("C_getdNtWgtColSums", as.integer(timeInd-1), as.integer(nsub), 
              as.integer(ntimes), as.numeric(y), as.numeric(wgt), 
              ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival")
  if (DEBUG) cat("End: C_getdNtWgtColSums\n")
  ret  <- tmp$ret
  ret

}

cc_getSumAAwgt <- function(Aphase2ByRow, AwgtByRow, DEBUG=0) {

  nsub <- nrow(Aphase2ByRow)
  q    <- ncol(Aphase2ByRow)
  ret  <- rep(0, q*q)
  if (DEBUG) cat("Begin: C_getSumAAwgt\n")
  tmp  <- .C("C_getSumAAwgt", as.numeric(t(Aphase2ByRow)), as.numeric(t(AwgtByRow)), 
             as.integer(nsub), as.integer(q), ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival")
  if (DEBUG) cat("End: C_getSumAAwgt\n")
  ret <- matrix(tmp$ret, nrow=q, ncol=q, byrow=TRUE)
  ret
}

cc_get_drond_U_eta <- function(XbyRow, AbyRow, timeInd, time0Ind, ntimes, y, S1ByRow, S0, expXwgt,
                               weights, DEBUG=0) {

  if (!length(time0Ind)) time0Ind <- -9
  n                <- nrow(AbyRow)
  q                <- ncol(AbyRow)
  p                <- ncol(XbyRow)
  ret_dS0tEtaByCol <- rep(0, ntimes*q)
  ret_dUtEtaByRow  <- rep(0, p*q)
  if (DEBUG) cat("Begin: C_get_drond_U_eta\n")
  
  tmp <- .C("C_get_drond_U_eta", as.integer(DEBUG), as.numeric(t(XbyRow)), as.numeric(t(AbyRow)), as.integer(p), as.integer(q), 
            as.integer(n), as.integer(ntimes), as.integer(timeInd-1), as.integer(time0Ind-1), 
            as.numeric(y), as.numeric(t(S1ByRow)),
            as.numeric(S0), as.numeric(expXwgt), as.numeric(weights), 
            ret_dS0tEtaByCol=as.numeric(ret_dS0tEtaByCol), 
            ret_dUtEtaByRow=as.numeric(ret_dUtEtaByRow), PACKAGE="CaseCohortCoxSurvival") 
  if (DEBUG) cat("End: C_get_drond_U_eta\n")

  ret_dS0tEtaByCol = matrix(tmp$ret_dS0tEtaByCol, nrow=q, ncol=ntimes, byrow=FALSE)
  ret_dUtEtaByRow  = matrix(tmp$ret_dUtEtaByRow, nrow=p, ncol=q, byrow=TRUE)

  list(drond.S0t.eta=ret_dS0tEtaByCol, drond.U.eta=ret_dUtEtaByRow)
}

cc_getBetaScore <- function(XbyRow, timeInd, time0Ind, ntimes, y, weights, 
                            expXwgt, S1t, S0t, DEBUG=0) {

  if (!length(time0Ind)) time0Ind <- -9
  nsub     <- length(timeInd)
  p        <- ncol(XbyRow)
  XeByRow  <- t(S1t/matrix(S0t, nrow=ntimes, ncol=p, byrow=FALSE))
  csumS0   <- cc_getdNtWgtColSums(timeInd, ntimes, y, weights)/S0t
  retByRow <- rep(0, nsub*p)
  if (DEBUG) cat("Begin: C_getBetaScore\n")
  tmp <- .C("C_getBetaScore", as.integer(nsub), as.integer(p), as.numeric(XeByRow), 
            as.numeric(t(XbyRow)), as.integer(timeInd-1), as.integer(time0Ind-1), 
            as.integer(ntimes), as.numeric(y), 
            as.numeric(weights), as.numeric(expXwgt), as.numeric(csumS0), 
            retByRow=as.numeric(retByRow), PACKAGE="CaseCohortCoxSurvival")
  if (DEBUG) cat("End: C_getBetaScore\n")

  retByRow <- matrix(tmp$retByRow, nrow=nsub, ncol=p, byrow=TRUE)
  retByRow
}

cc_infl_lambda0_tau12 <- function(infl1BetaByRow, S1tByRow, infl1EtaByRow, 
                          dS0tByCol, infl2BetaByRow, infl2EtaByRow,
                          lambda0, S0, timeInd, time0Ind, y, expXwgt, 
                          ntimes, Tau12Times, phase2, DEBUG=0) {

  if (!length(time0Ind)) time0Ind <- -9
  nsub        <- length(y)
  nsubA       <- nrow(infl1BetaByRow)
  p           <- ncol(S1tByRow)
  q           <- nrow(dS0tByCol)
  nTau12Times <- length(Tau12Times)
  ret1        <- rep(0, nsubA)
  ret2        <- rep(0, nsubA) # in C code, only nsub is needed
  if (DEBUG) cat("Begin: C_infl_lambda0_tau12\n")
  tmp  <- .C("C_infl_lambda0_tau12", as.numeric(t(infl1BetaByRow)), as.numeric(t(S1tByRow)), 
            as.numeric(t(infl1EtaByRow)), as.numeric(dS0tByCol), as.numeric(t(infl2BetaByRow)),
            as.numeric(t(infl2EtaByRow)), as.numeric(lambda0), as.numeric(S0), 
            as.integer(timeInd-1), as.integer(time0Ind-1), as.numeric(y), as.numeric(expXwgt), 
            as.integer(nsub), as.integer(nsubA), as.integer(ntimes), as.integer(p), as.integer(q),
            as.integer(Tau12Times-1), as.integer(nTau12Times), 
            ret1=as.numeric(ret1), ret2=as.numeric(ret2), PACKAGE="CaseCohortCoxSurvival") 
  if (DEBUG) cat("End: C_infl_lambda0_tau12\n")
  ret2[phase2] <- (tmp$ret2)[1:nsub]
  list(infl1.Lambda0.Tau1Tau2=tmp$ret1, infl2.Lambda0.Tau1Tau2=ret2)

}

cc_infl_lambda0_tau12_noCalib <- function(inflBetaByRow, S1tByRow, 
                          lambda0, S0, timeInd, time0Ind, y, expXwgt, 
                          ntimes, Tau12Times, DEBUG=0) {

  if (!length(time0Ind)) time0Ind <- -9
  nsub        <- length(y)
  p           <- ncol(S1tByRow)
  nTau12Times <- length(Tau12Times)
  ret         <- rep(0, nsub)
  if (DEBUG) cat("Begin: C_infl_lambda0_tau12_noCalib\n")
  tmp  <- .C("C_infl_lambda0_tau12_noCalib", as.numeric(t(inflBetaByRow)),
            as.numeric(t(S1tByRow)), as.numeric(lambda0), as.numeric(S0), 
            as.integer(timeInd-1), as.integer(time0Ind-1),
            as.numeric(y), as.numeric(expXwgt), 
            as.integer(nsub), as.integer(ntimes), as.integer(p),
            as.integer(Tau12Times-1), as.integer(nTau12Times), 
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival") 
  if (DEBUG) cat("End: C_infl_lambda0_tau12_noCalib\n")
  tmp$ret

}

cc_getPhase2Var <- function(subcohort, y, infl2ByCol, stratVec.cohort,
                     stratCnts.casecohort, stratCnts.cohort) {

  # strataVec.cohort must contain strata 0, 1, ...
  # vectors of counts must contain all strata
  # subcohort is the phase-2 sample indicators

  stratFlag  <- !is.null(stratVec.cohort)
  subset     <- (y %in% 0) & subcohort
  nsubset    <- sum(subset)
  p          <- ncol(infl2ByCol)
  s.m        <- stratCnts.casecohort
  s.n        <- stratCnts.cohort
  n          <- length(y)
  m          <- sum(s.m)

  if (stratFlag) {
    W   <- stratVec.cohort[subset]
  } else {
    W   <- -9999
  }
  s.m.n <- s.m/s.n
  s.n.m <- s.n/s.m
  tmp    <- !is.finite(s.n.m)
  if (any(tmp)) s.n.m[tmp] <- 1

  omegaConst <- s.n.m*(s.n - 1)/(s.m - 1)
  tmp        <- !is.finite(omegaConst)
  if (any(tmp)) omegaConst[tmp] <- 0
  covarConst <- s.m.n*(s.m - 1)/(s.n - 1) - s.m.n^2
  tmp        <- !is.finite(covarConst)
  if (any(tmp)) covarConst[tmp] <- 0
  omegaDiag  <- s.n.m
  covarDiag  <- s.m.n*(1 - s.m.n)
  ret        <- as.numeric(rep(0, p))

  # For calibrated = FALSE
  if (nrow(infl2ByCol) < length(subset)) subset <- y[subcohort] %in% 0 

  # Error check
  if (nrow(infl2ByCol) != length(subset)) stop("INTERNAL CODING ERROR")

  tmp <- .C("C_getPhase2Var", as.integer(stratFlag), as.integer(W), 
             as.numeric(omegaConst), as.numeric(covarConst), 
             as.numeric(omegaDiag), as.numeric(covarDiag), as.integer(nsubset), 
             as.numeric(infl2ByCol[subset, , drop=FALSE]), as.integer(p), 
             ret=ret, PACKAGE="CaseCohortCoxSurvival") 

  tmp$ret
}

cc_getS0GammaCasetimes <- function(B3ByCol, obstimes3, obstimes3_0, neventTimes2, expXwgt) {

  if (!length(obstimes3_0)) obstimes3_0 <- -9
  n3       <- nrow(B3ByCol)
  J3       <- ncol(B3ByCol)
  retByCol <- rep(-9999.9, J3*neventTimes2)
  if (n3 != length(expXwgt)) stop("INTERNAL CODING ERROR 1")
  if (n3 != length(obstimes3)) stop("INTERNAL CODING ERROR 2")

  tmp <- .C("C_getS0GammaCasetimes", as.numeric(B3ByCol), 
            as.integer(obstimes3-1), as.integer(obstimes3_0-1),
            as.integer(neventTimes2), as.numeric(expXwgt),
            as.integer(n3), as.integer(J3), retByCol=as.numeric(retByCol), 
            PACKAGE="CaseCohortCoxSurvival")
  retByCol <- matrix(tmp$retByCol, nrow=J3, ncol=neventTimes2, byrow=FALSE)
  retByCol
}

cc_infl2_lambda0t <- function(obstimes2, y2, infl2BetaByRow, S1tCaseByRow,
                      infl2GammaByRow, dS0tGammaCaseByCol, S0tCase,
                      lambda0, Tau12Times) {

  n2          <- length(obstimes2)
  p           <- ncol(infl2BetaByRow)
  J3          <- nrow(dS0tGammaCaseByCol)
  nTau12Times <- length(Tau12Times)

  if (n2 != length(y2)) stop("INTERNAL CODING ERROR 1")
  if (n2 != nrow(infl2BetaByRow)) stop("INTERNAL CODING ERROR 2")
  if (n2 != nrow(infl2GammaByRow)) stop("INTERNAL CODING ERROR 3")
  if (nrow(S1tCaseByRow) != length(S0tCase)) stop("INTERNAL CODING ERROR 4")
  if (ncol(S1tCaseByRow) != p) stop("INTERNAL CODING ERROR 5")
  if (ncol(dS0tGammaCaseByCol) != length(S0tCase)) stop("INTERNAL CODING ERROR 6")

  ret <- rep(-9999.9, n2)
  tmp <- .C("C_infl2_lambda0t", as.integer(obstimes2-1), as.numeric(y2), 
            as.numeric(t(infl2BetaByRow)), as.numeric(t(S1tCaseByRow)),
            as.numeric(t(infl2GammaByRow)), as.numeric(dS0tGammaCaseByCol), 
            as.numeric(S0tCase), as.numeric(lambda0), as.integer(Tau12Times-1),
            as.integer(nTau12Times), as.integer(n2), as.integer(p), as.integer(J3), 
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival") 
  tmp$ret
}

cc_infl3_lambda0t <- function(obstimes3, obstimes3_0, expXwgt, infl3BetaByRow, S1tCaseByRow,
                      infl3GammaByRow, dS0tGammaCaseByCol, S0tCase,
                      lambda0, Tau12Times, phase3) {

  # phase3 vector of indices for phase 3 subjects
  if (!length(obstimes3_0)) obstimes3_0 <- -9
  n2          <- nrow(infl3BetaByRow)
  n3          <- length(phase3)
  p           <- ncol(infl3BetaByRow)
  J3          <- nrow(dS0tGammaCaseByCol)
  nTau12Times <- length(Tau12Times)

  if (n3 != length(expXwgt)) stop("INTERNAL CODING ERROR 1")
  if (n2 != nrow(infl3GammaByRow)) stop("INTERNAL CODING ERROR 3")
  if (nrow(S1tCaseByRow) != length(S0tCase)) stop("INTERNAL CODING ERROR 4")
  if (ncol(S1tCaseByRow) != p) stop("INTERNAL CODING ERROR 5")
  if (ncol(dS0tGammaCaseByCol) != length(S0tCase)) stop("INTERNAL CODING ERROR 6")
  if (n3 != length(obstimes3)) stop("INTERNAL CODING ERROR 7")

  # Must be initialized to 0 !!!
  ret <- rep(0, n2) # Only n3 will be used in C code
  tmp <- .C("C_infl3_lambda0t", as.integer(obstimes3-1), as.integer(obstimes3_0-1), as.numeric(expXwgt), 
            as.numeric(t(infl3BetaByRow[phase3, , drop=FALSE])), as.numeric(t(S1tCaseByRow)),
            as.numeric(t(infl3GammaByRow[phase3, , drop=FALSE])), as.numeric(dS0tGammaCaseByCol), 
            as.numeric(S0tCase), as.numeric(lambda0), as.integer(Tau12Times-1),
            as.integer(nTau12Times), as.integer(n3), as.integer(p), as.integer(J3), 
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival") 
  ret[phase3] <- (tmp$ret)[1:n3]

  ret
}

cc_infl2_lambda0t_noEst <- function(obstimes2, y2, S0tCase, Tau12Times) {

  n2          <- length(obstimes2)
  nTau12Times <- length(Tau12Times)

  if (n2 != length(y2)) stop("INTERNAL CODING ERROR 1")

  ret <- rep(-9999.9, n2)
  tmp <- .C("C_infl2_lambda0t_noEst", as.integer(obstimes2-1), as.numeric(y2), 
            as.numeric(S0tCase), as.integer(Tau12Times-1), as.integer(nTau12Times),
            as.integer(n2), ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival")
  tmp$ret

}

cc_infl3_lambda0t_noEst <- function(obstimes3, obstimes3_0, expXwgt, infl3BetaByRow, 
                            S1tCaseByRow, S0tCase, lambda0, Tau12Times, phase3) {

  if (!length(obstimes3_0)) obstimes3_0 <- -9
  n2          <- nrow(infl3BetaByRow)
  n3          <- length(obstimes3)
  nTau12Times <- length(Tau12Times)
  p           <- ncol(infl3BetaByRow)

  if (n3 != length(expXwgt)) stop("INTERNAL CODING ERROR 1")
  if (n3 != length(phase3)) stop("INTERNAL CODING ERROR 2")
  if (p != ncol(S1tCaseByRow)) stop("INTERNAL CODING ERROR 3")

  ret <- rep(0, n2)  # only first n3 will be used in c code
  tmp <- .C("C_infl3_lambda0t_noEst", as.integer(obstimes3-1), as.integer(obstimes3_0-1), as.numeric(expXwgt), 
            as.numeric(t(infl3BetaByRow[phase3, , drop=FALSE])), as.numeric(t(S1tCaseByRow)), 
            as.numeric(S0tCase), as.numeric(lambda0), as.integer(Tau12Times-1), 
            as.integer(nTau12Times), as.integer(n3), as.integer(p), 
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival")
  ret[phase3] <- (tmp$ret)[1:n3]
  ret
}

cc_phase23VarEstF <- function(W, W3, y, weights, valVec, infl3ByCol, stratFlag, DEBUG=0) {

  # Vectors passed in are for phase2 data

  n <- length(y)
  p <- ncol(infl3ByCol)

  nW <- length(W)
  if (nW && (n != nW)) stop("INTERNAL CODING ERROR 1")
  if (n != length(W3)) stop("INTERNAL CODING ERROR 2")
  if (n != length(weights)) stop("INTERNAL CODING ERROR 3")
  if (n != nrow(infl3ByCol)) stop("INTERNAL CODING ERROR 6")

  if (!nW) W <- -9999
  transWgt <- 1 - 1/weights
  ret <- rep(-9999.0, p)

  tmp <- .C("C_phase23VarEstF", as.integer(DEBUG), as.integer(W), as.integer(W3), as.integer(y), 
            as.integer(n), as.numeric(transWgt), as.numeric(valVec), 
            as.numeric(infl3ByCol), as.integer(p), as.integer(stratFlag), 
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival") 
  tmp$ret

}


cc_phase23VarEstT <- function(W, W3, y, weights.phase2, weights.p2.phase2, weights.p3.phase2, 
                               valVec, infl2ByCol, infl3ByCol, stratFlag, DEBUG=0) {

  # Vectors passed in are for phase2 data

  n <- length(y)
  p <- ncol(infl2ByCol)

  nW <- length(W)
  if (nW && (n != nW)) stop("INTERNAL CODING ERROR 1")
  if (n != length(W3)) stop("INTERNAL CODING ERROR 2")
  if (n != length(weights.p2.phase2)) stop("INTERNAL CODING ERROR 3")
  if (n != length(weights.p3.phase2)) stop("INTERNAL CODING ERROR 4")
  if (n != nrow(infl2ByCol)) stop("INTERNAL CODING ERROR 5")
  if (n != nrow(infl3ByCol)) stop("INTERNAL CODING ERROR 6")
  if (p != ncol(infl3ByCol)) stop("INTERNAL CODING ERROR 7")
  if (n != length(weights.phase2)) stop("INTERNAL CODING ERROR 8")

  if (!nW)  W <- -9999
  transWgt  <- 1 - 1/weights.phase2
  transWgt2 <- 1/weights.p2.phase2*(1 - 1/weights.p3.phase2) 
  ret <- rep(-9999.0, p)
  tmp <- .C("C_phase23VarEstT", as.integer(DEBUG), as.integer(W), as.integer(W3), as.integer(y), 
            as.integer(n), as.numeric(transWgt), as.numeric(transWgt2), 
            as.numeric(valVec), as.numeric(infl2ByCol), as.numeric(infl3ByCol), 
            as.integer(p), as.integer(stratFlag), 
            ret=as.numeric(ret), PACKAGE="CaseCohortCoxSurvival") 
  tmp$ret
}

cc_prodCovWgtStratT_row <- function(W, W3, y, row, valVec, transWgt) {

  n   <- length(y)
  ret <- rep(-9999999.0, n) 
  tmp <- .C("C_prodCovWgtStratT_row", as.integer(W), as.integer(W3), 
            as.integer(y), as.integer(row-1), as.integer(n), 
            as.numeric(valVec), as.numeric(transWgt), 
            ret=as.numeric(ret))
  tmp$ret

}

cc_prodCovWgtStratF_row <- function(W3, y, row, valVec, transWgt) {

  n   <- length(y)
  ret <- rep(-9999999.0, n) 
  tmp <- .C("C_prodCovWgtStratF_row", as.integer(W3), 
            as.integer(y), as.integer(row-1), as.integer(n), 
            as.numeric(valVec), as.numeric(transWgt), 
            ret=as.numeric(ret))
  tmp$ret

}
