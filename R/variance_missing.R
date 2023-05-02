cc_var_miss <- function(data.p2, obj, infl2, infl3, estimated.weights=TRUE) {

  # In the variance.missingdata function, 
  #-	weights should be a vector with containing the design weights (= weights.phase2*weights.phase3) for all the case-cohort members (=phase-three sample). 
  #     As you mentioned, weights.phase3 may have been specified by the user, or estimated by the estimation.weights.phase3 function.
  #-	weights.phase2 should be a vector containing the design weights (= weights.phase2*weights.phase3) for all the individuals in the phase-two sample.
  #-	weights.p2.phase2 should be a vector containing the phase-two design weights (= weights.phase2) for all the individuals in the phase-two sample

  n                 <- obj$nrow.cohort     
  status            <- data.p2[, obj$status, drop=TRUE]
  weights.p2.phase2 <- data.p2[, obj$weights.phase2, drop=TRUE]
  weights           <- weights.p2.phase2*data.p2[, obj$weights.phase3, drop=TRUE]
  estWgtFlag        <- estimated.weights
  stratFlag         <- obj$stratFlag
  strat3Flag        <- obj$strat3Flag
  W                 <- NULL
  W3                <- NULL
  if (stratFlag)  W  <- data.p2[, obj$strata, drop=TRUE]
  if (strat3Flag) {
    W3 <- data.p2[, obj$strata.phase3, drop=TRUE]
  } else {
    W3 <- rep(1, nrow(data.p2))
  }



  if (estWgtFlag) {
    # weights not used, phase3.01 not used, use phase 2 sample,
    #  use weights.phase2, weights.p2.phase2
    weights.phase2 <- weights
    weights        <- NULL
    phase3.01      <- NULL
  } else {
    # use phase 3 sample, weights, weights.phase2, weights.p2.phase2
    # weights.phase2 and weights.p2.phase2 are for phase 2 subs
    # weights are for phase 3 subs
    phase3.01 <- data.p2[, obj$phase3, drop=TRUE]
    if (!is.null(W)) W <- W[phase3.01]
    W3             <- W3[phase3.01]
    weights.phase2 <- weights
    weights        <- weights[phase3.01]
    status         <- status[phase3.01]
  }

  tmp    <- obj$strata.counts
  if (!length(tmp)) stop("INTERNAL CODING ERROR 1")
  s.n    <- tmp$cnts.cohort
  s.m    <- tmp$cnts.casecohort
  s.m.n  <- s.m/s.n
  valVec <- 1 - s.m.n*(s.n - 1)/(s.m - 1)
  tmp    <- s.m <= 1
  if (any(tmp)) valVec[tmp] <- 0

  ret <- cc_variance.missingdata(n, W, W3, status, phase3.01, valVec,
                                  weights, weights.phase2, weights.p2.phase2, 
                                  infl2, infl3, stratified.p2=stratFlag,
                                  estimated.weights=estWgtFlag, DEBUG=obj$DEBUG) 
  ret
}

cc_variance.missingdata <- function(n, W, W3, status, phase3.01, valVec,
                                  weights, weights.phase2, weights.p2.phase2, 
                                  infl2, infl3, stratified.p2=FALSE,
                                  estimated.weights=FALSE, DEBUG=0) {

  if (estimated.weights == FALSE) {
    if (stratified.p2 == TRUE){
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))

      phase23.var <-  cc_phase23VarEstF(W, W3, status, weights, valVec, 
                                 infl3[phase3.01, , drop=FALSE], stratified.p2, DEBUG=DEBUG)
      var         <- diag(superpop.var) + phase23.var
      
    } else {
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 /  weights.p2.phase2,
                                                    infl3) +
                                      crossprod(infl3 / sqrt(weights.phase2)))
      
      phase23.var <-  cc_phase23VarEstF(NULL, W3, status, weights, valVec, 
                                        infl3[phase3.01, , drop=FALSE], stratified.p2, DEBUG=DEBUG)
      var         <- diag(superpop.var) + phase23.var
    }
  } else {
    
    weights.p3.phase2 <- weights.phase2 / weights.p2.phase2
    
    if (stratified.p2 == TRUE) {
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, 
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
      
      phase23.var <-  cc_phase23VarEstT(W, W3, status,
                                    weights.phase2, weights.p2.phase2, weights.p3.phase2, 
                                    valVec, infl2, infl3, stratified.p2, DEBUG=DEBUG)
      var         <- diag(superpop.var) + phase23.var
    } else {
    
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, 
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
   
      phase23.var <-  cc_phase23VarEstT(NULL, W3, status,
                          weights.phase2, weights.p2.phase2, weights.p3.phase2, 
                          valVec, infl2, infl3, stratified.p2, DEBUG=DEBUG)
      var        <- diag(superpop.var) + phase23.var
    }
  }
  var
}







## -----------------------------------------------------------------------------
## Function: variance.missingdata()
## -----------------------------------------------------------------------------
## Description: This function returns the variance estimate that influence-based
##              and follows the complete variance decomposition, for 
##              parameters such as log-relative hazard, cumulative baseline 
##              hazard or covariate specific pure-risk, when covariate 
##              information is missing for individuals in phase-two sample
## -----------------------------------------------------------------------------
## Arguments:
##
## n                  number of individuals in the whole cohort
##
## casecohort         if stratified = TRUE, data frame with W (phase-two 
##                    strata), strata.m (number of sampled individuals in the 
##                    strata in the second phase), strata.n (strata sizes in the 
##                    cohort), for each individual in the case cohort. If 
##                    stratified = FALSE, data frame with m (number of sampled 
##                    individuals) and n (cohort size), for each individual in 
##                    the case cohort
##                    
## casecohort.phase2  if stratified = TRUE, data frame with W (phase-two 
##                    strata), strata.m (number of sampled individuals in the 
##                    strata in the second phase), strata.n (strata sizes in the 
##                    cohort) and phase3 (phase-three sampling indicator), for 
##                    each individual in the phase-two sample. If stratified = 
##                    FALSE, data frame with m (number of sampled individuals in
##                    the second phase), n (cohort size) and unstrat.phase3 
##                    (phase-three sampling indicator), for each individual in 
##                    the phase two sample
##
## weights            vector with design weights for each individual in the case
##                    cohort data
##
## weights.phase2     vector with design weights for each individual in the 
##                    phase two sample
##
## weights.p2.phase2  vector with phase-two design weights for each individual 
##                    in the phase two sample
##
## infl2              phase-two influences on a parameter such as log-relative 
##                    hazard, cumulative baseline hazard or covariate specific 
##                    pure-risk  
##
## infl3              phase-three influences on a parameter such as log-relative 
##                    hazard, cumulative baseline hazard or covariate specific 
##                    pure-risk  
##
## stratified.p2      was the second phase of sampling stratified on W? Default 
##                    is FALSE
##
## estimated.weights  were the phase-three weights estimated? Default is FALSE
## -----------------------------------------------------------------------------


variance.missingdata <- function(n, casecohort, casecohort.phase2, 
                                  weights, weights.phase2, weights.p2.phase2, 
                                  infl2, infl3, stratified.p2 = NULL,
                                  estimated.weights = NULL) {

  if (is.null(stratified.p2)) {
    stratified.p2 <- FALSE
  }
  if (is.null(estimated.weights)) {
    estimated.weights <- FALSE
  }

  if (estimated.weights == FALSE) {
    if (stratified.p2 == TRUE){
      casecohort.W    <- casecohort$W
      casecohort.W3   <- casecohort$W3
      casecohort.n    <- casecohort$strata.n
      casecohort.m    <- casecohort$strata.m
      casecohort.m.n  <- casecohort.m / casecohort.n
      casecohort.n.m  <- casecohort.n / casecohort.m
      casecohort.status <- casecohort$status      

      product.covar.weight <- function (indiv) { 
        return((casecohort.W == casecohort.W[indiv]) * 
                 (casecohort.W3 == casecohort.W3[indiv]) * 
                 (casecohort.status[indiv] == 0) * 
                 (1 - casecohort.m.n * (casecohort.n - 1) / (casecohort.m - 1) ))
      }
      prod.covar.weight       <- do.call(rbind, lapply(1:nrow(casecohort), 
                                                       product.covar.weight))
      diag(prod.covar.weight) <- (1 - 1 / weights)

      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2,
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))


      phase23.var   <- with(casecohort.phase2, 
                            t(infl3[phase3 == 1,]) %*% prod.covar.weight %*% 
                              infl3[phase3 == 1,])

      var           <- diag(superpop.var + phase23.var)
      
    } else {
      casecohort.W3   <- casecohort$W3
      casecohort.n    <- casecohort$n
      casecohort.m    <- casecohort$m
      casecohort.m.n  <- casecohort.m / casecohort.n
      casecohort.n.m  <- casecohort.n / casecohort.m
      casecohort.status <- casecohort$status      

      product.covar.weight <- function (indiv) { 
        return((casecohort.W3 == casecohort.W3[indiv]) * 
                 (casecohort.status[indiv] == 0) * 
                 (1 - casecohort.m.n * (casecohort.n - 1) / (casecohort.m - 1) ))
      }
      prod.covar.weight       <- do.call(rbind, lapply(1:nrow(casecohort), 
                                                       product.covar.weight))
      diag(prod.covar.weight) <- (1 - 1 / weights)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 /  weights.p2.phase2,
                                                    infl3) +
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- with(casecohort.phase2, 
                            t(infl3[unstrat.phase3 == 1,]) %*% 
                              prod.covar.weight %*% infl3[unstrat.phase3 == 1,])
      var           <- diag(superpop.var + phase23.var)
    }
  } else {
    
    weights.p3.phase2 <- weights.phase2 / weights.p2.phase2
    infl <- infl2 + infl3
    
    if (stratified.p2 == TRUE) {
      phase2.W   <- casecohort.phase2$W
      phase2.W3   <- casecohort.phase2$W3
      phase2.n    <- casecohort.phase2$strata.n
      phase2.m    <- casecohort.phase2$strata.m
      phase2.m.n  <- phase2.m / phase2.n
      phase2.n.m  <- phase2.n / phase2.m
      phase2.status <- casecohort.phase2$status

      product.covar.weight <- function (indiv) {
        return((phase2.W == phase2.W[indiv]) * 
                 (phase2.W3 == phase2.W3[indiv]) * 
                 (phase2.status[indiv] == 0) * 
                 (1 - phase2.m.n * (phase2.n - 1) / (phase2.m - 1) ))
      }
      prod.covar.weight        <- do.call(rbind, 
                                           lapply(1:nrow(casecohort.phase2), 
                                                  product.covar.weight))
      diag(prod.covar.weight)  <- (1 - 1 / weights.phase2)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, 
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- t(infl) %*% prod.covar.weight %*% infl - t(infl2) %*% 
        diag(c(1 / weights.p2.phase2 * (1 - 1 / weights.p3.phase2))) %*% 
        (infl2 + 2 * infl3)
      var           <- diag(superpop.var + phase23.var)
    } else {
      phase2.W3   <- casecohort.phase2$W3
      phase2.n    <- casecohort.phase2$n
      phase2.m    <- casecohort.phase2$m
      phase2.m.n  <- phase2.m / phase2.n
      phase2.n.m  <- phase2.n / phase2.m
      phase2.status <- casecohort.phase2$status

      product.covar.weight <- function (indiv) {
        return((phase2.W3 == phase2.W3[indiv]) * 
                 (phase2.status[indiv] == 0) * 
                 (1 - phase2.m.n * (phase2.n - 1) / (phase2.m - 1) ))
      }
      prod.covar.weight        <- do.call(rbind, 
                                           lapply(1:nrow(casecohort.phase2), 
                                                  product.covar.weight))
      diag(prod.covar.weight)  <- (1 - 1 / weights.phase2)
      
      superpop.var  <- n / (n-1) * (crossprod(infl2 / sqrt(weights.p2.phase2)) +
                                      2 * crossprod(infl2 / weights.p2.phase2, 
                                                    infl3) + 
                                      crossprod(infl3 / sqrt(weights.phase2)))
      phase23.var   <- t(infl) %*% prod.covar.weight %*% infl - t(infl2) %*% 
        diag(c(1 / weights.p2.phase2 * (1 - 1 / weights.p3.phase2))) %*% 
        (infl2 + 2 * infl3)
      var           <- diag(superpop.var + phase23.var)
    }
  }
  return(variance = var)
}

