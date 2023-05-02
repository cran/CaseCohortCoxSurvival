getCalibratedWgts <- function(data, obj, A) {

  if (!obj$calibrated) return(NULL)
  if (is.null(A)) stop("INTERNAL CODING ERROR")
  wgts <- cc_calibration(data, obj, A)  

  ret <- cc_addCalibWgts(wgts, data, obj)
  ret
}

cc_calibration <- function(data, obj, A) {

  if (obj$print) cat("Computing calibrated weights\n")
  if (nrow(A) != nrow(data)) stop("INTERNAL CODING ERROR")
  subset   <- getSubsetVec(data, obj) 
  weights  <- data[subset, obj$weights.phase2, drop=TRUE]
  A.phase2 <- A[subset, , drop=FALSE]
  total    <- colSums(A)
  op       <- obj$weights.op
  eta0     <- rep(0, ncol(A))
  
  tmp  <- calibration(A.phase2, weights, total, eta0=eta0, 
                       niter.max=op$niter.max, epsilon.stop=op$epsilon.stop)
  wgts <- tmp$calibrated.weights

  wgts
}

cc_addCalibWgts <- function(wgts, data, obj) {

  # Add weights to data
  wgtv                   <- getRandomVariable("wgt.calib_")
  data[, wgtv]           <- NA
  subset                 <- data[, obj$phase2, drop=TRUE]
  data[subset, wgtv]     <- wgts
  obj$weights.calibrated <- wgtv
  list(data=data, obj.list=obj)
}

## -----------------------------------------------------------------------------
## Function: calibration()
## -----------------------------------------------------------------------------
## Description: This function calibrates the design weights using the raking 
##              procedure. It matches the weighted totals of the auxiliary 
##              variables in the stratified case-cohort (with calibrated 
##              weights), to the un-weighted auxiliary variables totals in the 
##              whole cohort. The Newton Raphson method is used to solve the 
##              optimization problem
## -----------------------------------------------------------------------------
## Arguments:
##
##  A.phase2        matrix with the values of the auxiliary variables to be used 
##                  for the calibration of the weights in the stratified 
##                  case-cohort (phase-two data)
##
##  design.weights  design weights to be calibrated
##
##  total           un-weighted auxiliary variables totals in the whole cohort
##
##  eta0            vector of initial/possible values for eta, to be used as 
##                  seed in the iterative procedure. Default is (0,...,0)
##
##  niter.max       maximum number of iterations for the Newton Raphson method. 
##                  Default is 10^4 iterations
##
##  epsilon.stop    threshold for the difference between the estimated weighted 
##                  total and the total in the whole cohort. If this difference 
##                  is less than the value of epsilon.stop, no more iterations 
##                  will be performed. Default is 10^(-10)
## -----------------------------------------------------------------------------

calibration <- function (A.phase2, design.weights, total, eta0 = NULL, 
                         niter.max = NULL, epsilon.stop = NULL) {  
  A.phase2            <- as.matrix(A.phase2)
  
  n                   <- nrow(A.phase2)
  q                   <- ncol(A.phase2)
  
  if ((is.null(eta0)) || (length(eta0) != q)) {
    eta0          <- rep(0, q)
  }
  if ((is.null(epsilon.stop)) || (epsilon.stop <= 0)) {
    epsilon.stop  <- 10^(-10)
  }
  if ((is.null(niter.max)) || (niter.max <= 0)) {
    niter.max     <- 10^4
  }
  
  conv <- FALSE
  for (niter in 1:niter.max) {
    calibrated.weights  <- design.weights * exp(A.phase2 %*% eta0)
    A.calweighted       <- sweep(A.phase2, 1, calibrated.weights, "*")
    f                   <- colSums(A.calweighted) - total
    
    epsilon             <- max(abs(f))
    if (epsilon < epsilon.stop){
      conv <- TRUE
      break
    }
   
    #AA.calweighted      <- array(NA, dim = c(q, q, n))
    #for (i in 1:n) {
    #  AA.calweighted[,, i] <- tcrossprod(A.phase2[i,], A.calweighted[i,])
    #  }
    #drond.f.eta0        <- apply(X = AA.calweighted, MARGIN = c(1,2), FUN = sum)
    drond.f.eta0        <- cc_getSumAAwgt(A.phase2, A.calweighted)

    inv                 <- try(solve(drond.f.eta0))
    if ("try-error" %in% class(inv)) stop("ERROR: singular matrix obtained in calibration algorithm")
    eta0                <- eta0 - c(f %*% inv)
  }
  if (!conv) stop("ERROR: calibration algorithm did not converge")
  calibrated.weights    <- design.weights * exp(A.phase2 %*% eta0)
  if (any(!is.finite(calibrated.weights))) stop("ERROR: non-finite calibrated weights obtained")
  A.calweighted         <- sweep(A.phase2, 1, calibrated.weights, "*")
  
  return(list(eta.hat = eta0, calibrated.weights = calibrated.weights, 
              estimated.total = colSums(A.calweighted)))
}

## -----------------------------------------------------------------------------
## Function: calibration()
## -----------------------------------------------------------------------------
## Description: This function calibrates the design weights using the raking 
##              procedure. It matches the weighted totals of the auxiliary 
##              variables in the stratified case-cohort (with calibrated 
##              weights), to the un-weighted auxiliary variables totals in the 
##              whole cohort. The Newton Raphson method is used to solve the 
##              optimization problem
## -----------------------------------------------------------------------------
## Arguments:
##
##  A.phase2        matrix with the values of the auxiliary variables to be used 
##                  for the calibration of the weights in the stratified 
##                  case-cohort (phase-two data)
##
##  design.weights  design weights to be calibrated
##
##  total           un-weighted auxiliary variables totals in the whole cohort
##
##  eta0            vector of initial/possible values for eta, to be used as 
##                  seed in the iterative procedure. Default is (0,...,0)
##
##  niter.max       maximum number of iterations for the Newton Raphson method. 
##                  Default is 10^4 iterations
##
##  epsilon.stop    threshold for the difference between the estimated weighted 
##                  total and the total in the whole cohort. If this difference 
##                  is less than the value of epsilon.stop, no more iterations 
##                  will be performed. Default is 10^(-10)
## -----------------------------------------------------------------------------

calibration.orig <- function (A.phase2, design.weights, total, eta0 = NULL, 
                         niter.max = NULL, epsilon.stop = NULL) {
  
  A.phase2            <- as.matrix(A.phase2)
  
  n                   <- nrow(A.phase2)
  q                   <- ncol(A.phase2)
  
  if ((is.null(eta0)) || (length(eta0) != q)) {
    eta0          <- rep(0, q)
    message("eta0 has not been filled out. Default is 0. ")
  }
  if ((is.null(epsilon.stop)) || (epsilon.stop <= 0)) {
    epsilon.stop  <- 10^(-10)
    message("epsilon.stop has not been filled out. Default is 10^(-10). ")
  }
  if ((is.null(niter.max)) || (niter.max <= 0)) {
    niter.max     <- 10^4
    message("niter.max has not been filled out. Default is 10^4 iterations.")
  }
  
  for (niter in 1:niter.max) {
    calibrated.weights  <- design.weights * exp(A.phase2 %*% eta0)
    A.calweighted       <- sweep(A.phase2, 1, calibrated.weights, "*")
    f                   <- colSums(A.calweighted) - total
    
    epsilon             <- max(abs(f))
    if(epsilon < epsilon.stop){
      break
      }
    
    AA.calweighted      <- array(NA, dim = c(q, q, n))
    for (i in 1:n) {
      AA.calweighted[,, i] <- tcrossprod(A.phase2[i,], A.calweighted[i,])
      }
    drond.f.eta0        <- apply(X = AA.calweighted, MARGIN = c(1,2), FUN = sum)
    eta0                <- eta0 - c(f %*% solve(drond.f.eta0))
  }
  
  calibrated.weights    <- design.weights * exp(A.phase2 %*% eta0)
  A.calweighted         <- sweep(A.phase2, 1, calibrated.weights, "*")
  
  return(list(eta.hat = eta0, calibrated.weights = calibrated.weights, 
              estimated.total = colSums(A.calweighted)))
}

