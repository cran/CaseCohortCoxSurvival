print.casecohortcoxsurv <- function(x, ...) {

  dig  <- 4
  xnms <- names(x)
  beta <- x$beta
  lam  <- x$Lambda0
  col1 <- c(beta, lam)
  bnms <- names(beta)
  lnm  <- "Lambda0"
  if (lnm %in% bnms) lnm <- paste0(".", lnm, ".")
  rnms <- c(bnms, lnm)
  if ("beta.var" %in% xnms) {
    ncol <- 3
    p3   <- 0
  } else {
    ncol <- 5
    p3   <- 1
  }
  mat           <- matrix(data=NA, nrow=length(col1), ncol=ncol)
  mat[, 1]      <- col1
  ii            <- 1
  rownames(mat) <- rnms
  cnms          <- rep("", ncol)
  cnms[1]       <- "Parm"
  if (!p3) {
    ii        <- ii + 1
    lv        <- x[["Lambda0.var", exact=TRUE]]
    if (is.null(lv)) lv <- NA
    mat[, ii] <- c(x$beta.var, lv)
    cnms[ii]  <- "SE"
    ii        <- ii + 1
    mat[, ii] <- c(x$beta.robustvar, x$Lambda0.robustvar) 
    cnms[ii]  <- "SE.robust"
  } else {
    if ("beta.var.estimated" %in% xnms) {
      ii        <- ii + 1
      mat[, ii] <- c(x$beta.var.estimated, x$Lambda0.var.estimated)
      cnms[ii]  <- "SE.est"
      ii        <- ii + 1
      mat[, ii] <- c(x$beta.robustvar.estimated, x$Lambda0.robustvar.estimated) 
      cnms[ii]  <- "SE.robust.est"
    }
    if ("beta.var.design" %in% xnms) {
      ii        <- ii + 1
      mat[, ii] <- c(x$beta.var.design, x$Lambda0.var.design)
      cnms[ii]  <- "SE.design"
      ii        <- ii + 1
      mat[, ii] <- c(x$beta.robustvar.design, x$Lambda0.robustvar.design)
      cnms[ii]  <- "SE.robust.design"
    }
  }
  if (ii < ncol(mat)) {
    mat  <- mat[, 1:ii, drop=FALSE]
    cnms <- cnms[1:ii]
  }
  colnames(mat) <- cnms 
  mat[, 2:ii]   <- sqrt(mat[, 2:ii])
  mat           <- round(mat, digits=dig)
  print(mat)

  invisible(x)

} 
