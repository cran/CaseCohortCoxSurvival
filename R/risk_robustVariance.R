estimatePureRisk <- function(obj, x) {

  check_object(obj)
  new.data <- check_risk.data(x, obj$risk.obj$obj.list$covars.orig, obj$risk.obj$data) 
  ret      <- estimatePureRisk_main(new.data, obj) 
  ret
}

estimatePureRisk_main <- function(risk.data, obj) {

  var.data <- getVarianceForData(risk.data, obj$risk.obj$data, obj$risk.obj$obj.list, 
                                 obj$risk.obj$infl)
  ret <- list()
  tmp <- var.data[["rv", exact=TRUE]]
  if (!is.null(tmp)) ret[["var"]] <- tmp
  tmp <- var.data[["rv.T", exact=TRUE]]
  if (!is.null(tmp)) ret[["var.estimated"]] <- tmp
  tmp <- var.data[["rv.F", exact=TRUE]]
  if (!is.null(tmp)) ret[["var.design"]] <- tmp
  ret
}

getRobustVarForEst <- function(infl) {

  lstnms <- c("infl", "infl.T", "infl.F") 
  prmnms <- c("infl.", "infl2.", "infl3.")
  prmlen <- nchar(prmnms)
  m      <- length(prmnms)

  ret <- list()
  nms <- names(infl)
  for (nm in nms) {
    res <- list()
    lst <- infl[[nm]]
    all <- names(lst)
    for (i in 1:m) {
      tmp  <- substr(all, 1, prmlen[i]) == prmnms[i]
      all2 <- all[tmp]
      if (!length(all2)) next
      for (v in all2) res[[v]] <- robustvariance(lst[[v]])
    }
    ret[[nm]] <- res 
  }
  ret
}

setUpRiskData_catcovars <- function(data, covars, obj) {

  catvars <- obj[["catvars", exact=TRUE]]
  if (!length(catvars)) return(data)
  tmp <- covars %in% catvars
  if (!any(tmp)) return(data)
  vars      <- covars[tmp]
  varlevels <- obj[["varlevels", exact=TRUE]]
  if (!length(varlevels)) stop("INTERNAL CODING ERROR 1")
  ii <- 0
  for (v in vars) {
    vec  <- data[, v, drop=TRUE]
    tmp  <- !is.na(vec)
    vec  <- vec[tmp]
    uvec <- unique(vec)
    if (length(uvec) > 1) next
    levs <- varlevels[[v, exact=TRUE]]
    if (!length(levs)) stop("INTERNAL CODING ERROR 2") 
    tmp  <- !(levs %in% uvec)
    levs <- levs[tmp]
    if (!length(levs)) stop("INTERNAL CODING ERROR 3")
    # Add a dummy row to the data
    ii            <- ii + 1
    row           <- data[1, , drop=FALSE]
    row[1, v]     <- levs[1]
    rownames(row) <- paste0("...row..", ii)
    data          <- rbind(data, row)
  }
  for (v in vars) data[, v] <- factor(data[, v, drop=TRUE])
  data
}

setUpRiskData <- function(data, covars, beta, obj) {

  rn   <- rownames(data)
  if (is.null(rn)) rn <- as.character(1:nrow(data))
  beta.nms <- names(beta)
  nr0      <- nrow(data)

  # Watch out for categorical covariates which are constant in data because
  #   model.matrix will throw an error
  data <- setUpRiskData_catcovars(data, covars, obj)
  rn   <- rownames(data)

  # Initialize return matrix to 0
  ret <- matrix(data=0, nrow=nrow(data), ncol=length(beta))
  rownames(ret) <- rn
  colnames(ret) <- beta.nms

  # Be careful with categorical covariates
  mat  <- getModelMatrix(data, covars, rem.intercept=TRUE)
  cx   <- colnames(mat)
  tmp  <- cx %in% beta.nms
  miss <- cx[!tmp]
  cx   <- cx[tmp]
  if (length(miss)) {
   str <- paste0("'", miss, "'")
   str <- paste0(str, collapse=", ")
   msg <- paste0("The risk data model matrix has columns ", str, 
                 " that are not covariates. They will be dropped")
   warning(msg)
  }
  if (!length(cx)) {
    warning("No common variables between covars and risk data model matrix")
    return(NULL)
  }
  mat <- mat[, cx, drop=FALSE]
  ret[rownames(mat), cx] <- mat

  # Remove rows if dummy rows were added
  if (nrow(ret) > nr0) ret <- ret[1:nr0, , drop=FALSE]
  ret
}

myPasteCols <- function(x, sep="*") {

  nc  <- ncol(x)
  ret <- x[, 1, drop=TRUE]
  if (nc > 1) {
    for (i in 2:nc) ret <- paste(ret, x[, i, drop=TRUE], sep=sep)
  }
  ret
 
}

