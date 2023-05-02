
isVector <- function(x) {
  ret <- is.vector(x) && !("list" %in% class(x))
  ret
}

isList <- function(x) {
  ret <- is.list(x) && ("list" %in% class(x))
  ret
}


isString <- function(x) {
  ret <- (length(x) == 1) && is.character(x)
  ret
}

getDataVec <- function(data, col) {

  unlist(data[, col, drop=TRUE])

}

getGlmFormula <- function(yvar, covars) {

  if (length(covars)) {
    str <- paste0(covars, collapse=" + ")
    ret <- paste0(yvar, " ~ ", str) 
  } else {
    ret <- paste0(yvar, " ~ 1")
  }
  ret
}

getFormula <- function(vars) {

  if (length(vars)) {
    str <- paste0(vars, collapse=" + ")
    ret <- paste0(" ~ ", str) 
  } else {
    ret <- " ~ 1"
  }
  ret
}

getModelMatrix <- function(data, vars, rem.intercept=TRUE) {

  form <- as.formula(getFormula(vars))
  mat  <- model.matrix(form, data=data)
  if (rem.intercept) mat <- mat[, -1, drop=FALSE]
  mat

}

getCoxFormula <- function(yvar, timevar, covars) {

  if (length(timevar) == 1) {
    str0 <- paste0("Surv(", timevar, ", ", yvar, ")")
  } else {
    str0 <- paste0("Surv(", timevar[1], ", ", timevar[2], ", ", yvar, ")")
  }
  if (length(covars)) {
    str <- paste0(covars, collapse=" + ")
    ret <- paste0(str0, " ~ ", str) 
  } else {
    ret <- paste0(str0, " ~ 1")
  }
  ret
}

getRandomVariable <- function(prefix, len=4) {

  cc  <- c(letters, as.character(0:9))
  vec <- sample(cc, len, replace=TRUE)
  str <- paste0(vec, collapse="")
  ret <- paste0(prefix, ".", str)
  ret
}

getNonMissing <- function(vec) {

  type <- getVarType(vec, only.finite=0) 
  if (type == "cat") {
    ret <- !is.na(vec)
  } else {
    ret <- is.finite(vec)
  }
  ret
}

removeMiss <- function(data, ncc.subset, vars1, vars2, min.nrow=5, print=1) {

  nr0     <- nrow(data)
  if (is.null(ncc.subset)) ncc.subset <- rep(TRUE, nr0)

  # Remove rows with missing values in vars1
  n1 <- length(vars1)
  if (n1) {
    tmp <- getNonMissing(data[, vars1[1], drop=TRUE])
    if (n1 > 1) {
      for (v in vars1[-1]) {
        tmp <- tmp & getNonMissing(data[, v, drop=TRUE])
      }
    }
    data       <- data[tmp, , drop=FALSE]
    ncc.subset <- ncc.subset[tmp]
  }

  # Remove missing from phase 2 covariates
  if (length(vars2)) {
    tmp <- rep(TRUE, nrow(data))
    for (v in vars2) {
      tmp <- ncc.subset & !getNonMissing(data[, v, drop=TRUE])
      if (any(tmp)) {
        data       <- data[!tmp, , drop=FALSE]
        ncc.subset <- ncc.subset[!tmp]
      }
    } 
  }
  nr <- nrow(data)
  m  <- nr0 - nr
  if (nr < min.nrow) stop("ERROR: after removing missing values, the data contains too few rows")
  if (print && m) {
    if (m > 1) {
      message(paste0("NOTE: ", nr0-nr, " subjects have been removed due to missing values"))
    } else {
      message(paste0("NOTE: 1 subject has been removed due to missing values"))
    }
  }

  list(data=data, ncc.subset=ncc.subset)

}

removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

}

default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} 

getQuotedVarStr <- function(vars, qchar="'", sep=", ") {

  if (!length(vars)) return("")
  vars <- paste0(qchar, vars, qchar)
  ret  <- paste0(vars, collapse=sep)
  ret
}

keepValsInVecList <- function(x, keep) {

  # x is a vector or named list of vectors

  len <- length(x)
  if (!len) return(NULL)

  if (isList(x)) {
    nms <- names(x)
    if (is.null(nms)) stop("INTERNAL CODING ERROR 1")
    for (nm in nms) {
      vals <- x[[nm, exact=TRUE]]
      if (length(vals)) {
        tmp  <- vals %in% keep
        vals <- vals[tmp]
      }
      if (!length(vals)) vals <- NULL
      x[[nm]] <- vals
    }
  } else if (isVector(x)) {
    tmp <- x %in% keep
    x   <- x[tmp]
    if (!length(x)) x <- NULL
  } else {
    stop("INTERNAL CODING ERROR 2")
  }
  x
}

unfactor <- function(fac, fun=NULL) {

  # fac   Factor
  # fun   Function like as.character or as.numeric, etc

  if (is.factor(fac)) {
    ret <- levels(fac)[fac]
  } else {
    ret <- fac
  }

  if (!is.null(fun)) ret <- fun(ret)

  ret

} # END: unfactor
