check_string <- function(x, name, valid) {

  len <- length(x)
  if (!is.character(x) || (len != 1)) {
    str <- paste0(valid, collapse=", ")
    stop(paste0("ERROR: ", name, " must be one of ", str))
  }
  if (!(x %in% valid)) {
    str <- paste0(valid, collapse=", ")
    stop(paste0("ERROR: ", name, " must be one of ", str))
  }
  NULL
}

check_time.scale <- function(x, name="time.scale") {

  valid <- c("age", "time-on-study")
  if (!isString(x)) stop("ERROR: ", name, " must be ", valid[1], " or ", valid[2])
  x <- tolower(removeWhiteSpace(x))
  if (!(x %in% valid)) stop("ERROR: ", name, " must be ", valid[1], " or ", valid[2])
  x

}

check_overlap <- function(set1, set2, name1, name2) {

  if (length(intersect(set1, set2))) {
    msg <- paste0("ERROR: ", name1, " and ", name2, " cannot have variables in common")
    stop(msg)
  }
  NULL
}

check_outcome <- function(x, data, name="status") {

  check_var(x, name, colnames(data), exist=1) 
  check_outcome_vec(data[, x, drop=TRUE], name)
  NULL
}

check_outcome_vec <- function(x, name) {

  type <- getVarType(x, only.finite=1)
  if (type != "binary") stop(paste0("ERROR: ", name, " must be a binary variable"))
  NULL
 
}

check_time_vars <- function(data, vars, name="time.vars") {

  check_vars(vars, name, colnames(data), maxlen=2, minlen=1)
  n  <- length(vars)
  x1 <- data[, vars[1], drop=TRUE]
  check_numeric(x1, vars[1], pos=1, len=0) 
  if (n > 1) {
    x2 <- data[, vars[2], drop=TRUE]
    check_numeric(x2, vars[2], pos=1, len=0)
    if (vars[1] == vars[2]) stop(paste("ERROR: duplicate variable name in ", name))
    tmp <- x1 > x2
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      msg <- paste0("ERROR: there are subjects with ", vars[1], " > ", vars[2])
      stop(msg)
    } 
  }
  NULL
}

check_weights <- function(var, data, name="weights") {
 
  cx <- colnames(data)
  check_var(var, name, cx, exist=0)  
  if (length(var)) {
    x <- data[, var, drop=TRUE]
    check_numeric(x, name, pos=1, len=0, maxval=Inf) 
  }
  NULL
}

check_strata <- function(var, data, name="strata") {
 
  cx <- colnames(data)
  check_var(var, name, cx, exist=0)  
  if (length(var)) {
    x <- data[, var, drop=TRUE]
    if (length(unique(x)) > nrow(data)/2) {
      stop("ERROR with stratification variable ", var)
    }
  }
  NULL
}


check_data <- function(x, name="data") {

  if (!is.data.frame(x)) stop(paste0("ERROR: ", name, " must be a data frame"))
  if (!nrow(x)) stop(paste0("ERROR: ", name, " contains no rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", name, " contains no columns"))
  NULL

}

check_risk.data <- function(x, covars, data, name="risk.data") {

  if (is.null(x)) return(NULL)
  if (!length(covars)) return(NULL)
  if (is.list(x)) {
    x <- try(as.data.frame(x, stringsAsFactors=FALSE), silent=TRUE)
    if ("try-error" %in% class(x)) {
      stop(paste0("ERROR: ", name, " could not be coerced into a data frame"))
    }
  }
  if (!is.data.frame(x)) {
    stop(paste0("ERROR: ", name, " must be a data frame or a list that can be coerced to a data frame"))
  }
  if (!nrow(x)) stop(paste0("ERROR: ", name, " contains no rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", name, " contains no columns"))
  tmp  <- !(covars %in% colnames(x))
  if (any(tmp)) {
    miss <- covars[tmp]
    str  <- paste0(miss, collapse=", ")
    msg  <- paste0("ERROR: variables '", str, "' were not found in ", name)
    stop(msg)
  }

  # A column cannot be all NA
  for (v in covars) {
    rvec <- x[, v, drop=TRUE]
    if (all(is.na(rvec))) {
      msg <- paste0("ERROR: column ", v, " in ", name, " contains all missing values")
      stop(msg)
    }

    # Check type of variable against data
    dvec  <- data[, v, drop=TRUE]
    tmp1  <- is.character(rvec) && is.numeric(dvec)
    tmp2  <- is.character(dvec) && is.numeric(rvec)
    if (tmp1 || tmp2) {
      msg <- paste0("ERROR: variable ", v, " has a different type in data and in ", name)
      stop(msg)
    }
  }  

  x
}

check_subset <- function(x, data, name="phase2") {

  check_var(x, name, colnames(data), exist=0) 
  vec <- getDataVec(data, x)
  check_integer(vec, name, len=length(vec), valid=0:1, rem.miss=TRUE)
  NULL

}

check_aux.stat.vars <- function(x, xvars, zvars) {

  nm  <- "aux.vars"
  nmx <- "timeDep.covars"
  nmz <- "timeIndep.covars"

  n <- length(x)
  if (!n) return(NULL)
  if (!is.character(x)) stop(paste0("ERROR: ", nm, " must be NULL or a character vector of column names"))
  valid <- c(xvars, zvars)
  tmp   <- !(x %in% valid)
  if (any(tmp)) {
    msg <- paste0("ERROR: ", nm, " must be variables listed in ", nmx, " or ", nmz)
    stop(msg)
  }
  NULL

}

check_var <- function(x, name, valid, exist=1) {

  if (!exist) {
    if (!length(x)) return(NULL)
  }
  if (!isString(x)) stop(paste0("ERROR: ", name, " must be a column name in the data"))
  if (!(x %in% valid)) stop(paste0("ERROR: ", name, "=", x, " is not a valid column name in the data"))
 
  NULL

}

check_vars <- function(x, name, valid, maxlen=0, minlen=0) {

  len <- length(x)
  if (!len) return(NULL)
  if (maxlen && (len > maxlen)) stop(paste0("ERROR: ", name, " must have length <= ", maxlen))
  if (minlen && (len < minlen)) stop(paste0("ERROR: ", name, " must have length >= ", minlen))

  if (!is.character(x)) stop(paste0("ERROR: ", name, " must be column names in the data"))
  tmp <- !(x %in% valid)
  if (any(tmp)) {
    miss <- x[tmp]
    str  <- paste0(miss, collapse=", ")  
    stop(paste0("ERROR: ", name, " contains the invalid column name(s) ", str))
  }
  NULL

}

check_varnames <- function(varnames, data.cols, data.name="data") {

  if (!length(varnames)) return(NULL)
  tmp <- !(varnames %in% data.cols)
  if (any(tmp)) {
    miss <- paste0(varnames[tmp], collapse=", ")
    msg  <- paste0("ERROR: ", data.name, " does not contain column(s) ", miss)
    stop(msg)
  }
  NULL
}

check_numeric <- function(x, name, pos=1, len=1, maxval=Inf) {

  if (!is.numeric(x)) stop(paste0("ERROR: ", name, " must be numeric"))
  if (len && (len != length(x))) stop(paste0("ERROR: ", name, " must have length ", len))
  if (pos) {
    tmp <- x <= 0
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      if (len == 1) {
        stop(paste0("ERROR: ", name, " must be a positive value")) 
      } else {
        stop(paste0("ERROR: ", name, " must have positive values")) 
      }
    }  
  }
  if (is.finite(maxval)) {
    tmp <- x > maxval
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      stop(paste0("ERROR: ", name, " must be <= ", maxval)) 
    }
  } 

  NULL
}

check_integer <- function(x, name, len=1, minval=-Inf, maxval=Inf, valid=NULL, rem.miss=FALSE) {

  if (len && (len != length(x))) stop(paste0("ERROR: ", name, " must have length ", len))
  if (!is.numeric(x)) stop(paste0("ERROR: ", name, " must be integer"))
  if (rem.miss) {
    tmp <- is.finite(x)
    x   <- x[tmp]
    if (!length(x)) stop(paste0("ERROR: ", name, " contains all non-finite values"))
  }
  tmp <- (x != floor(x))
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) stop(paste0("ERROR: ", name, " must be integer"))  
 
  if (is.finite(minval)) {
    tmp <- x < minval
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      stop(paste0("ERROR: ", name, " must be >= ", minval))  
    }  
  }
  if (is.finite(maxval)) {
    tmp <- x > maxval
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      stop(paste0("ERROR: ", name, " must be <= ", maxval)) 
    }
  } 

  if (length(valid)) {
    tmp <- !(x %in% valid)
    if (any(tmp)) {
      stop(paste0("ERROR: ", name, " contains invalid values")) 
    }
  } 

  NULL
}

check_binary <- function(x, nm, len=1) {

  n <- length(x)
  if (!n) stop(paste0("ERROR: ", nm, " must have values TRUE or FALSE"))
  if (len && (len != n)) stop(paste0("ERROR: ", nm, " must have length ", len))
  if (!all(x %in% 0:1)) stop(paste0("ERROR: ", nm, " must have values TRUE or FALSE"))
  NULL
}

check_status <- function(x, data, name="status") {

  check_var(x, name, colnames(data), exist=1) 
  vec <- getDataVec(data, x)
  check_integer(vec, name, len=length(vec), valid=0:1, rem.miss=TRUE)
  NULL
}

check_covars <- function(x, data, name="covars", len=0) {

  check_vars(x, name, colnames(data))
  if (length(x)) {
    if (any(duplicated(x))) stop(paste0("ERROR: ", name, " contains duplicated names"))
    #for (v in x) {
    #  vec <- getDataVec(data, v)
    #  check_numeric(vec, v, pos=0, len=0, maxval=Inf)
    #}
  } else if (len) {
    stop(paste0("ERROR: at least one covariate must be specified"))
  }
  NULL
}

check_phase2 <- function(x, data, name="phase2") {

  if (!length(x)) return(NULL)
  check_var(x, name, colnames(data), exist=0)  
  if (length(x)) {
    vec <- getDataVec(data, x)
    tmp <- vec %in% 0:1
    if (!all(tmp)) stop(paste0("ERROR: ", name, " must be a binary (0/1) variable"))
  }
  NULL
}

check_phase23 <- function(x2, x3, data, status, name2="subcohort", name3="phase3") {

  check_phase2(x2, data, name=name2)
  if (!is.null(x2) && !is.null(x3)) {
    check_phase2(x3, data, name=name3)
    vec2 <- getDataVec(data, x2)
    y    <- getDataVec(data, status)
    vec2 <- (vec2 %in% 1) | (y %in% 1)

    vec3 <- getDataVec(data, x3)
    tmp  <- (vec3 %in% 1) & !(vec2 %in% 1)
    if (any(tmp)) stop("ERROR: there are phase-3 subjects not in phase-2")   
  }
  NULL
}

check_proxy.vars <- function(x, data, name="proxy.vars") {

  lenx <- length(x)
  if (!lenx) return(NULL)
  if (!is.vector(x) && !is.list(x)) {
    stop(paste0("ERROR: ", name, " must be NULL, a list, or a vector"))
  }
  check_vars(unlist(x), name, colnames(data))
  
  NULL

}

check_predicted.covars <- function(x, data, name="predicted.cox.phase2") {

  lenx <- length(x)
  if (!lenx) return(NULL)
  if (!is.list(x)) {
    stop(paste0("ERROR: ", name, " must be NULL or a named list with names being the phase-two covariates"))
  }
  # Each element of the list must have length 1
  for (i in 1:lenx) {
    if (length(x[[i]]) != 1) {
      msg <- paste0("ERROR: each element of ", name, " must be a single variable")
      stop(msg)  
    }
  }
  check_vars(unlist(x), name, colnames(data))
  
  NULL
}


check_auxvars <- function(x, data, name="aux.vars") {

  if (length(x)) {
    check_vars(x, name, colnames(data))
    if (any(duplicated(x))) stop(paste0("ERROR: ", name, " contains duplicated names"))
    for (v in x) {
      # Check for missing values
      vec <- getDataVec(data, v)
      if (any(is.na(vec))) {
        msg <- paste0("ERROR: ", name, "=", v, " contains missing values")
        stop(msg)
      }
    }
  }
  NULL


}

check_beta_coef <- function(fit) {

  beta <- fit$coefficients
  if (any(is.na(beta))) {
    msg <- paste0("ERROR: not all parameters were estimated, check covariates")
    stop(msg)
  }
  NULL
}

checkOptionListNames <- function(op, valid, name) {
    if (!length(op)) return(NULL)

    # Names cannot be ""
    nms <- trimws(names(op))
    tmp <- nchar(nms) < 1
    if (any(tmp)) {
      stop(paste("ERROR: the list ", name, " is not valid", sep="")) 
    }
    tmp <- !(nms %in% valid)
    if (any(tmp)) {
      err <- paste(nms[tmp], collapse=",", sep="")
      if (length(err) > 1) {
        stop(paste("ERROR: ", err, " are not valid option names for ", name, sep=""))
      } else {
        stop(paste("ERROR: ", err, " is not a valid option name for ", name, sep=""))
      }
    }  

    NULL

} # END: checkOptionListNames

check_calib.op.eta0 <- function(x, name) {

  if (!length(x)) return(NULL)
  check_numeric(x, name, pos=0, len=0, maxval=Inf)
  NULL
}

check_calib.op.eps <- function(x, name) {

  if (!length(x)) return(NULL)
  check_numeric(x, name, pos=1, len=1, maxval=Inf)
  NULL
}

check_calib.op.maxit <- function(x, name) {

  if (!length(x)) return(NULL)
  check_integer(x, name, len=1, minval=1, maxval=Inf, valid=NULL, rem.miss=FALSE)
  NULL
}

check_calibrated.op <- function(x, name="weights.op") {

  if (!length(x)) x <- list()
  if (!is.list(x)) stop(paste0("ERROR: ", name, " must be a list or NULL"))
  valid <- c("niter.max", "epsilon.stop")
  def   <- list(1e4, 1e-10)
  if (length(x)) {
    checkOptionListNames(x, valid, name)
    #check_calib.op.eta0(x[["eta0", exact=TRUE]], "op$eta0")
    check_calib.op.eps(x[["epsilon.stop", exact=TRUE]], "op$epsilon.stop")
    check_calib.op.maxit(x[["niter.max", exact=TRUE]], "op$niter.max")
  }
  x <- default.list(x, valid, def)
  x
}

check_opString <- function(x, name, valid, nchar=1) {

  valid        <- tolower(valid)
  cv           <- substr(valid, 1, 1)
  names(valid) <- cv
  err          <- paste0("'", valid, "'")
  err          <- paste0(err, collapse=", ")
  if (!length(x)) x <- valid[1]
  if (!isString(x)) stop(paste0("ERROR: ", name, " must be one of ", err))
  if (length(x) != 1) stop(paste0("ERROR: ", name, " must be one of ", err))
  x  <- tolower(trimws(x))
  c1 <- substr(x, 1, 1)
  if (!(c1 %in% cv)) stop(paste0("ERROR: ", name, " must be one of ", err))
  x <- valid[c1]
  names(x) <- NULL
  x
}

check_weights.phase3 <- function(x, name="weights.phase3") {

  valid <- c("both", "design", "estimated")
  names(valid) <- c("b", "d", "e")
  if (!length(x)) x <- valid[1]
  if (!isString(x)) stop(paste0("ERROR: ", name, " must be one of 'both', 'design', or 'estimated'"))
  if (length(x) != 1) stop(paste0("ERROR: ", name, " must be one of 'both', 'design', or 'estimated'"))
  x  <- tolower(trimws(x))
  c1 <- substr(x, 1, 1)
  if (!(c1 %in% names(valid))) stop(paste0("ERROR: ", name, " must be one of 'both', 'design', or 'estimated'"))
  x <- valid[c1]
  x
}

check_aux.method <- function(x, name="aux.method") {

  valid <- c("breslow", "shin")
  if (!length(x)) x <- valid[1]
  x <- check_opString(x, name, valid, nchar=1)
  x
}

check_phase2.strata.counts <- function(x, data, stratav, name="subcohort.strata.counts") {

  if (!length(x)) return(NULL)
  if (!is.list(x)) stop(paste0("ERROR: ", name, " must be NULL or a list"))
  nms <- names(x)
  if (!length(nms)) stop(paste0("ERROR: ", name, " must contain the strata names and counts"))
  svec <- unique(data[, stratav, drop=TRUE])
  lenx <- length(x)
  if (lenx != length(svec)) {
    stop(paste0("ERROR: length of ", name, " must equal the number of strata"))
  }
  tmp  <- !(svec %in% nms)
  miss <- nms[tmp]
  if (length(miss)) {
    str <- paste0("'", miss, "'")
    str <- paste0(str, collapse=", ")
    msg <- paste0("ERROR: ", name, " is missing strata ", str)
    stop(msg)
  }
  for (i in 1:length(x)) {
    tmp <- x[[i]]
    if ((length(tmp) != 1) || !is.numeric(tmp) || (tmp < 1)) {
      stop(paste0("ERROR: ", name, " must contain positive counts"))
    }
  }
  NULL
}

check_object <- function(x, name="obj") {

  if (!("casecohortcoxsurv" %in% class(x))) {
    stop(paste0("ERROR: ", name, " must be an object of class 'casecohortcoxsurv'"))
  }
  if (!is.list(x)) stop(paste0("ERROR: ", name, " must be an object of class 'casecohortcoxsurv'"))
  robj <- x[["risk.obj", exact=TRUE]]
  if (is.null(robj)) stop(paste0("ERROR: object 'risk.obj' not found in ", name))
  req <- c("data", "infl", "obj.list")
  for (nm in req) {
    if (is.null(robj[[nm, exact=TRUE]])) {
      stop(paste0("ERROR: risk.obj$", nm, " not found"))
    }
  }

  NULL
}
