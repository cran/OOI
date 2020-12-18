validate_input <- function(X, Z, X.location, Z.location, wgt,
                           allow.dif.rows = FALSE){ #do we allow for different rows in X & Z?
  if(is.null(Z) & is.null(X.location) & is.null(Z.location)){
    stop("Z or locations are needed")
  }
  if((!is.null(X.location) & is.null(Z.location)) |
     (is.null(X.location) & !is.null(Z.location))){
    stop("both Z.location and X.location are needed")
  }
  if(any(is.na(X)))
    stop("X contains NA values, which aren't allowed")
  n <- nrow(X)
  validate_type(Z, "Z", n, allow.dif.rows)
  validate_type(Z.location, "Z.location", n, allow.dif.rows)
  validate_type(X.location, "X.location", n, allow.dif.rows)
  validate_colnames(X, "x")
  validate_colnames(Z, "z")
  if((!inherits(wgt, "numeric") | n != length(wgt) | any(is.na(wgt))) & !allow.dif.rows){
    stop(paste("wgt must be numeric,",
               "with the same number of examples as X.",
               "missing values aren't allowed"))
  }
}

validate_formula <- function(formula, Z, X.location){
  #validate formula
  if(is.null(formula)){
    if(is.null(Z)){
      message("formula is missing. default formula (~ x_*d) is used")
      return(~ x_ * d)
    } else {
      if(is.null(X.location)){
        message("formula is missing. default formula (~ x_*z_) is used")
        return(~ x_ * z_)
      } else {
        message("formula is missing. default formula (~ x_*z_ + x_*d) is used")
        return(~ x_ * z_ + x_*d)
      }
    }
  } else {
    terms <- labels(terms(formula))
    fchar <- substr(terms, 1, 1)
    if(is.null(X.location)){
      if(any(fchar == "d"))
        stop("d appears in the formula but locations are missing")
    } else {
      if(!any(fchar == "d"))
        stop("locations were provided but d doesn't appear in the formula")
    }
    return(formula)
  }
}

validate_colnames <- function(df, char){
  if(!is.null(df)){
    col_names <- colnames(df)
    fchar <- substr(col_names, 1, 1)
    if(any(fchar != char) | is.null(col_names)){
      stop(paste("All column names in",
                 ifelse(char == "x", "X", "Z"),
                 "should start with", char))
    }
  }
}

validate_type <- function(df, df.name, exp.length, allow.dif.rows){
  if(!is.null(df)){
    #files read by read_dta sometimes have a strange class
    if(!inherits(df, c("matrix","data.frame")) | length(class(df)) > 1){
      stop(paste(df.name, "should be either matrix or data.frame"))
    }
    if((nrow(df) != exp.length) & !allow.dif.rows){
      stop(paste("X and", df.name, "don't have the same number of rows"))
    }
    if(any(is.na(df))){
      stop(paste(df.name, "contains NA values, which aren't allowed"))
    }
  }
}


validate_dist <- function(d.fun, d.order, x, z){
  dist_output <- d.fun(x, z)
  dim_fun <- length(dist_output)
  if(is.null(d.order)){
    d.order = rep(1, dim_fun)
  } else {
    dim_order <- length(d.order)
    if(dim_fun != dim_order){
      stop(paste0("dist.fun returns ", dim_fun,
                 "-dimensional result, while dist.order is ", dim_order, "-dimensional"))
    }
  }
  return(d.order)
}

