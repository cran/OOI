#' Outside Option Index
#'
#' calculates the 'outside option index' (defined as
#' \eqn{-\sum P(Z|X) * log(P(Z|X) / P(Z))})
#'  for workers, using employer-employee data.
#'
#' @param formula a formula describing the model to be fitted in order to
#'                estimate P(Z|X) / P(Z). This formula uses a syntax similar to STATA, and so
#'                "x_" refers to all variables with the prefix "x", while
#'                "z_" refers to all variables with the prefix "z". Similarly, "d" refers
#'                to the distance polynomial (see the example below).
#' @param X matrix or data frame with workers characteristics. Note that all column names
#'          should start with "x" (necessary for the inner function 'coef_reshape').
#' @param Z an optional matrix or data frame with jobs characteristics. Note that all column names
#'          should start with "z" (necessary for the inner function 'coef_reshape').
#' @param X.location an optional matrix or data frame with location for workers. Could be
#'                   geographical location (i.e., geo-coordinates) or any other
#'                   feature that can be used in order to measure distance between
#'                   worker and job using 'dist.fun'. Currently the package supports only numeric
#'                   inputs.
#' @param Z.location same as 'X.location' but for jobs.
#' @param wgt an optional numeric vector of weights.
#' @param pred logical. If TRUE (default), predicts the ooi for the provided data.
#' @param method a method for estimating P(Z|X) / P(Z). Currently not in use.
#' @param sim.factor a variable that determines how much fake data to simulate
#'                   (relative to real data).
#' @param dist.fun a distance function to calculate the distance between X.location and
#'                 Z.location. Users interested in using more than one distance metric
#'                 should provide a function that returns for each row of X.location and
#'                 Z.location a vector with all the necessary metrics. Also - the function
#'                 should use columns by their index and not by their names.
#'                 The default function is \code{\link{geo_dist}}, which is suitable
#'                 for data with geo-coordinates.
#' @param dist.order a numeric vector specifying for each distance metric
#'                   an order of the distance polynomial.
#' @param seed the seed of the random number generator.
#'
#' @return An "ooi" object. This object is a list containing
#' the following components:
#'  \item{coeffs}{coefficients from the estimated logit.}
#'  \item{coeffs_sd}{coefficients SE.}
#'  \item{pseudo_r2}{McFadden's pseudo-R squared for the estimated logit.}
#'  \item{standardized_coeffs}{standardized coefficients.}
#'  \item{ooi}{the Outside Option Index.}
#'  \item{hhi}{the Herfindahl-Hirschman Index, an alternative measure for outside options.}
#'  \item{job_worker_prob}{the log probability of each worker to work at his *specific* job (rahter than
#'                         to work at a job with his specific z)}
#'  \item{orig_arg}{a list containing the original arguments (necessary
#'  for \code{\link{predict.ooi}}).}
#'
#' @examples
#' #generate data
#' #worker and job characteristics:
#' n <- 100
#' men <- rbinom(n, 1, 0.5)
#' size <- 1 + rgeom(n, 0.1)
#' size[men == 0] <- size[men == 0] + 2
#' worker_resid <- data.frame(r = round(runif(n, 0, 20), 1))
#' job_location <- data.frame(l = round(runif(n, 20, 40), 1))
#' #prepare data
#' #define distance function:
#' dist_metric <- function(x, y){abs(y - x)}
#' X <- data.frame(men = men)
#' Z <- data.frame(size = size)
#' #add "x" / "z" to column names:
#' X <- add_prefix(X, "x.")
#' Z <- add_prefix(Z, "z.")
#' #estimate P(Z|X) / P(Z) and calculate the ooi:
#' ooi_object <- OOI(formula = ~ x_*z_ + x_*d + z_*d, X = X, Z = Z,
#'                   X.location = worker_resid, Z.location = job_location,
#'                   sim.factor = 3, dist.fun = dist_metric, dist.order = 3)
#' #we can extract the ooi using predict():
#' ooi <- predict(ooi_object)
#' summary(ooi)
#' @export
#' @importFrom stats as.formula binomial contrasts glm model.matrix predict runif terms vcov
#' @importFrom utils setTxtProgressBar txtProgressBar


OOI <- function(formula = NULL,
                X,
                Z = NULL,
                X.location = NULL,
                Z.location = NULL,
                wgt = rep(1, nrow(X)),
                pred = TRUE,
                method = "logit",
                sim.factor = 1,
                dist.fun = geo_dist,
                dist.order = NULL,
                seed = runif(1, 0, .Machine$integer.max)){
  validate_input(X, Z, X.location, Z.location, wgt)
  formula <- validate_formula(formula, Z, X.location)
  if(!is.null(X.location)){
    dist.order <- validate_dist(dist.fun, dist.order, x = X.location[1, ,drop = FALSE],
                                z = Z.location[1, ,drop = FALSE])
  }
  #prepare formula for estimation
  var_names <- c(colnames(X), colnames(Z))
  formula <- prep_form(formula, var_names, dist.order)
  #prepare data for logit estimation
  est_data <- prep_data(X, Z, wgt, sim.factor, seed)
  #merge with distance
  if(!is.null(X.location)){
    if(all(apply(X.location, 2, is.numeric))){
      #rounding is necessary for gen_dist()
      X.location = round(X.location, 3)
    }
    x_loc <- X.location[est_data$worker_id, , drop = FALSE]
    z_loc <- Z.location[est_data$job_id, , drop = FALSE]
    D <- calc_dist(x_loc, z_loc, dist.fun, dist.order)
    est_data <- cbind(est_data, D)
  }
  #estimate logit
  logit <- glm(as.formula(formula), family = binomial(link='logit'),
               data = est_data, weights = est_data$w)
  coeffs <- logit$coefficients
  if(logit$rank < length(coeffs)){
    warning(paste("Logit was estimated on a singular matrix.",
                  "Replacing NA coefficients with 0"))
    coeffs <- replace(coeffs, is.na(coeffs), 0)
  }
  job_worker_prob <- get_probs(logit, est_data$y == 1, wgt)
  coeffs_se <- round(sqrt(diag(vcov(logit))), 4)
  pseudo_r2 <- round(1 - (logit$deviance / logit$null.deviance), 3)
  #calculate standardized coefficients
  est_data <- expand_matrix(est_data[est_data$y == 1,])
  stand_coeffs <- standardize(coeffs[-1], est_data, wgt)
  #prepare output
  orig_arg <- list(X = X,
                    Z = Z,
                    X.location = X.location,
                    Z.location = Z.location,
                    wgt = wgt,
                    dist.fun = dist.fun,
                    dist.order = dist.order)
  output <- list(coeffs = coeffs,
                 coeffs_se = coeffs_se,
                 pseudo_r2 = pseudo_r2,
                 standardized_coeffs = stand_coeffs,
                 orig_arg = orig_arg,
                 job_worker_prob = job_worker_prob,
                 formula = formula)
  #predict OOI
  if(pred){
    #reshape coefficients (necessary for prediction)
    coef_matrices <- coef_reshape(coeffs)
    tmp <- predict_ooi(coef_matrices, X, Z, X.location, Z.location,
                       wgt, dist.fun, dist.order)
    output[["ooi"]] <- tmp$ooi
    output[["hhi"]] <- tmp$hhi
  }
  class(output) <- "ooi"
  return(output)
}


#' Predict Outside Option Index
#'
#' predicts the OOI for new coefficients (for counterfactual analysis)
#' and/or new data.
#'
#' @param object an ooi object.
#' @param new.coef a new *named* vector of coefficients. Check the coefficients produced by
#'                 the main function to see the right format for this vector.
#' @param new.X a new X matrix / data frame.
#' @param new.Z a new Z matrix / data frame.
#' @param new.X.location a new X.location matrix / data frame.
#' @param new.Z.location a new Z.location matrix / data frame.
#' @param new.wgt a new vector of weights
#' @param hhi whether to predict the HHI (Herfindahl-Hirschman Index, an alternative measure for
#'            outside options) instead of the OOI. default is FALSE.
#' @param both whether to return a list with both HHI and OOI when suppling new inputs (default is FALSE).
#'             Necessary especially when predicting takes a lot of time.
#' @param ... further arguments passed to or from other methods.
#'
#' @return If there are no new arguments, returns the original results (ooi/hhi). Otherwise,
#'         returns a vector of ooi/hhi (or a list of both) calculated using the new arguments.
#' @examples
#' #generate data
#' #worker and job characteristics:
#' n <- 100
#' men <- rbinom(n, 1, 0.5)
#' size <- 1 + rgeom(n, 0.1)
#' size[men == 0] <- size[men == 0] + 2
#' worker_resid <- data.frame(r = round(runif(n, 0, 20), 1))
#' job_location <- data.frame(l = round(runif(n, 20, 40), 1))
#' #prepare data
#' #define distance function:
#' dist_metric <- function(x, y){abs(y - x)}
#' X <- data.frame(men = men)
#' Z <- data.frame(size = size)
#' #add "x" / "z" to column names:
#' X <- add_prefix(X, "x.")
#' Z <- add_prefix(Z, "z.")
#' #estimate P(Z|X) / P(Z) and calculate the ooi:
#' ooi_object <- OOI(formula = ~ x_*z_ + x_*d + z_*d, X = X, Z = Z,
#'                   X.location = worker_resid, Z.location = job_location,
#'                   sim.factor = 3, dist.fun = dist_metric, dist.order = 3)
#' #we can extract the ooi using predict():
#' ooi <- predict(ooi_object)
#' #or the hhi:
#' ooi <- predict(ooi_object, hhi = TRUE)
#' #we can also estimate the ooi with different coefficients:
#' coeffs <- ooi_object$coeffs
#' coeffs[names(coeffs) == "x.men"] <- 0
#' new_ooi <- predict(ooi_object, new.coef = coeffs)
#' #or new data:
#' Z2 <- data.frame(z.size = 1 + rgeom(n, 0.1))
#' new_ooi <- predict(ooi_object, new.Z = Z2)
#' @export

predict.ooi <- function(object,
                        new.coef = NULL,
                        new.X = NULL,
                        new.Z = NULL,
                        new.X.location = NULL,
                        new.Z.location = NULL,
                        new.wgt = NULL,
                        hhi = FALSE,
                        both = FALSE, ...){
  #if there are no new inputs, return pre-computed results.
  new_inp <- list(new.coef, new.X, new.Z, new.X.location,
                  new.Z.location, new.wgt)
  if(all(unlist(lapply(new_inp, is.null)))){
    if(hhi) object$hhi else object$ooi
  } else {
    #prepare inputs
    x <- object$orig_arg
    coeffs <- if(is.null(new.coef)) object$coeffs else new.coef
    X <- if(is.null(new.X)) x$X else new.X
    Z <- if(is.null(new.Z)) x$Z else new.Z
    X.location <- if(is.null(new.X.location)) x$X.location else new.X.location
    Z.location <- if(is.null(new.Z.location)) x$Z.location else new.Z.location
    wgt <- if(is.null(new.wgt)) x$wgt else new.wgt
    #reshape coefficients (necessary for prediction)
    coef_matrices <- coef_reshape(coeffs)
    #predict OOI
    validate_input(X, Z, X.location, Z.location, wgt, allow.dif.rows = TRUE)
    tmp <- predict_ooi(coef_matrices, X, Z, X.location, Z.location,
                       wgt, x$dist.fun, x$dist.order)
    if(both){
      tmp
    } else {
      if(hhi) tmp$hhi else tmp$ooi
    }
  }
}

#' Add prefix
#'
#' Adds a prefix to the column names of a matrix / data.frame.
#'
#' @param df a data.frame or a matrix.
#' @param prefix a prefix to be added.
#'
#' @return a matrix / data.frame with new column names.
#' @export
add_prefix <- function(df, prefix){
  col_names <- colnames(df)
  col_names <- paste0(prefix, col_names)
  colnames(df) <- col_names
  return(df)
}
