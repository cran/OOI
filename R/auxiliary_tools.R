#' Geographical distance
#'
#' calculates geo distance between *two* points.
#'
#' @param x.loc a 2-length vector. The first value is for longitude,
#'              the second for latitude.
#' @param z.loc a 2-length vector. The first value is for longitude,
#'              the second for latitude.
#' @return distance in miles.
#' @export
geo_dist <- function(x.loc, z.loc){
  x.loc <- as.numeric(x.loc)
  z.loc <- as.numeric(z.loc)
  x.long <- x.loc[1]
  x.lat <- x.loc[2]
  z.long <- z.loc[1]
  z.lat <- z.loc[2]
  tmp <- ((z.lat-x.lat)/2 * pi/180)^2
  tmp <- tmp + cos(z.lat * pi/180) * cos(x.lat * pi/180) * sin((z.long-x.long)/2 * pi/180)^2
  tmp <- 2 * atan2(tmp^0.5, (1-tmp)^0.5)
  return(tmp * 3959) #distance in miles
}

#calculates 1:1 distance
calc_dist <- function(X.location, Z.location, fun = geo_dist, dist.order){
  nr_X <- nrow(X.location)
  nr_Z <- nrow(Z.location)
  if(nr_X != nr_Z){
    stop("X.location and Z.location should have same number of rows")
  } else {
    dist_dim <- length(dist.order)
    distance <- matrix(ncol = dist_dim, nrow = nr_X)
    for(i in 1:nr_X){
      distance[i,] <- fun(X.location[i, ], Z.location[i, ])
    }
  }
  colnames(distance) <- paste0("d", 1:dist_dim, "1")
  distance <- as.data.frame(distance)
  #add high order distance
  for(i in 1:dist_dim){
    if(dist.order[i] > 1){
      for(j in 2:dist.order[i]){
        distance$tmp <- distance[,paste0("d", i, "1")]^j
        colnames(distance)[names(distance) == "tmp"] <- paste0("d", i, j)
      }
    }
  }
  #Arrange columns in alphabetical order (e.g., d11 d12 d21)
  distance = distance[, colnames(distance)[order(colnames (distance))] , drop = FALSE]
  return(distance)
}

#adds simulated matches to real data
prep_data <- function(X, Z = NULL, wgt = rep(1, nrow(X)),
                      sim.factor = 1, seed = NULL){
  set.seed(seed)
  n <- nrow(X)
  n.fake = round(n * sim.factor)
  real_data <- data.frame(worker_id = 1:n, job_id = 1:n, y = rep(1, n))
  #simulate matches
  worker_fake_id <- sample(n, size = n.fake, replace = TRUE, prob = wgt)
  job_fake_id <- sample(n, size = n.fake, replace = TRUE, prob = wgt)
  fake_data <- data.frame(worker_id = worker_fake_id,
                          job_id = job_fake_id, y = rep(0, n.fake))
  res <- rbind(real_data, fake_data)
  res$y <- as.factor(res$y)
  #merge with X Z & weights
  res$w[res$y == 1] <- wgt
  res$w[res$y == 0] <- mean(wgt) #weights for fake matches
  if(is.null(Z)){
    res <- cbind(res, X[res$worker_id, , drop = FALSE])
  } else {
    res <- cbind(res, X[res$worker_id, , drop = FALSE], Z[res$job_id, , drop = FALSE])
  }
  return(res)
}

#prepares formula for glm
#this function also validates that without distance, X*Z must be included
prep_form <- function(formula, var.names, dist.order){
  terms_labels <- labels(terms(formula))
  n <- length(terms_labels)
  #flag for distance
  d_included <- FALSE
  #flag for X*Z
  xz_included <- FALSE
  #initialize formula
  form <- "~ 1"
  for(i in 1:n){
    term <- terms_labels[i]
    #Is this an interaction term?
    inter_term <- grepl(":", term)
    if(!inter_term){
      if(term == "d"){d_included <- TRUE}
      form <- paste(form, ext_names(term, var.names, dist.order), sep = " + ")
    } else {
      #divide to left vars & right vars
      left_vars <- ext_names(div_inter(term)$left, var.names, dist.order)
      right_vars <- ext_names(div_inter(term)$right, var.names, dist.order)
      fchar_l <- substr(left_vars, 1, 1)
      fchar_r <- substr(right_vars, 1, 1)
      if(fchar_r == "x" & fchar_l == "z" | fchar_r == "z" & fchar_l == "x")
        xz_included <- TRUE
      #paste back
      tmp <- paste0("(", left_vars, ")", ":", "(", right_vars, ")")
      form <- paste(form, tmp, sep = " + ")
    }
  }
  if(!d_included & !xz_included)
    stop("either distance ('d') or interaction between X & Z must be included")
  form <- paste("y", form, sep = " ")
  return(form)
}

#divides interaction term into left term and right term
div_inter <- function(inter_term){
  point_pos <- gregexpr(pattern = ":", inter_term)[[1]][1]
  res <- list(left = substr(inter_term, 1, point_pos - 1),
              right = substr(inter_term, point_pos + 1, nchar(inter_term)))
  return(res)
}


#extracts original variable names from term
ext_names <- function(term, var.names, dist.order){
  #this term ends with "_"?
  n.char <- nchar(term)
  last_char <- substr(term, n.char, n.char)
  if(last_char == "_"){
    #take all variables starting with the expression to the left
    term_init <-  substr(term, 1, n.char - 1)
    names_init <- substr(var.names, 1, n.char - 1)
    res <- var.names[term_init == names_init]
    #Is this a distance term?
  } else if (term == "d"){
    #generate d#i#j where i is the i-th distance metric and j is the power
    dist_dim <- length(dist.order)
    dist_terms <- paste0("d", 1:dist_dim)
    res <- apply(cbind(dist_terms, dist.order), 1,
                        function(x){paste0(x[1], seq(1, x[2], 1))})
    res <- unlist(res)
  } else {
    res <- term
  }
  return(paste(res, collapse = " + "))
}


#standardizes coefficients
standardize <- function(coeffs, dat, wgt){
  #calculate sd for relevant variables
  coef_names <- names(coeffs)
  inter_pos <- grepl(":", coef_names)
  rel_vars <- coef_names[!inter_pos] #variables without interaction
  sd <- apply(dat[, rel_vars], 2,
              function(x, w = wgt){sqrt(modi::weighted.var(x, w))})
  coeffs[rel_vars] <- coeffs[rel_vars] * sd
  #for interaction terms, we need to mulitply by the SD of each variable
  inter_pos <- which(inter_pos)
  for(i in inter_pos){
    term <- names(coeffs[i])
    var1 <- div_inter(term)$left
    var2 <- div_inter(term)$right
    coeffs[i] <- coeffs[i] * sd[var1] * sd[var2]
  }
  return(coeffs)
}

#converts data.frame to matrix and expands factors to a set of dummy variables
#(including reference category)
expand_matrix <- function(df){
  if(is.matrix(df) | is.null(df)){
    return(df)
  }
  factors_ind <- sapply(df, is.factor)
  if(sum(factors_ind) == 0){
    return(as.matrix(df))
  }
  contrasts_arg <- lapply(data.frame(df[,factors_ind]),
                          contrasts, contrasts = FALSE)
  names(contrasts_arg) <- colnames(df)[factors_ind]
  df <- model.matrix( ~ .-1, data = df, contrasts.arg = contrasts_arg)
  return(df)
}


#cbind that can handle empty data frames
cbind_null <- function(df1, df2){
  if(is.null(df1)){
    return(df2)
  } else if(is.null(df2)){
    return(df1)
  } else {
    return(cbind(df1, df2))
  }
}

#extract the log probability of each worker to work at his specific job from the estimated logit model
#logit - a logit model (object)
#indices - indices of the original data (and not the simulated). logical
get_probs <- function(logit, indices, wgt){
  wgt <- wgt / sum(wgt)
  logp <- predict(logit)[indices]
  p <- exp(logp)
  p <- p * wgt
  sum_over_p <- sum(p)
  #normalize p to averaged to 1
  logp <- logp - log(sum_over_p)
  return(logp)
}




