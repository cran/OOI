#prepares coeffs for prediction process
coef_reshape <- function(coeffs){
  coef_names <- names(coeffs)
  inter_coef <- coef_names[grepl(":", coef_names)]
  #prepare coef matrix to calculate X*A*Z
  #rows are for x, columns for z
  xz_mat <- prep_matrix(inter_coef, "x", "z", coeffs)
  #prepare coef matrix to calculate X*A*D
  #rows are for x, columns for d
  xd_mat <- prep_matrix(inter_coef, "x", "d", coeffs)
  #prepare coef matrix to calculate D*A*Z
  #rows are for d, columns for z
  dz_mat <- prep_matrix(inter_coef, "d", "z", coeffs)
  #simple terms
  sim_coef <- coef_names[!grepl(":", coef_names)]
  fchar <-  substr(sim_coef, 1, 1)
  z_coef <- coeffs[sim_coef[fchar == "z"]]
  if (length(z_coef) == 0)
    z_coef <- NULL
  d_coef <- coeffs[sim_coef[fchar == "d"]]
  if (length(d_coef) == 0)
    d_coef <- NULL
  return(list(xz_mat = xz_mat, xd_mat = xd_mat, dz_mat = dz_mat,
              z_coef = z_coef, d_coef = d_coef))
}

#prepares matrix of coefficients for interaction terms
prep_matrix <- function(inter_coef, term1, term2, coeffs){
  left_vars <- unlist(lapply(inter_coef, function(x){div_inter(x)$left}))
  right_vars <- unlist(lapply(inter_coef, function(x){div_inter(x)$right}))
  #extract first character
  fchar_left <- substr(left_vars, 1, 1)
  fchar_right <- substr(right_vars, 1, 1)
  #position of relevant interactions
  pos <- (fchar_right == term1 & fchar_left == term2) |
          (fchar_right == term2 & fchar_left == term1)
  if(any(pos)){
    vars1 <- unique(c(left_vars[pos & fchar_left == term1],
                      right_vars[pos & fchar_right == term1]))
    vars2 <- unique(c(left_vars[pos & fchar_left == term2],
                      right_vars[pos & fchar_right == term2]))
    #rows are for term1, columns for term2
    mat <- matrix(, nrow = length(vars1), ncol = length(vars2),
                dimnames = list(vars1, vars2))
    #fill mat with coeffs
    for(r in 1:nrow(mat)){
      for(c in 1:ncol(mat)){
        ind <- (left_vars == vars1[r] & right_vars == vars2[c]) |
          (left_vars == vars2[c] & right_vars == vars1[r])
        mat[r, c] <- coeffs[inter_coef[ind]]
      }
    }
    return(mat)
  } else {
    return(NULL)
  }
}

#this function returns df with two columns: unique - indicate which row in the input is unique
#and bin - bin id (duplicated rows get the same bin id).
#the df is in the same order of original data.
bin_data <- function(dat){
  dat <- as.data.frame(dat) #necessary for order()
  #add id column - so the output will be in the same order of the input
  dat$id <- 1:nrow(dat)
  #sort dat without id (so the cumsum will get us what we want)
  dat_no_id <- dat[ ,!(colnames(dat) %in% "id"), drop = F]
  ord <- do.call(order, dat_no_id)
  dat <- dat[ord, , drop = FALSE]
  dat_no_id <- dat_no_id[ord, , drop = FALSE]
  dat$unique <- !duplicated(dat_no_id)
  dat$bin <- cumsum(dat$unique)
  dat <- dat[order(dat$id), , drop = FALSE]
  return(data.frame(unique = dat$unique, bin = dat$bin))
}

#calculates ooi for the logit method, by looping over workers district
#(district is just grouping of X.location). This is nice because d(i,j) is the
#same for all workers from the same district.

#we are trying to predict log(P(Z|X) / P(Z)), by the following formula:
#log(P(Z|X) / P(Z)) = X1*A1*Z1' + X2*A2*D' + B1*Z2 + B2*D' + A3*"Z3*D"
#and then to calculate the ooi.
predict_ooi <- function(coef.mat, X,
                        Z = NULL,
                        X.location = NULL,
                        Z.location = NULL,
                        wgt = rep(1, nrow(X)),
                        dist.fun = geo_dist,
                        dist.order = 2) {
  #data binning, to speed up prediction
  full_X <- cbind_null(X, X.location)
  tmp <- bin_data(full_X)
  X <- X[tmp$unique, , drop = FALSE]
  if(!is.null(X.location)){
    X.location <- as.matrix(X.location)
    Z.location <- as.matrix(Z.location)
  }
  X.location <- X.location[tmp$unique, , drop = FALSE]
  binX <- tmp$bin[tmp$unique]
  res <- data.frame(id = 1:nrow(tmp), binX = tmp$bin) #to merge the results back later
  full_Z <- cbind_null(Z, Z.location)
  tmp <- bin_data(full_Z)
  Z <- Z[tmp$unique, , drop = F]
  Z.location <- Z.location[tmp$unique, , drop = FALSE]
  binZ <- tmp$bin[tmp$unique]
  wgt_list <- split(wgt, tmp$bin)
  wgt <- sapply(wgt_list, sum) #those are the new weights
  wgt <- wgt[binZ] #keep order
  #################
  n_X <- nrow(X) #currently it has to be exist
  n_Z <- ifelse(is.null(Z), nrow(Z.location), nrow(Z)) #one of the two should exist
  one <- rep(1, n_X)
  X <- expand_matrix(X)
  Z <- expand_matrix(Z)
  A1 <- coef.mat$xz_mat
  if(is.null(A1)){
    #initialize with zeros, so they wouldnt affect calculations.
    X1 <- matrix(rep(0, n_X), ncol = 1)
    A1_Z1 <- matrix(rep(0, n_Z), nrow = 1)
  } else {
    X1 <- X[,rownames(A1)]
    Z1 <- Z[,colnames(A1)]
    A1_Z1 <- A1 %*% t(Z1)
  }
  B1 <- coef.mat$z_coef
  if(!is.null(B1)){
    Z2 <- Z[,names(B1)]
    B1_Z2 <- B1 %*% t(Z2)
    A1_Z1 <- rbind(A1_Z1, B1_Z2)
    rownames(A1_Z1)[nrow(A1_Z1)] <- "cons"
    #add column of 1 (later we calculate 1 %*% B1_Z2 and so adding B1_Z2)
    X1 <- cbind(X1, one)
  }
  #generate district table
  dis_table <- gen_dist(X.location, n_X)
  districts <- unique(dis_table$dis)
  #create a progressbar object
  message("Predicting OOI:", "\n")
  pb <- txtProgressBar(min = 0, max = length(districts),
                       initial = 0, style = 3)
  stepi <- 1
  if(is.null(X.location)){
    for(i in districts){
      workers <- dis_table$worker[dis_table$dis == i]
      logp <- X1[workers,] %*% A1_Z1
      tmp <- calc_ooi(logp, wgt)
      dis_table$ooi[dis_table$dis == i] <- tmp$ooi
      dis_table$hhi[dis_table$dis == i] <- tmp$hhi
      #print and update progress
      setTxtProgressBar(pb,stepi)
      stepi <- stepi + 1
    }
  } else {
    B2 <- coef.mat$d_coef
    A2 <- coef.mat$xd_mat
    A3 <- coef.mat$dz_mat
    #prepare matrices
    are_null <- c(is.null(A2), is.null(A3))
    if(!are_null[1]){ #we need to calculate X2*A2*D' and B2*D'
      X2 <- X[,row.names(A2)]
      #combine B2*D' into A2*D'
      A2 <- rbind(A2, B2)
      X2 <- cbind(X2, one)
      X2_A2 <- X2 %*% A2
      X <- cbind(X1, X2_A2)
    }
    if(!are_null[2]){ #we need to calculate A3*"Z3*D
      Z3 <- Z[,colnames(A3), drop = FALSE]
    }
    #start calculations
    if(all(!are_null)){
      for(i in districts){
        workers <- dis_table$worker[dis_table$dis == i]
        Xi <- X[workers, ,drop = FALSE]
        D <- gen_dist_mat(workers, X.location, Z.location, dist.fun, dist.order)
        A1_Z1i <- A1_Z1
        for(p in 1:sum(dist.order)){
          DpZ3 <- as.vector(D[,p]) * Z3
          A1_Z1i["cons",] <- A1_Z1i["cons", ] + DpZ3 %*% A3[p, ]
        }
        logp <- Xi %*% rbind(A1_Z1i, t(D))
        tmp <- calc_ooi(logp, wgt)
        dis_table$ooi[dis_table$dis == i] <- tmp$ooi
        dis_table$hhi[dis_table$dis == i] <- tmp$hhi
        #print and update progress
        setTxtProgressBar(pb,stepi)
        stepi <- stepi + 1
      }
    } else if (!are_null[1]){
      for(i in districts){
        workers <- dis_table$worker[dis_table$dis == i]
        Xi <- X[workers,]
        D <- gen_dist_mat(workers, X.location, Z.location, dist.fun, dist.order)
        logp <- Xi %*% rbind(A1_Z1, t(D))
        tmp <- calc_ooi(logp, wgt)
        dis_table$ooi[dis_table$dis == i] <- tmp$ooi
        dis_table$hhi[dis_table$dis == i] <- tmp$hhi
        #print and update progress
        setTxtProgressBar(pb,stepi)
        stepi <- stepi + 1
      }
    } else if (!are_null[2]){
      for(i in districts){
        workers <- dis_table$worker[dis_table$dis == i]
        Xi <- X1[workers,]
        D <- gen_dist_mat(workers, X.location, Z.location, dist.fun, dist.order)
        A1_Z1i <- A1_Z1
        for(p in 1:sum(dist.order)){
          DpZ3 <- as.vector(D[,p]) * Z3
          A1_Z1i["cons",] <- A1_Z1i["cons",] + DpZ3 %*% A3[p,]
        }
        A1_Z1i["cons",] <- A1_Z1i["cons",] + D %*% B2
        logp <- Xi %*% A1_Z1i
        tmp <- calc_ooi(logp, wgt)
        dis_table$ooi[dis_table$dis == i] <- tmp$ooi
        dis_table$hhi[dis_table$dis == i] <- tmp$hhi
        #print and update progress
        setTxtProgressBar(pb,stepi)
        stepi <- stepi + 1
      }
    } else {
      for(i in districts){
        workers <- dis_table$worker[dis_table$dis == i]
        Xi <- X1[workers, ,drop = FALSE]
        D <- gen_dist_mat(workers, X.location, Z.location, dist.fun, dist.order)
        A1_Z1i <- A1_Z1
        if(is.null(B1)){
          Xi <- cbind(Xi, rep(1, nrow(Xi)))
          A1_Z1i <- rbind(A1_Z1i, t(D %*% B2))
        } else {
          A1_Z1i["cons",] <- A1_Z1i["cons",] + D %*% B2
        }
        logp <- Xi %*% A1_Z1i
        tmp <- calc_ooi(logp, wgt)
        dis_table$ooi[dis_table$dis == i] <- tmp$ooi
        dis_table$hhi[dis_table$dis == i] <- tmp$hhi
        #print and update progress
        setTxtProgressBar(pb,stepi)
        stepi <- stepi + 1
      }
    }
  }
  #merge back the results
  dis_table$binX <- binX
  res <- merge(x = res, y = dis_table, by = "binX", all.x = T)
  res <- res[order(res$id), , drop = FALSE] #keep the original order
  return(list(ooi = res$ooi, hhi = res$hhi))
}

#calculates ooi from predicted log(P(Z|X) / P(Z)).
#by definition, ooi = P(Z|X) * log(P(Z|X) / P(Z)). so we need to calculate P(Z) and then
#ooi = P(Z) * exp(log(P(Z|X) / P(Z))) * log(P(Z|X) / P(Z)).
#P(Z) are just wgt (normalized to sum to 1).
#the last step is to normalize P(Z|X) to sum to 1.

calc_ooi <- function(logp, wgt){
  wgt <- wgt / sum(wgt)
  one <- rep(1, length(wgt))
  p <- exp(logp)
  p <- t(t(p) * wgt)
  #now p is actually P(Z|X), and we to normalize it to sum to 1 for each worker.
  sum_p <- as.vector(p %*% one)
  p <- p / sum_p
  logp <- logp - log(sum_p)
  #Sum and get the indices
  ooi <- -(p * logp) %*% one
  hhi <- p^2 %*% one
  return(list(ooi = ooi, hhi = hhi))
}

#preforms grouping on X.loc and returns district table. if X.loc is missing, generates
#fake district table.
gen_dist <- function(X.loc = NULL, n){
  if(is.null(X.loc)){
    n_dis <- max(round(n / 50), 2)
    dis <- cut(1:n, n_dis, labels = as.character(1:n_dis))
  } else {
    dis <- apply(X.loc, 1, paste0, collapse = "_")
  }
  dis_table <- cbind.data.frame(worker = 1:n, dis = dis,
                                ooi = rep(NA, n), hhi = rep(NA, n))
  return(dis_table)
}

#generates distance matrix for workers from the same district
gen_dist_mat <- function(workers, X.location, Z.location,
                         dist.fun = geo_dist, dist.order){
  X_loc <- unique(X.location[workers,])
  #replicate X_loc to be compatible with calc_dist
  n <- nrow(Z.location)
  X_loc <- matrix(rep(X_loc, n), nrow = n, byrow = TRUE)
  D <- calc_dist(X_loc, Z.location, dist.fun, dist.order)
  D <- as.matrix(D)
  return(D)
}


