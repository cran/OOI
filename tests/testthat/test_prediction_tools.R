library(OOI)

set.seed(123)

test_that("coef_reshape returns correct coefficient matrices", {
  coeffs <- round(rexp(12), 2)
  names(coeffs) <- c("x1", "x2", "z1", "z2", "d", "x1:z1",
                     "x1:z2", "x2:z1", "x2:z2", "d:x1", "d:x2", "d:z1")
  b <- coeffs[3:4]
  A <- matrix(coeffs[6:9], ncol = 2, byrow = TRUE)
  A_dist <- coeffs[10:11]
  dz_coef <- coeffs[12]
  output <- coef_reshape(coeffs)
  expect_true(all(output$xz_mat == A) & all(output$z_coef == b) &
                all(output$xd_mat == A_dist) & all(output$dz_mat == dz_coef))
})

#generate data with duplicates (to check that the binning procedure working fine)
n <- 100
m <- 4
w <- rexp(n, 0.5)
dist_order = 2
Xnames <- paste0("x.", 1:m)
Znames <- paste0("z.", 1:m)
X <- matrix(rgeom(n*m, 0.5), ncol = m, dimnames = list(NULL, Xnames))
Z <- matrix(rgeom(n*m, 0.5), ncol = m, dimnames = list(NULL, Znames))
X_loc <- matrix(round(runif(2*n, 40, 42), 0), ncol = 2)
Z_loc <- matrix(round(runif(2*n, 40, 42), 0), ncol = 2)
est_data <- prep_data(X, Z, wgt =  w, sim.factor = 5, seed = 4)
#merge with distance
x_loc <- X_loc[est_data$worker_id,]
z_loc <- Z_loc[est_data$job_id,]
D <- calc_dist(x_loc, z_loc, dist.order = dist_order)
est_data <- cbind(est_data, D)

test_that("predict_ooi returns the same results as predict; formula = ~ x_ * z_ + x_ * d", {
  #prepare formula for estimation
  formula <- prep_form(~ x_ * z_ + x_ * d, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
               data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = TRUE)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})

test_that("predict_ooi returns the same results as predict; formula = ~ x_*z_ + x_*d + d*z_", {
  #prepare formula for estimation
  formula <- prep_form(~ x_*z_ + x_*d + d*z_, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
                                data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = TRUE)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})

test_that("predict_ooi returns the same results as predict; formula = ~ x_*z_ + d*z_", {
  #prepare formula for estimation
  formula <- prep_form(~ x_*z_ + d*z_, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
                                data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = TRUE)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})

test_that("predict_ooi returns the same results as predict; formula = ~ x_ * d", {
  #prepare formula for estimation
  formula <- prep_form(~ x_ * d, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
                                data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = TRUE)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})

test_that("predict_ooi returns the same results as predict; formula = ~ x_ + z_ + d_", {
  #prepare formula for estimation
  formula <- prep_form(~ x_ + z_ + d, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
                                data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = TRUE)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})

test_that("predict_ooi returns the same results as predict; formula = ~ z_ + d", {
  #prepare formula for estimation
  formula <- prep_form(~ z_ + d, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
                                data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = T)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})

test_that("predict_ooi returns the same results as predict; formula = ~ x_ + d", {
  #prepare formula for estimation
  formula <- prep_form(~ x_ + d, c(Xnames, Znames), dist.order = dist_order)
  #estimate logit
  logit <- suppressWarnings(glm(as.formula(formula), family = binomial(link='logit'),
                                data = est_data, weights = est_data$w))
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- predict_ooi(coef_matrices, X, Z, X_loc, Z_loc, wgt = w)$ooi
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = m, byrow = TRUE,
                 dimnames = list(NULL, Xnames))
    X_loc_i <- matrix(rep(X_loc[i,], n), nrow = n, byrow = TRUE)
    D <- calc_dist(X_loc_i, Z_loc, dist.order = dist_order)
    df <- cbind.data.frame(Xi, Z, D)
    f <- t(as.matrix(predict(logit, newdata = df)))
    w <- w / sum(w)
    ooi_predict <- c(ooi_predict, calc_ooi(f, w)$ooi)
  }
  expect_true(max(abs(ooi_predict - ooi_package)) < 0.001)
})
