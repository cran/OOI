cran_skip = TRUE #skip those tests when submitting to cran

library(OOI)
set.seed(123)

#simulate data (for which we can calculate things theoretically)
n <- 2000
men_rate <- 0.1
women_rate <- 0.2
men <- rbinom(n, 1, 0.5)
X_loc <- matrix(runif(n, 0, 100), ncol = 1)
dist <- rep(NA, n)
men_inc <- men == 1
dist[men_inc] <- rexp(n = sum(men_inc), rate = men_rate) #distance for men
dist[!men_inc] <- rexp(n = sum(!men_inc), rate = women_rate) #distance for women
direction <- sample(c(1, -1), size = n, TRUE)
Z_loc <- matrix(rep(NA, n), ncol = 1)
Z_loc[men_inc,] <- X_loc[men_inc,] + dist[men_inc] * direction[men_inc]
Z_loc[!men_inc,] <- X_loc[!men_inc,] + dist[!men_inc] * direction[!men_inc]
X <- matrix(men, ncol = 1, dimnames = list(NULL, "x.men"))

#define simple distance function
dis_function <- function(x, y){abs(x - y)}

#choose workers who are far enough from the edges (for them p(z) = 1/n is reasonable)
q25 <- quantile(X_loc[,1], probs = 0.3)
q75 <- quantile(X_loc[,1], probs = 0.7)
central <- (X_loc[,1] > q25) & (X_loc[,1] < q75)

test_that("OOI returns correct output", {
  skip_if(cran_skip)
  ooi_obj <- suppressWarnings(OOI(~ x_ * d, X = X, X.location = X_loc, Z.location = Z_loc,
                                  dist.fun = dis_function, dist.order = 1, sim.factor = 2))
  ooi <- ooi_obj$ooi
  ooi_men <- ooi[men_inc & central]
  ooi_women <- ooi[!men_inc & central]
  #theoretical results
  theo_ooi_men <- 1 - log(50 * men_rate)
  theo_ooi_women <- 1 - log(50 * women_rate)
  ooi_men_hat <- mean(ooi_men)
  ooi_women_hat <- mean(ooi_women)
  expect_true(abs(ooi_men_hat - theo_ooi_men) < 0.1 &
                abs(ooi_women_hat - theo_ooi_women) < 0.1)
})

test_that("OOI returns the right job_worker probabilities", {
  skip_if(cran_skip)
  ooi_obj <- suppressWarnings(OOI(~ x_ * d, X = X, X.location = X_loc, Z.location = Z_loc,
                                  dist.fun = dis_function, dist.order = 1, sim.factor = 2))
  logp_hat <- ooi_obj$job_worker_prob
  logp <- log(0.5 * men_rate * exp(-dist * men_rate)) #men
  logp[!men_inc] <- log(0.5 * women_rate * exp(-dist[!men_inc] * women_rate)) #women
  p <- exp(logp)
  wgt <- 1/n
  p <- p * wgt
  sum_over_p <- sum(p)
  #normalize p to sum to 1
  logp <- logp - log(sum_over_p)
  diff <- abs(logp_hat - logp)
  expect_true(mean(diff) < 0.1)
})

test_that("OOI returns correct output for multi-dimensional distance", {
  skip_if(cran_skip)
  set.seed(1312)
  #add another distance dimension - city. workers and jobs are always in
  #the same city
  city <- rbinom(n, 1, 0.5)
  X_loc <- cbind(X_loc, city); Z_loc <- cbind(Z_loc, city)
  dis_function <- function(x, y){c(abs(x[1] - y[1]), 1*(x[2] != y[2]))}
  ooi_obj <- suppressWarnings(OOI(~ x_ * d , X = X, X.location = X_loc, Z.location = Z_loc,
                                  dist.fun = dis_function, dist.order = c(1,1), sim.factor = 1))

  ooi <- ooi_obj$ooi
  ooi_men <- ooi[men_inc & central]
  ooi_women <- ooi[!men_inc & central]
  #theoretical results
  theo_ooi_men <- 1 - log(50 * men_rate) - log(2)
  theo_ooi_women <- 1 - log(50 * women_rate) - log(2)
  ooi_men_hat <- mean(ooi_men)
  ooi_women_hat <- mean(ooi_women)
  expect_true(abs(ooi_men_hat - theo_ooi_men) < 0.1 &
                abs(ooi_women_hat - theo_ooi_women) < 0.1)
})

#cuurently we dont use this test:

# #an auxiliary function that samples from the unit circle
# samp_from_unit <- function(n){
#   t <- 2*pi*runif(n)
#   u <- runif(n) + runif(n)
#   r <- ifelse(u > 1, 2-u, u)
#   res <- data.frame(x = r*cos(t), y = r*sin(t))
#   res
# }
#
# #currently this test doesnt pass
# test_that("OOI returns correct output for multi-dimensional distance", {
#   n <- 1000
#   men_rate <- 2
#   women_rate <- 4
#   men <- rbinom(n, 1, 0.5)
#   #sample workers location from the unit circle
#   X_loc <- samp_from_unit(n)
#   dist <- rep(NA, n)
#   men_inc <- men == 1
#   dist[men_inc] <- rexp(n = sum(men_inc), rate = men_rate) #distance for men
#   dist[!men_inc] <- rexp(n = sum(!men_inc), rate = women_rate) #distance for women
#   angle <- 2*pi*runif(n)
#   Z_loc <- X_loc + cbind(dist*cos(angle), dist*sin(angle))
#   dis_function <- function(x, z){c(abs(x[1] - z[1]), abs(x[2] - z[2]))}
#   X <- matrix(men, ncol = 1, dimnames = list(NULL, "x.men"))
#   ooi_obj <- suppressWarnings(OOI(~ x_ * d, X = X, X.location = X_loc, Z.location = Z_loc,
#                                   dist.fun = dis_function, dist.order = c(2,2), sim.factor = 3))
#   ooi <- ooi_obj$ooi
#   #choose workers who are close enough to th center
#   central <- apply(abs(X_loc) < 0.5, 1, all)
#   ooi_men <- ooi[men_inc & central]
#   ooi_women <- ooi[!men_inc & central]
#   #theoretical results
#   theo_ooi_men <- 1 - log(0.5 * men_rate)
#   theo_ooi_women <- 1 - log(0.5 * women_rate)
#   ooi_men_hat <- mean(ooi_men, na.rm = T)
#   ooi_women_hat <- mean(ooi_women, na.rm = T)
#   expect_true(abs(ooi_men_hat - theo_ooi_men) < 0.1 &
#                 abs(ooi_women_hat - theo_ooi_women) < 0.1)
# })


test_that("OOI returns the same results for matrices and data frames with factors", {
  skip_if(cran_skip)
  #simulate data (matrices and data frames)
  n <- 50
  men <-rbinom(n, 1, 0.5)
  native <- rbinom(n, 1, 0.5)
  size <- rbinom(n, 3, 0.5)
  wage <- rnorm(n, 100, 10)
  X_df <- data.frame(x.men = factor(men, c(0,1), c("yes", "no")),
                     x.native = factor(native, c(0,1), c("yes", "no")))
  X_mat <- as.matrix(cbind(men, native))
  X_mat <- add_prefix(X_mat, "x.")
  Z_df <- data.frame(z.size = factor(size, c(0,1,2,3), c("A", "B", "C", "D")),
                     z.wage = wage)
  Z_mat <- cbind(wage, A = 1*(size == 1), B = 1*(size == 0),
                 C = 1*(size == 3))
  Z_mat <- add_prefix(Z_mat, "z.")
  X_loc <- matrix(runif(2*n, 40, 42), ncol = 2)
  Z_loc <- matrix(runif(2*n, 40, 42), ncol = 2)
  #compare results:
  mat_results <- suppressWarnings(OOI(~  x_*z_ + z_*d + x_*d, X = X_mat,
                                      Z = Z_mat, X_loc, Z_loc, sim.factor = 3,
                                      seed = 2))
  df_results <- suppressWarnings(OOI(~  x_*z_ + z_*d + x_*d, X = X_df,
                                     Z = Z_df, X_loc, Z_loc, sim.factor = 3,
                                     seed = 2))
  expect_true(max(abs(mat_results$ooi - df_results$ooi)) < 0.01)
})
