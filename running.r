rm(list=ls())
library(fbttbchol)

k <- 50
n <- 50
x1 <- seq(0, 1, length.out=k)
x2 <- seq(0, 1, length.out=n)
coords <- as.matrix(expand.grid(x1, x2))

covlist <- make_covlist(coords, 10)

CC <- covlist[[1]]

f_test <- f_chol(CC, k, n)
f_test_inv <- f_invchol(CC, k, n)

arma_test <- arma_chol(CC)
arma_test_inv <- arma_invchol(CC)

max(f_test - arma_test)
max(f_test_inv - arma_test_inv)


rbenchmark::benchmark(
  fichol <- f_invchol_list(coords, k, n, covlist),
  aichol <- arma_invchol_list(coords, covlist), replications=1)

