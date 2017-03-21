# test_Ckmeans.1d.dp.R
#
# Joe Song
# Created: May 3, 2016

library(testthat)
library(Ckmeans.1d.dp.mod)

context("Checking on several examples")

test_that("Weighted input", {

  x <- c(-1, 2, 4, 5, 6)
  y <- c( 4, 3, 1, 1, 1)

  result <- Ckmeans.1d.dp.mod(x, as.list(x), as.list(y), k=3)
  expect_equal(result$size[which(result$size$nclusters==3), 'size'], c(4,3,3))
  expect_equal(result$cluster[which(result$cluster$nclusters==3), 'cluster_index'], c(1,2,3,3,3))
  expect_equal(result$centers[which(result$centers$nclusters==3), 'centers'], c(-1, 2, 5))

  x <- c(1, 2, 3, 4, 5, 6)
  z <- as.list(c(-.9, 1, 1.1, 1.9, 2, 2.1))
  w <- as.list(c( 3,  1,   2,   2, 1, 1))

  result <- Ckmeans.1d.dp.mod(x, z, w, k=3)
  expect_equal(result$size[which(result$size$nclusters==3), 'size'], c(3,3,4))
  expect_equal(result$cluster[which(result$cluster$nclusters==3), 'cluster_index'], c(1,2,2,3,3,3))

})


test_that("Given the number of clusters", {

  x <- c(-1, -1, -1, -1, 2, 2, 2, 4, 5, 6)

  result <- Ckmeans.1d.dp.mod(x, as.list(x), as.list(rep(1,length(x))), 3)
  expect_equal(result$size[which(result$size$nclusters==3), 'size'], c(4,3,3))

  cluster.truth <- c(1,1,1,1,2,2,2,3,3,3)
  expect_equal(result$cluster[which(result$cluster$nclusters==3), 'cluster_index'], cluster.truth)

  centers.truth <- c(-1, 2, 5)
  expect_equal(result$centers[which(result$centers$nclusters==3), 'centers'], centers.truth)
  
})
