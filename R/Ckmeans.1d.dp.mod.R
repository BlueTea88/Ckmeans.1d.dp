#' Ckmeans.1d.dp.mod
#'
#' Modified one-dimensional K-means function.
#'
#' @param x vector of one-dimensional values to cluster
#' @param w list of vectors representing cell weights (one vector per x value)
#' @param z list of vectors used to calculate variance objective (one vector per x value)
#' @param k number of clusters
#'
#' @details
#' Modified function from the Ckmeans.1d.dp R package (attributed to Haizhou Wang and Joe Song).
#'
#' Segments the one-dimensional vector x into clusters without breaking the continuity
#' of x values such that only neighbouring values can be grouped together.
#'
#' Each x value has an associate vector of weights and z values which are used in the
#' calculation of the weighted sum of squares objective.  The clusters are then optimised
#' by minimising the weighted sum of squares (according to the vector of weights and z)
#' such that the continuity of x values is preserved.
#'
#' @return
#' A list of results containing:
#' \itemize{
#'   \item \code{cluster} optimised cluster indices at varying number of clusters
#'   \item \code{centers} cluster centers at varying number of clusters
#'   \item \code{size} the size of each cluster for varying number of clusters change
#' }
#'
#' @useDynLib Ckmeans.1d.dp.mod
#' @export
Ckmeans.1d.dp.mod <- function(x, z, w, k=2)
{
  # Check to see if k is less than 0.
  k <- as.integer(ceiling(k))
  if( k <= 1 ) stop ("Specify k greater than 1.")

  # Check to see if cluster level bigger than the unique number of the input vector
  n.unique <- length(unique(x))

  if(n.unique < 2) {
    stop(paste("Input vector", deparse(substitute(x)), "needs to have at least length 2!\n"))
  } else if(n.unique < k){
    stop("Number of clusters is greater than the unique number of elements in the input vector.")
  }

  # Check that z and w are lists
  if (class(z) != 'list') stop("z input needs to be a list.")
  if (class(w) != 'list') stop("w input needs to be a list.")

  # Check x values
  if (!all(x == sort(x))) stop("x values need to be sorted.")
  if (any(is.na(x))) stop("NA values found for x.")

  # Check y values
  if (any(sapply(z, is.na))) stop("NA values found for z.")
  if (length(z) != length(x)) stop("Expect same length between x and z.")
  if (length(unique(sapply(z, length))) > 1) stop("Expect z vectors to be all the same length.")

  # Check w values
  if (any(sapply(w, is.na))) stop("NA values found for w.")
  if (length(w) != length(x)) stop("Expect same length between x and w.")
  if (length(unique(sapply(w, length))) > 1) stop("Expect w vectors to be all the same length.")
  if (length(w[[1]]) != length(z[[1]])) stop("Expect same dimensions between w and z.")

  # Form data which will be passed to external C++ function
  clusters <- vector("integer", length(x)*k)
  center <- vector("double", 0.5*(1+k)*(k))  # length 1 + 2 + ... + k
  size <- vector("double", 0.5*(1+k)*(k))

  z_long <- c()
  for (i in z) z_long <- c(z_long, i)

  w_long <- c()
  for (i in w) w_long <- c(w_long, i)

  # Call external C++ function
  result <- .C("Ckmeans_1d_dp", PACKAGE='Ckmeans.1d.dp.mod',
               x_data=as.double(x), x_length=as.integer(length(x)),
               z_data=as.double(z_long), z_length=as.integer(length(z[[1]])),
               weight=as.double(w_long), k_in=as.integer(k),
               cluster=as.integer(clusters), centers=as.double(center),
               size=as.double(size))

  # Amend results to a more readable output
  nk1 <- numeric(0)
  nk2 <- numeric(0)
  for (i in 1:k){
    nk1 <- c(nk1, rep(i,length(x)))
    nk2 <- c(nk2, rep(i,i))
  }

  temp <- list()
  temp$cluster <- data.frame(nclusters=nk1, cluster_index=result$cluster)
  temp$centers <- data.frame(nclusters=nk2, centers=result$centers)
  temp$size <- data.frame(nclusters=nk2, size=result$size)
  temp$x <- x
  temp$z <- z
  temp$w <- w

  return(temp)
}
