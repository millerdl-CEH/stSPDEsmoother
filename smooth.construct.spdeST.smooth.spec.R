# Analysis for "Scalar-on-Function Regression with space-time SPDE Smoothing"
# Erin Bryce and David L Miller
# 2026

# smoother setup functions for mgcv

smooth.construct.spdeST.smooth.spec <- function(object, data, knots) {
  dim <- length(object$term)
  if (dim != 3) stop("Need 3D smooth: x1, x2, t")

  x <- cbind(data[[object$term[1]]],
             data[[object$term[2]]],
             data[[object$term[3]]])
  loc.space <- x[, 1:2]
  time <- x[, 3]

  mesh <- object$xt$mesh
  mesh.time <- object$xt$mesh.time

  A <- inla.spde.make.A(mesh, loc = loc.space, group = time)
  object$X <- as.matrix(A)

  fem <- inla.mesh.fem(mesh)
  Q.time <- mesh.time$qt

  object$S <- list()
  object$S[[1]] <- as.matrix(kronecker(fem$c1, Q.time))
  object$S[[2]] <- as.matrix(kronecker(2 * fem$g1, Q.time))
  object$S[[3]] <- as.matrix(kronecker(fem$g2, Q.time))

  object$L <- matrix(c(2, 2, 2, 4, 2, 0), ncol = 2)
  object$rank <- rep(ncol(object$X), 3)
  object$null.space.dim <- 0
  object$df <- ncol(object$X)

  object$mesh <- mesh
  object$mesh.time <- mesh.time$qt

  class(object) <- "spdeST.smooth"
  return(object)
}

Predict.matrix.spdeST.smooth <- function(object, data){

  x <- matrix(0, nr = length(data[[1]]), nc = 3)
  x[, 1] <- data[[object$term[1]]]
  x[, 2] <- data[[object$term[2]]]
  x[, 3] <- data[[object$term[3]]]

  Xp <- inla.spde.make.A(object$mesh, x[, 1:2, drop=FALSE], group = x[, 3])
  return(as.matrix(Xp))
}
