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
  # inla.spde.make.A truncates non-integer group values (2.9999 -> 2), and
  # bam(discrete=TRUE) can hand us discretised times that aren't exact
  # integers, so round to the nearest time index first
  time <- round(x[, 3])

  mesh <- object$xt$mesh
  mesh.time <- object$xt$mesh.time

  Q.time <- mesh.time$qt

  # n.group fixes the number of time blocks to the size of Q.time, so the
  # columns of A line up with the penalty even if not every time is observed
  A <- inla.spde.make.A(mesh, loc = loc.space, group = time,
                        n.group = nrow(Q.time))
  object$X <- as.matrix(A)

  fem <- inla.mesh.fem(mesh)

  # INLA orders the columns of A with the mesh index fastest within each time
  # block (column = node + n.mesh * (time - 1)), so the precision is
  # Q.time (x) Q.space, i.e., time goes first in the Kronecker product
  object$S <- list()
  object$S[[1]] <- as.matrix(kronecker(Q.time, fem$c1))
  object$S[[2]] <- as.matrix(kronecker(Q.time, 2 * fem$g1))
  object$S[[3]] <- as.matrix(kronecker(Q.time, fem$g2))

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

  # n.group must match the fitted model, otherwise predicting at a subset of
  # times gives too few columns
  Xp <- inla.spde.make.A(object$mesh, x[, 1:2, drop=FALSE], group = round(x[, 3]),
                         n.group = nrow(object$mesh.time))
  return(as.matrix(Xp))
}
