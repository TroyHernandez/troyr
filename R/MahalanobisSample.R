
#' @export
MahalanobisSample <- function(x, n){
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  V.sx <- svd(xs, nu = 0)$v
  rng.x1 <- apply(xs %*% V.sx, 2, range)
  z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1],
                                               max = M[2]), nn = n)
  z <- tcrossprod(z1, V.sx) + m.x
  z
}
