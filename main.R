#' Kernel Function
#'
#' This function computes the fourth order kernel operation.
#'
#' @param dx The input Value.
#' @return The fourth order kernel value.
#' @export
kernal_wrapper <- function(dx) {
  .Call('_SSIndex3_kernal', PACKAGE = 'SSIndex3', dx)
}

#' Shape Parameter Estimator
#'
#' This function conduct . Maximize
#'
#' @param beta Shape Parameter.
#' @param Z Description of Z.
#' @param T2 Description of T2.
#' @param C2 Description of C2.
#' @param Y2 Description of Y2.
#' @param h Description of h.
#' @param ht Description of ht.
#' @return Description of the result.
#' @export
M2_wrapper <- function(beta, Z, T2, C2, Y2, h, ht) {
  .Call('_SSIndex3_M2', PACKAGE = 'SSIndex3', beta, Z, T2, C2, Y2, h, ht)
}

#' Shape Function Estimator
#'
#' This function calculates shapes based on the given parameters (describe its purpose here).
#'
#' @param n Description of n.
#' @param m Description of m.
#' @param midx Description of midx.
#' @param tij Description of tij.
#' @param yi Description of yi.
#' @param xb Description of xb.
#' @param x Description of x.
#' @param t Description of t.
#' @param h Description of h.
#' @param w Description of w.
#' @param med Description of med.
#' @param result Description of result.
#' @return A description of what the function returns.
#' @export
shapeFun3_wrapper <- function(n, m, midx, tij, yi, xb, x, t, h, w, med, result) {
  .Call('_SSIndex3_shapeFun3', PACKAGE = 'SSIndex3', n, m, midx, tij, yi, xb, x, t, h, w, med, result)
}

#' Size Parameter Estimator
#'
#' @param n 
#' @param xr 
#' @param mFhat 
#' @param w 
#'
#' @return
#' @export
#'
#' @examples
size <- function(n, xr, mFhat, w) {
  .Call('_SSIndex3_size', PACKAGE = 'SSIndex3', n, xr, mFhat, w)
}