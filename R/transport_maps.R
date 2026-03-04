library(distributional)
library(tibble)
library(dplyr)

# ==============================================================================
# Dispatch layer: eval_cdf / eval_quantile
#
# These two S3 generics are the sole extension point for adding new distribution
# types (e.g., simulation-based).  Any object that implements these two methods
# works transparently with everything downstream.
# ==============================================================================

#' Evaluate the CDF of a distribution at one or more points
#' @param d A distribution object (distributional or custom shift_dist)
#' @param x Numeric vector of evaluation points
#' @return Numeric vector of probabilities in (0, 1)
eval_cdf <- function(d, x) UseMethod("eval_cdf")

#' Evaluate the quantile (inverse-CDF) of a distribution at one or more probs
#' @param d A distribution object
#' @param p Numeric vector of probabilities in (0, 1)
#' @return Numeric vector of quantiles
eval_quantile <- function(d, p) UseMethod("eval_quantile")

# Methods for <distribution> objects from the distributional package
# Note: cdf() is exported by distributional; quantile() dispatches via stats::quantile
eval_cdf.distribution <- function(d, x) unlist(distributional::cdf(d, x))
eval_quantile.distribution <- function(d, p) unlist(stats::quantile(d, p))

# ==============================================================================
# Custom distribution constructor (simulation-ready path)
#
# Wrap any pair of (cdf_fn, quantile_fn) closures into a shift_dist object that
# satisfies both generics.  Enables simulation-derived or non-standard
# distributions to work with all shift-representability tests without changes.
# ==============================================================================

#' Create a custom distribution from cdf / quantile functions
#' @param cdf_fn      Function(x) -> probabilities in (0,1)
#' @param quantile_fn Function(p) -> quantiles
#' @param label       Optional character label shown in outputs
#' @return A `shift_dist` object
make_dist <- function(cdf_fn, quantile_fn, label = NULL) {
  structure(
    list(cdf = cdf_fn, quantile = quantile_fn, label = label),
    class = "shift_dist"
  )
}

eval_cdf.shift_dist <- function(d, x) d$cdf(x)
eval_quantile.shift_dist <- function(d, p) d$quantile(p)

print.shift_dist <- function(x, ...) {
  label <- if (!is.null(x$label)) x$label else "<custom>"
  cat("<shift_dist:", label, ">\n")
  invisible(x)
}

# ==============================================================================
# distribution_family(): the primary user-facing entry point
#
# Accepts any mix of distributional / shift_dist objects, validates, and names
# them automatically.  The reference distribution (default: first one) sets the
# common scale for all transport maps.
# ==============================================================================

#' Assemble a family of distributions for shift-representability testing
#'
#' @param ... Distribution objects (distributional or shift_dist), or a single
#'   named list of them.
#' @param reference Integer index of the reference distribution X_1 (default 1).
#' @return A `shift_family` object
#'
#' @examples
#' fam <- distribution_family(
#'   noise  = dist_normal(0, 1),
#'   signal = dist_normal(1, 1)
#' )
distribution_family <- function(..., reference = 1L) {
  dists <- list(...)

  if (length(dists) < 2) stop("A shift family requires at least 2 distributions.")
  names(dists) <- names(dists) %||% paste0("X", seq_along(dists))
  structure(dists, reference = as.integer(reference), class = "shift_family")
}

print.shift_family <- function(x, ...) {
  nms <- names(x)
  ref_label <- nms[attr(x, "reference")]
  cat(sprintf("<shift_family> with %d distributions (reference: %s)\n", length(nms), ref_label))
  cat(paste(nms, collapse = ", "))
  invisible(x)
}

# ==============================================================================
# Transport maps
#
# psi_n(x) = F_1^{-1}( F_n(x) )
# All maps are returned as plain R closures, making composition trivial:
#   psi_n(psi_m(x_grid))  is just  maps[["n"]](maps[["m"]](x_grid))
# ==============================================================================

#' Build a single transport map from d_n to d_ref
#'
#' Returns the closure  psi(x) = F_ref^{-1}(F_n(x)).
#' Probabilities are clipped to (1e-12, 1-1e-12) for numerical stability.
make_transport_map <- function(d_n, d_ref) {
  force(d_n)
  force(d_ref)
  function(x) {
    p <- eval_cdf(d_n, x)
    p <- pmax(1e-12, pmin(1 - 1e-12, p))
    eval_quantile(d_ref, p)
  }
}

#' Return all reference transport maps psi_n for a shift_family
#'
#' The reference map (identity) is excluded from the list.
#'
#' @param dists A `shift_family` object
#' @return Named list of closures; element `"Xn"` is  x -> F_1^{-1}(F_n(x))
transport_maps <- function(dists) {
  ref <- attr(dists, "reference")
  d_ref <- dists[[ref]]
  maps <- lapply(dists[-ref], function(d) make_transport_map(d, d_ref))
  setNames(maps, names(dists)[-ref])
}

#' Evaluate all transport maps over an x_grid, returning a tidy tibble
#'
#' The reference distribution is included as the identity map for plotting.
#'
#' @param dists  A `shift_family` object
#' @param x_grid  Numeric vector of evaluation points
#' @return Tibble with columns: x, map_id (factor), psi_x
eval_transport_maps <- function(dists, x_grid) {
  maps <- transport_maps(dists)
  lapply(maps, do.call, args = list(x_grid))
}
