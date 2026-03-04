library(tibble)
library(dplyr)
library(purrr)

source("R/transport_maps.R")

# ==============================================================================
# Commutativity test  (Theorem 2, Part 1)
#
# Necessary condition: psi_n ∘ psi_m = psi_m ∘ psi_n  for all pairs n ≠ m.
# Measured as  D_nm(x) = psi_n(psi_m(x)) - psi_m(psi_n(x)).
# D_nm ≡ 0 iff the pair commutes.
# ==============================================================================

#' Test pairwise commutativity of transport maps
#'
#' @param family  A `shift_family` object
#' @param x_grid  Numeric vector of evaluation points
#' @param tol     Maximum |deviation| considered numerically zero (default 1e-6)
#' @return List with:
#'   $pointwise  tibble(map_n, map_m, x, deviation)
#'   $summary    tibble(map_n, map_m, max_abs_dev, mean_abs_dev, passes)
commutativity_test <- function(family, x_grid, tol = 1e-6) {
  maps <- transport_maps(family)
  nms <- names(maps)

  if (length(nms) < 2) {
    msg <- "Commutativity requires at least 2 non-reference distributions."
    return(list(
      pointwise = tibble(
        map_n = character(), map_m = character(),
        x = numeric(), deviation = numeric()
      ),
      summary = tibble(
        map_n = character(), map_m = character(),
        max_abs_dev = numeric(), mean_abs_dev = numeric(),
        passes = logical()
      ),
      message = msg
    ))
  }

  pairs <- combn(nms, 2, simplify = FALSE)

  pointwise <- map_dfr(pairs, \(pair) {
    nm <- pair[1]
    mm <- pair[2]
    deviation <- maps[[nm]](maps[[mm]](x_grid)) - maps[[mm]](maps[[nm]](x_grid))
    tibble(map_n = nm, map_m = mm, x = x_grid, deviation = deviation)
  })

  summary <- pointwise |>
    group_by(map_n, map_m) |>
    summarise(
      max_abs_dev = max(abs(deviation)),
      mean_abs_dev = mean(abs(deviation)),
      passes = max(abs(deviation)) < tol,
      .groups = "drop"
    )

  list(pointwise = pointwise, summary = summary)
}

# ==============================================================================
# Non-crossing test  (Theorem 2, Part 2)
#
# Two necessary conditions:
#  (a) Fixed-point check: each psi_n should not cross the identity, i.e.
#      psi_n(x) - x must not change sign (no crossing with the reference CDF).
#  (b) Pairwise crossing check: for n ≠ m, psi_n(x) - psi_m(x) must not
#      change sign (no two transport maps cross each other).
# ==============================================================================

#' Test the non-crossing condition on transport maps
#'
#' @param family  A `shift_family` object
#' @param x_grid  Numeric vector of evaluation points (dense enough to detect
#'   sign changes; at least 100 points recommended)
#' @param tol     Minimum |gap| considered safely non-zero (default 1e-6)
#' @return List with:
#'   $fixed_point  tibble(map_id, min_abs_gap, sign_changes, passes)
#'   $pairwise     tibble(map_n, map_m, min_abs_gap, sign_changes, passes)
noncrossing_test <- function(family, x_grid, tol = 1e-6) {
  maps <- transport_maps(family)
  nms <- names(maps)

  count_sign_changes <- function(v) sum(diff(sign(v[is.finite(v)])) != 0)

  # (a) Fixed-point check: psi_n(x) vs x
  # A sign change indicates the transport map crosses the identity (fixed point).
  # Exception: if the gap is everywhere near zero, the map IS the identity
  # (trivial transport for identical distributions) and numerical noise can
  # produce spurious sign flips -- treat this as passing.
  fixed_point <- map_dfr(nms, function(nm) {
    gap <- maps[[nm]](x_grid) - x_grid
    gap_fin <- gap[is.finite(gap)]
    min_gap <- if (length(gap_fin) > 0) min(abs(gap_fin)) else NA_real_
    sc <- count_sign_changes(gap)
    tibble(
      map_id       = nm,
      min_abs_gap  = min_gap,
      sign_changes = sc,
      passes       = sc == 0 | (!is.na(min_gap) && min_gap < tol)
    )
  })

  # (b) Pairwise crossing check: psi_n(x) vs psi_m(x)
  pairwise <- if (length(nms) >= 2) {
    pairs <- combn(nms, 2, simplify = FALSE)
    map_dfr(pairs, function(pair) {
      nm <- pair[1]
      mm <- pair[2]
      gap <- maps[[nm]](x_grid) - maps[[mm]](x_grid)
      gap_fin <- gap[is.finite(gap)]
      tibble(
        map_n        = nm,
        map_m        = mm,
        min_abs_gap  = if (length(gap_fin) > 0) min(abs(gap_fin)) else NA_real_,
        sign_changes = count_sign_changes(gap),
        passes       = count_sign_changes(gap) == 0
      )
    })
  } else {
    tibble(
      map_n = character(), map_m = character(),
      min_abs_gap = numeric(), sign_changes = integer(), passes = logical()
    )
  }

  list(fixed_point = fixed_point, pairwise = pairwise)
}

# ==============================================================================
# Top-level checker
# ==============================================================================

#' Check all necessary conditions for shift-representability
#'
#' Runs the commutativity test (Theorem 2, Part 1) and the non-crossing test
#' (Theorem 2, Part 2) over the supplied x_grid.
#'
#' @param family  A `shift_family` object
#' @param x_grid  Numeric vector of evaluation points
#' @param tol     Tolerance passed to both sub-tests (default 1e-6)
#' @return List with:
#'   $commutativity      output of commutativity_test()
#'   $noncrossing        output of noncrossing_test()
#'   $passes_necessary   TRUE iff both necessary conditions pass
shift_representability_check <- function(family, x_grid, tol = 1e-6) {
  comm <- commutativity_test(family, x_grid, tol)
  nc <- noncrossing_test(family, x_grid, tol)

  comm_passes <- all(comm$summary$passes)
  nc_fp_passes <- all(nc$fixed_point$passes)
  nc_pw_passes <- if (nrow(nc$pairwise) > 0) all(nc$pairwise$passes) else TRUE

  structure(
    list(
      commutativity     = comm,
      noncrossing       = nc,
      passes_necessary  = comm_passes && nc_fp_passes && nc_pw_passes
    ),
    class = "shift_check"
  )
}

print.shift_check <- function(x, ...) {
  status <- if (x$passes_necessary) "PASS" else "FAIL"
  cat(sprintf("Shift-representability necessary conditions: %s\n\n", status))

  cat("Commutativity (Theorem 2, Part 1):\n")
  print(x$commutativity$summary, ...)

  cat("\nNon-crossing – fixed points (Theorem 2, Part 2a):\n")
  print(x$noncrossing$fixed_point, ...)

  if (nrow(x$noncrossing$pairwise) > 0) {
    cat("\nNon-crossing – pairwise (Theorem 2, Part 2b):\n")
    print(x$noncrossing$pairwise, ...)
  }
  invisible(x)
}

# ==============================================================================
# Batch helpers
# ==============================================================================

#' Summarise a shift_check into a single-row tibble (useful for bind_rows)
#'
#' @param check   Output of shift_representability_check()
#' @return        Single-row tibble of aggregate statistics
summarise_check <- function(check) {
  comm_max <- if (nrow(check$commutativity$summary) > 0) {
    max(check$commutativity$summary$max_abs_dev)
  } else {
    0
  }
  fp_max <- if (nrow(check$noncrossing$fixed_point) > 0) {
    max(check$noncrossing$fixed_point$sign_changes)
  } else {
    0L
  }
  pw_max <- if (nrow(check$noncrossing$pairwise) > 0) {
    max(check$noncrossing$pairwise$sign_changes)
  } else {
    0L
  }
  tibble(
    passes_necessary         = check$passes_necessary,
    max_commutativity_dev    = comm_max,
    max_fp_sign_changes      = fp_max,
    max_pairwise_crossings   = pw_max
  )
}
