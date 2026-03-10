# Transport-invariance identification from exemplar CDFs
#
# Core functions implementing the invariance criterion and optimization
# described in meta/transport-invariance-identification.md

#' Compute the commutator H_theta ∘ phi_i - phi_i ∘ H_theta evaluated on a mesh
#'
#' @param theta Parameter vector for the candidate family
#' @param u_mesh Numeric vector of mesh points in (0,1)
#' @param x_mesh Precomputed Q_1(u_mesh)
#' @param phi_values List of precomputed phi_i(u_mesh) for i = 2,...,k
#' @param exemplar_cdfs List of CDF functions for exemplars (including reference)
#' @param Q_1 Quantile function for reference
#' @param candidate_cdf_factory Function(theta) -> CDF function
#' @return List of numeric vectors, one per non-reference exemplar, each of
#'   length length(u_mesh) giving the commutator residuals
compute_commutator <- function(theta,
                               u_mesh,
                               x_mesh,
                               phi_values,
                               exemplar_cdfs,
                               Q_1,
                               candidate_cdf_factory) {
  F_theta <- candidate_cdf_factory(theta)
  H_theta_u <- F_theta(x_mesh)

  lapply(seq_along(phi_values), function(i) {
    phi_i_u <- phi_values[[i]]
    F_i <- exemplar_cdfs[[i + 1]]
    lhs <- F_theta(Q_1(phi_i_u))
    rhs <- F_i(Q_1(H_theta_u))
    lhs - rhs
  })
}

#' Compute the transport-invariance objective J(theta) over all exemplars
#'
#' @param theta Parameter vector for the candidate family
#' @param u_mesh Numeric vector of mesh points in (0,1)
#' @param x_mesh Precomputed Q_1(u_mesh)
#' @param phi_values List of precomputed phi_i(u_mesh) for i = 2,...,k
#' @param exemplar_cdfs List of CDF functions for exemplars (including reference)
#' @param Q_1 Quantile function for reference
#' @param candidate_cdf_factory Function(theta) -> CDF function
#' @param weights Quadrature weights (default: uniform)
#' @param loss Loss function (default: squared)
#' @return Scalar objective value
compute_objective <- function(theta,
                              u_mesh,
                              x_mesh,
                              phi_values,
                              exemplar_cdfs,
                              Q_1,
                              candidate_cdf_factory,
                              weights = NULL,
                              loss = function(z) z^2) {
  if (is.null(weights)) weights <- rep(1, length(u_mesh))

  residual_list <- compute_commutator(
    theta, u_mesh, x_mesh, phi_values, exemplar_cdfs, Q_1, candidate_cdf_factory
  )

  total_loss <- sapply(residual_list, \(x) sum(weights * loss(x)))
  sum(total_loss)
}

#' Plot the commutator profile for a candidate distribution
#'
#' Visualizes H_theta ∘ phi_i(u) - phi_i ∘ H_theta(u) as a function of u for
#' each non-reference exemplar. A flat line at zero indicates transport
#' invariance.
#'
#' @param theta Parameter vector for the candidate
#' @param exemplar_cdfs List of CDF functions (first is reference)
#' @param Q_1 Quantile function of the reference
#' @param candidate_cdf_factory Function(theta) -> CDF function
#' @param m Number of mesh points
#' @param eps Boundary padding
#' @param ... Additional arguments passed to plot()
#' @examples
#' # Two gamma exemplars with same shape (3) but different rates
#' exemplar_cdfs <- list(
#'   function(x) pgamma(x, shape = 3, rate = 1),
#'   function(x) pgamma(x, shape = 3, rate = 2.5)
#' )
#' Q_1 <- function(u) qgamma(u, shape = 3, rate = 1)
#'
#' # Candidate family: Gamma(shape, rate)
#' candidate_cdf_factory <- function(theta) {
#'   function(x) pgamma(x, shape = theta[1], rate = theta[2])
#' }
#'
#' # Correct shape -> residuals near zero
#' plot_commutator(c(3, 1.5), exemplar_cdfs, Q_1, candidate_cdf_factory,
#'   main = "Correct shape (transport-invariant)"
#' )
#'
#' # Wrong shape -> non-zero residuals
#' plot_commutator(c(5, 1.5), exemplar_cdfs, Q_1, candidate_cdf_factory,
#'   main = "Wrong shape (not transport-invariant)"
#' )
plot_commutator <- function(theta,
                            exemplar_cdfs,
                            Q_1,
                            candidate_cdf_factory,
                            m = 1000,
                            eps = 1e-3,
                            ...) {
  pre <- precompute(exemplar_cdfs, Q_1, m = m, eps = eps)
  residual_list <- compute_commutator(
    theta, pre$u_mesh, pre$x_mesh, pre$phi_values,
    exemplar_cdfs, Q_1, candidate_cdf_factory
  )

  k <- length(residual_list)
  y_range <- range(unlist(residual_list))
  colors <- if (k <= 8) seq_len(k) + 1 else rainbow(k)

  plot(pre$u_mesh, residual_list[[1]],
    type = "l", col = colors[1], lwd = 2,
    xlab = "u", ylab = "commutator residual",
    ylim = y_range, ...
  )
  if (k > 1) {
    for (i in 2:k) {
      lines(pre$u_mesh, residual_list[[i]], col = colors[i], lwd = 2)
    }
  }
  abline(h = 0, lty = 2, col = "gray")
  legend("topright", paste("exemplar", seq_len(k) + 1),
    col = colors, lwd = 2
  )
}

#' Plot a heatmap of the transport-invariance objective over a 2D parameter grid
#'
#' Evaluates J(theta) on a grid defined by two parameter vectors and displays
#' the result as a log-scaled heatmap using [image()]. Useful for visualizing
#' the objective landscape and confirming that invariant members form a valley
#' or ridge in parameter space.
#'
#' @param theta1_grid Numeric vector of values for the first parameter
#' @param theta2_grid Numeric vector of values for the second parameter
#' @param exemplar_cdfs List of CDF functions (first is reference)
#' @param Q_1 Quantile function of the reference
#' @param candidate_cdf_factory Function(theta) -> CDF function
#' @param param_names Character vector of length 2 for axis labels
#' @param m Number of mesh points
#' @param eps Boundary padding
#' @param loss Loss function
#' @param ... Additional arguments passed to image()
#' @return Invisibly returns the objective matrix (rows = theta1, cols = theta2)
plot_objective_heatmap <- function(theta1_grid,
                                   theta2_grid,
                                   exemplar_cdfs,
                                   Q_1,
                                   candidate_cdf_factory,
                                   param_names = c("theta1", "theta2"),
                                   m = 1000,
                                   eps = 1e-3,
                                   loss = function(z) z^2,
                                   overlay = NULL,
                                   ...) {
  pre <- precompute(exemplar_cdfs, Q_1, m = m, eps = eps)

  obj_matrix <- outer(theta1_grid, theta2_grid, Vectorize(function(t1, t2) {
    compute_objective(
      theta = c(t1, t2),
      u_mesh = pre$u_mesh,
      x_mesh = pre$x_mesh,
      phi_values = pre$phi_values,
      exemplar_cdfs = exemplar_cdfs,
      Q_1 = Q_1,
      candidate_cdf_factory = candidate_cdf_factory,
      loss = loss
    )
  }))

  log_obj <- log10(obj_matrix + .Machine$double.eps)
  palette <- hcl.colors(64, "YlOrRd", rev = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
  par(mar = c(5, 4, 4, 1))

  image(theta1_grid, theta2_grid, log_obj,
    xlab = param_names[1], ylab = param_names[2],
    col = palette, useRaster = TRUE,
    ...
  )
  # contour(theta1_grid, theta2_grid, log_obj, add = TRUE, col = "gray30")
  if (is.function(overlay)) overlay()

  # Color scale legend
  par(mar = c(5, 0.5, 4, 3))
  z_range <- range(log_obj, finite = TRUE)
  z_seq <- seq(z_range[1], z_range[2], length.out = 64)
  image(1, z_seq, matrix(z_seq, nrow = 1),
    col = palette, useRaster = TRUE,
    xaxt = "n", yaxt = "n", xlab = "", ylab = ""
  )
  axis(4, las = 1)
  mtext("log10 J", side = 4, line = 2)

  layout(1)
  par(mar = c(5, 4, 4, 2))

  invisible(obj_matrix)
}

#' Precompute mesh and transport maps for the optimization
#'
#' @param exemplar_cdfs List of CDF functions (first is reference)
#' @param Q_1 Quantile function of the reference
#' @param m Number of mesh points
#' @param eps Boundary padding in (0, 0.5)
#' @return List with u_mesh, x_mesh, phi_values
precompute <- function(exemplar_cdfs, Q_1, m = 1000, eps = 1e-2) {
  u_mesh <- seq(eps, 1 - eps, length.out = m)
  x_mesh <- Q_1(u_mesh)

  phi_values <- lapply(exemplar_cdfs[-1], function(F_i) F_i(x_mesh))

  list(u_mesh = u_mesh, x_mesh = x_mesh, phi_values = phi_values)
}

#' Search for a transport-invariant member via multi-start optimization
#'
#' @param exemplar_cdfs List of CDF functions (first element is the reference)
#' @param Q_1 Quantile function of the reference exemplar
#' @param candidate_cdf_factory Function(theta) -> CDF function
#' @param theta_lower Lower bounds for theta
#' @param theta_upper Upper bounds for theta
#' @param n_starts Number of random multi-start initializations
#' @param m Number of mesh points
#' @param eps Boundary padding
#' @param loss Loss function
#' @param optim_method Optimization method (passed to optim)
#' @return List with best theta, objective value, all results, diagnostics
find_invariant_member <- function(exemplar_cdfs,
                                  Q_1,
                                  candidate_cdf_factory,
                                  theta_lower,
                                  theta_upper,
                                  n_starts = 20,
                                  m = 1000,
                                  eps = 1e-3,
                                  loss = function(z) z^2,
                                  optim_method = "L-BFGS-B") {
  pre <- precompute(exemplar_cdfs, Q_1, m = m, eps = eps)

  obj_fn <- function(theta) {
    compute_objective(
      theta = theta,
      u_mesh = pre$u_mesh,
      x_mesh = pre$x_mesh,
      phi_values = pre$phi_values,
      exemplar_cdfs = exemplar_cdfs,
      Q_1 = Q_1,
      candidate_cdf_factory = candidate_cdf_factory,
      loss = loss
    )
  }

  d <- length(theta_lower)

  if (d == 1L) {
    opt <- optimize(obj_fn,
      lower = theta_lower, upper = theta_upper,
      tol = .Machine$double.eps^0.75
    )
    best <- list(
      par = opt$minimum, value = opt$objective, convergence = 0L,
      message = "optimize (Brent)"
    )
    results <- list(best)
  } else {
    results <- vector("list", n_starts)
    for (s in seq_len(n_starts)) {
      theta_init <- theta_lower + runif(d, 0.01, 1) * (theta_upper - theta_lower)
      res <- tryCatch(
        optim(
          par = theta_init,
          fn = obj_fn,
          method = optim_method,
          lower = theta_lower,
          upper = theta_upper,
          control = list(
            maxit = 1000,
            factr = 1e-12,
            pgtol = 1e-15,
            ndeps = rep(1e-6, d)
          )
        ),
        error = function(e) list(par = theta_init, value = Inf, convergence = -1)
      )
      results[[s]] <- res
    }
  }

  objectives <- vapply(results, function(r) r$value, numeric(1))
  best_idx <- which.min(objectives)
  best <- results[[best_idx]]

  # Uniqueness diagnostic: how many near-optimal solutions converge to same theta
  tol_obj <- max(best$value * 10, 1e-10)
  near_optimal <- which(objectives <= tol_obj)
  theta_matrix <- do.call(rbind, lapply(results[near_optimal], function(r) r$par))
  theta_spread <- if (length(near_optimal) > 1) apply(theta_matrix, 2, sd) else rep(0, d)

  list(
    theta_star = best$par,
    objective = best$value,
    convergence = best$convergence,
    all_objectives = objectives,
    near_optimal_theta_sd = theta_spread,
    n_near_optimal = length(near_optimal),
    all_results = results,
    precomputed = pre
  )
}
