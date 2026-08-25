# =============================================================================
# DimSEM_standardized() and DimSEM_ridgeReport()
# =============================================================================
#
# Two things DimSEM_parameterEstimates() does not currently give you:
#
#   DimSEM_standardized()  -- draw-wise standardized coefficients, on the
#                             same footing as blavaan's Std.all column
#   DimSEM_ridgeReport()   -- how much of each coefficient's posterior width
#                             comes from the alpha/beta ridge rather than
#                             from the data
#
# Standardization in DimSEM's parameterization. The engine fixes each facet's
# disturbance variance to 1 and lets the loadings carry the scale, so
#
#     G      = X*alpha + tau_g * u,        var(u) = 1
#     eta_k  = X*beta_k + delta_k*G + z_k, var(z_k) = 1
#     var(G)     = alpha^2 var(X) + tau_g^2
#     var(eta_k) = t_k^2 var(X) + delta_k^2 tau_g^2 + 1,   t_k = beta_k + delta_k*alpha
#
# with delta_k == 1 under hyper_loadings = "unity" and tau_g == 1 under
# "free". Standardized quantities are then
#
#     beta_k^std  = beta_k * sd(X) / sd(eta_k)
#     alpha^std   = alpha  * sd(X) / sd(G)
#     delta_k^std = delta_k * sd(G) / sd(eta_k)
#     t_k^std     = t_k    * sd(X) / sd(eta_k)   ( = beta^std + delta^std*alpha^std )
#
# All computed PER DRAW, so the reported intervals are correct rather than
# plug-ins on posterior means.
#
# REQUIRES the fit to still carry its draws. cmdstanr reads draws lazily from
# temporary CSV files; a fit that was saveRDS()'d and reloaded in a new
# session retains only the variables that were touched before saving. If
# these functions report missing variables, that is the cause -- refit with a
# current version of the package (which now materializes the draws before
# returning) or call fit$draws() before saving.
# =============================================================================


.dsx_draws <- function(object, vars) {
  d <- tryCatch(object$fit$draws(variables = vars), error = function(e) NULL)
  if (is.null(d)) {
    stop("Variable(s) not available in this fit: ", paste(vars, collapse = ", "),
         ".\nThe draws were probably lost to cmdstanr's lazy CSV loading when ",
         "the object was saved and reloaded. Refit, or call fit$draws() before ",
         "saveRDS().", call. = FALSE)
  }
  m <- posterior::as_draws_matrix(d)
  m
}

DimSEM_standardized <- function(object, data, prob = 0.95) {
  stopifnot(inherits(object, "DimSEM_estimate"))
  .dimsem_require_namespace <- get(".dimsem_require_namespace",
                                   envir = asNamespace("DimSEM"))
  .dimsem_require_namespace("posterior")

  map    <- object$param_map
  unity  <- identical(map$hyper_loadings, "unity")
  covs   <- object$covariates
  facs   <- map$factors
  m      <- length(facs)
  if (is.null(covs)) stop("No covariates in this fit.", call. = FALSE)

  # Covariate SD on the scale the model actually saw. If the fit standardized
  # its covariates, that scale is 1.
  sdx <- vapply(covs, function(cv) {
    if (isTRUE(object$covariate_scaling$standardized[[cv]] %||% FALSE)) 1
    else stats::sd(data[[cv]], na.rm = TRUE)
  }, numeric(1))

  bf <- .dsx_draws(object, "beta_f_cor")     # Q x K, column-major
  bg <- .dsx_draws(object, "beta_g_cor")     # Q x G
  tg <- if (unity) .dsx_draws(object, "tau_g") else NULL
  dl <- if (unity) NULL else .dsx_draws(object, "delta")
  n  <- nrow(bf)

  idx <- function(mat, q, k, nk) mat[, (k - 1L) * length(covs) + q]

  out <- list()
  for (q in seq_along(covs)) {
    beta  <- vapply(seq_len(m), function(k) idx(bf, q, k, m), numeric(n))
    alpha <- idx(bg, q, 1L, 1L)
    delta <- if (unity) matrix(1, n, m) else
      vapply(seq_len(m), function(k) dl[, k], numeric(n))
    tau   <- if (unity) as.numeric(tg[, 1]) else rep(1, n)

    tot   <- beta + delta * alpha
    vx    <- sdx[q]^2
    v_eta <- tot^2 * vx + delta^2 * tau^2 + 1
    v_G   <- alpha^2 * vx + tau^2

    b_std <- beta * sdx[q] / sqrt(v_eta)
    a_std <- alpha * sdx[q] / sqrt(v_G)
    d_std <- delta * sqrt(v_G) / sqrt(v_eta)
    t_std <- tot  * sdx[q] / sqrt(v_eta)

    qs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
    smry <- function(x, lab) data.frame(
      parameter = lab, est = mean(x), se = stats::sd(x),
      ci.lower = unname(stats::quantile(x, qs[1])),
      ci.upper = unname(stats::quantile(x, qs[2])),
      row.names = NULL, stringsAsFactors = FALSE)

    rows <- rbind(
      smry(a_std, paste0("G1 ~ ", covs[q], " (std)")),
      do.call(rbind, lapply(seq_len(m), function(k)
        smry(b_std[, k], paste0(facs[k], " ~ ", covs[q], " (std)")))),
      do.call(rbind, lapply(seq_len(m), function(k)
        smry(t_std[, k], paste0(facs[k], " ~ ", covs[q], " (std total)")))),
      do.call(rbind, lapply(seq_len(m), function(k)
        smry(d_std[, k], paste0("G1 =~ ", facs[k], " (std)"))))
    )
    out[[covs[q]]] <- rows
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}


# --- How much of each posterior width is the ridge? ---------------------------
#
# The flat direction is  alpha -> alpha + c,  beta_k -> beta_k - delta_k*c.
# Its posterior SD is set by the priors alone:
#     sd(c) = (1/s_alpha^2 + sum_k delta_k^2 / s_beta^2)^(-1/2)
# and it contributes delta_k * sd(c) to each beta_k and sd(c) to alpha. The
# ridge SHARE below is that contribution as a fraction of the observed
# posterior variance: near 1 means the coefficient's uncertainty is the prior
# talking, near 0 means it is the data.
#
# Also reported: corr(alpha, beta_k) across draws -- strongly negative is the
# ridge's signature -- and se(beta)/se(total), since the totals are what the
# likelihood identifies.
DimSEM_ridgeReport <- function(object, s_alpha = 0.5, s_beta = 0.5) {
  stopifnot(inherits(object, "DimSEM_estimate"))
  map   <- object$param_map
  unity <- identical(map$hyper_loadings, "unity")
  facs  <- map$factors
  m     <- length(facs)

  bf <- .dsx_draws(object, "beta_f_cor")
  bg <- .dsx_draws(object, "beta_g_cor")
  bt <- .dsx_draws(object, "beta_total")
  dl <- if (unity) NULL else .dsx_draws(object, "delta")

  beta  <- vapply(seq_len(m), function(k) bf[, k], numeric(nrow(bf)))
  tot   <- vapply(seq_len(m), function(k) bt[, k], numeric(nrow(bt)))
  alpha <- bg[, 1]
  delta <- if (unity) matrix(1, nrow(bf), m) else
    vapply(seq_len(m), function(k) dl[, k], numeric(nrow(dl)))

  d_bar <- colMeans(delta)
  sd_c  <- (1 / s_alpha^2 + sum(d_bar^2) / s_beta^2)^(-0.5)

  data.frame(
    parameter   = c(facs, "G1"),
    est         = c(colMeans(beta), mean(alpha)),
    se          = c(apply(beta, 2, stats::sd), stats::sd(alpha)),
    se_total    = c(apply(tot, 2, stats::sd), NA_real_),
    ridge_sd    = c(d_bar * sd_c, sd_c),
    ridge_share = pmin(1, (c(d_bar * sd_c, sd_c) /
                           c(apply(beta, 2, stats::sd), stats::sd(alpha)))^2),
    cor_alpha   = c(vapply(seq_len(m), function(k)
                      stats::cor(alpha, beta[, k]), numeric(1)), NA_real_),
    row.names = NULL
  )
}


# =============================================================================
# USAGE
# =============================================================================
# source("DimSEM_standardized.R")
#
# std_free <- DimSEM_standardized(fit_free, dat_ds)
# std_stz  <- DimSEM_standardized(fit_stz,  dat_ds)
# print(std_free, row.names = FALSE, digits = 3)
#
# DimSEM_ridgeReport(fit_free)   # prior identification: expect high shares
# DimSEM_ridgeReport(fit_stz)    # sum-to-zero: shares should be ~0
#
# READING IT
#   ridge_share near 1  -> that coefficient's interval is the prior's width,
#                          not the data's. The total effect is still sound.
#   ridge_share near 0  -> the coefficient is data-identified.
#   cor_alpha strongly negative -> alpha and beta trade off one-for-one, the
#                          ridge's signature.
