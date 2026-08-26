# =============================================================================
# DimSEM_ridgeReport() -- corrected
# =============================================================================
#
# BUG IN THE PREVIOUS VERSION. `ridge_share` was computed as
#
#     pmin(1, (delta_k * sd_c / se(beta_k))^2),
#     sd_c = (1/s_alpha^2 + sum(delta^2)/s_beta^2)^(-1/2)
#
# where sd_c is the posterior SD of the FLAT DIRECTION -- the alpha/beta ridge.
# That quantity exists only when the decomposition is left unconstrained. In a
# sum-to-zero fit there is no flat direction, so sd_c is a counterfactual: the
# width the ridge WOULD have had. Dividing the observed SE into it and capping
# at 1 then inverts the meaning entirely -- a "share" of 1.0000 signals that
# the observed SE is SMALLER than the counterfactual ridge, which is the
# constraint working, not failing.
#
# THIS VERSION keys the report to the identification actually used:
#
#   prior        -> ridge_share as before. The flat direction is real, the
#                   priors set its width, and the share is the fraction of the
#                   coefficient's posterior variance attributable to it.
#
#   sum_to_zero  -> no share is reported at all. Instead `ridge_avoided` gives
#                   sd_c(counterfactual) / se(beta), i.e. how many times
#                   NARROWER the constrained coefficient is than an
#                   unconstrained one would have been. Values above 1 quantify
#                   the benefit; values below 1 simply mean that coefficient
#                   had little ridge contamination to begin with (which is
#                   itself informative -- it happens where delta_k is small).
#
# In both cases the metric-invariant quantity se(beta)/se(total) is reported,
# since the totals are what the likelihood identifies regardless of convention.
# =============================================================================

DimSEM_ridgeReport <- function(object, s_alpha = 0.5, s_beta = 0.5) {
  stopifnot(inherits(object, "DimSEM_estimate"))
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior is required.", call. = FALSE)
  }
  map   <- object$param_map
  unity <- identical(map$hyper_loadings, "unity")
  facs  <- map$factors
  m     <- length(facs)

  # Which identification was actually used? Fall back to "prior" only if the
  # fit predates the field, and say so in the attribute.
  strategy <- object$covariate_identification %||%
    object$stan$identification$strategy %||% NA_character_
  constrained <- isTRUE(strategy == "sum_to_zero")

  gm <- function(v) {
    d <- tryCatch(object$fit$draws(variables = v), error = function(e) NULL)
    if (is.null(d)) stop("`", v, "` is not available in this fit.",
                         call. = FALSE)
    posterior::as_draws_matrix(d)
  }

  bf <- gm("beta_f_cor"); bg <- gm("beta_g_cor"); bt <- gm("beta_total")
  n  <- nrow(bf)
  beta  <- vapply(seq_len(m), function(k) bf[, k], numeric(n))
  tot   <- vapply(seq_len(m), function(k) bt[, k], numeric(n))
  alpha <- as.numeric(bg[, 1])

  # delta: 1 under unity; otherwise recovered draw-wise from the identity
  # beta_total = beta_f_cor + delta * beta_g_cor, which the generated
  # quantities block enforces by construction.
  delta <- if (unity) matrix(1, n, m) else sweep(tot - beta, 1, alpha, "/")

  d_bar <- colMeans(delta)
  sd_c  <- (1 / s_alpha^2 + sum(d_bar^2) / s_beta^2)^(-0.5)
  se_b  <- c(apply(beta, 2, stats::sd), stats::sd(alpha))
  ridge <- c(d_bar * sd_c, sd_c)
  se_t  <- c(apply(tot, 2, stats::sd), NA_real_)

  out <- data.frame(
    parameter = c(facs, "G1"),
    est       = c(colMeans(beta), mean(alpha)),
    se        = se_b,
    se_total  = se_t,
    se_ratio  = se_b / se_t,
    cor_alpha_beta = c(vapply(seq_len(m), function(k)
                         stats::cor(alpha, beta[, k]), numeric(1)), NA_real_),
    row.names = NULL
  )

  if (constrained) {
    # No flat direction exists; report the width that was avoided instead.
    out$ridge_sd_counterfactual <- ridge
    out$ridge_avoided <- ridge / se_b
    attr(out, "note") <- paste0(
      "Identification: sum_to_zero. There is no alpha/beta ridge in this ",
      "fit, so no ridge share is reported. `ridge_sd_counterfactual` is the ",
      "width the flat direction would have had WITHOUT the constraint, and ",
      "`ridge_avoided` is that width divided by the observed SE: values > 1 ",
      "measure how much narrower the constrained coefficient is. Values < 1 ",
      "mean that coefficient carried little ridge contamination anyway, ",
      "which happens where delta_k is small.")
  } else {
    out$ridge_sd    <- ridge
    out$ridge_share <- pmin(1, (ridge / se_b)^2)
    attr(out, "note") <- paste0(
      "Identification: ", strategy, ". The alpha/beta split is left to the ",
      "priors, so the flat direction is real and `ridge_share` is the ",
      "fraction of each coefficient's posterior variance contributed by it. ",
      "Shares near 1 mean the interval is the prior's width, not the data's.")
  }
  attr(out, "strategy") <- strategy
  attr(out, "delta_posterior_mean") <- d_bar
  out
}

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x


# =============================================================================
# USAGE
# =============================================================================
# source("DimSEM_ridgeReport.R")
#
# rC <- DimSEM_ridgeReport(fit_wstz)                       # sum_to_zero
# print(rC, row.names = FALSE, digits = 3); cat(attr(rC, "note"), "\n")
#
# rA <- DimSEM_ridgeReport(readRDS("dimsem_32_party_free_prior.rds"))   # prior
# print(rA, row.names = FALSE, digits = 3)
#
# The pair is the argument in one table: run A's `ridge_share` shows how much
# of each interval was prior width; run C's `ridge_avoided` shows how much of
# that width the constraint removed. Both should be largest where delta_k is
# largest, because the ridge displaces each coefficient by delta_k * c.
