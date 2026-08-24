# =============================================================================
# blavaan cross-check for the DimSEM Bayesian hierarchical model
# =============================================================================
#
# Purpose. Estimate the SAME hierarchical measurement + structural model that
# DimSEM_estimate(engine = "bayes", hyper_loadings = "unity") fits, but with
# the canned blavaan/bsem machinery used in the original replication script,
# so that the two estimators can be compared on identical structure.
#
# What is held fixed vs. the original replication code
# ----------------------------------------------------
#   HELD FIXED (mirrored):
#     * measurement syntax        Fat_k =~ item + item + ...
#     * hyper syntax              ideology =~ a*Fat_1 + a*Fat_2 + ...
#                                 (equality-constrained loadings + free hyper
#                                 variance; this is an exact reparameterization
#                                 of DimSEM's hyper_loadings = "unity")
#     * structural syntax         Fat_k ~ prior('normal(0, 0.5)')*<cov>
#                                 ideology ~ prior('normal(0, 0.5)')*<cov>
#     * dpriors()                 nu / alpha / lambda / beta = normal(0, 2)
#     * bsem() call               n.chains, burnin, sample, adapt_delta,
#                                 max_treedepth, cores, mcmcfile, fixed.x,
#                                 orthogonal
#     * the two vestigial lavaan DWLS fits are reproduced (see note below)
#
#   CHANGED, deliberately:
#     * The partition is taken from a DimSEM_propose() object instead of being
#       re-derived by the inline EGA + resolution-annealing loop. That is the
#       entire point of the comparison: same items, same communities, same
#       model -- different estimator. Set `refit_ega = TRUE` to run the
#       original loop instead and confirm the two partitions agree.
#
# Note on the original code's lavaan step. In the replication script the two
# lavaan::sem() DWLS fits extract `pars_base` / `pars_hyper`, but both objects
# are immediately overwritten with zeros, so only their DIMENSION survives (it
# sets how many prior strings are pasted). The fits are therefore vestigial:
# every prior is centred at 0 regardless. They are reproduced here for
# fidelity and can be skipped with `run_lavaan_starts = FALSE` -- the fitted
# model is bit-identical either way.
#
# Requires: blavaan, lavaan, DimSEM, posterior. EGAnet and dplyr only if
# refit_ega = TRUE.
# =============================================================================


# --- 1. Build the model syntax from a DimSEM proposal -------------------------

dimsem_to_blavaan_syntax <- function(proposal, covariate,
                                     hyper_name = "ideology",
                                     factor_prefix = "Fat_",
                                     beta_prior = "normal(0, 0.5)",
                                     equality_constrained = TRUE) {
  stopifnot(inherits(proposal, "DimSEM_proposal"))
  tab <- proposal$selected$partition_table
  tab <- tab[!is.na(tab$dimension), , drop = FALSE]   # drop audited orphans
  if (nrow(tab) == 0L) stop("Proposal contains no assigned items.")

  dims <- sort(unique(tab$dimension))
  measure <- lapply(dims, function(d) tab$item[tab$dimension == d])
  # Names mirror the replication script (Fat_1, Fat_2, ...); the mapping to
  # DimSEM's own labels is Fat_<d> <-> F<d>, i.e. identical ordering.
  names(measure) <- paste0(factor_prefix, dims)

  fat_measure <- paste(
    vapply(seq_along(measure), function(i) {
      paste0(names(measure)[i], " =~ ", paste(measure[[i]], collapse = " + "))
    }, character(1)),
    collapse = " \n "
  )

  # ideology =~ a*Fat_1 + a*Fat_2 + ...  -- equality-constrained loadings with
  # a freely estimated hyper variance. Only the product (loading x hyper SD)
  # enters the covariance structure, so this is the same model as fixing all
  # loadings to 1 and freeing the variance (DimSEM's "unity").
  rhs <- if (equality_constrained) paste0("a*", names(measure)) else names(measure)
  hyper_measure <- paste0(hyper_name, " =~ ", paste(rhs, collapse = " + "))
  model_measure <- paste(fat_measure, " \n ", hyper_measure)

  # Structural part: covariate regressed on every facet AND on the hyper
  # factor simultaneously -- the configuration whose alpha/beta split is
  # identified only by the priors under this parameterization.
  prior_tag <- paste0("prior('", beta_prior, "')*")
  fat_struct <- paste0(names(measure), " ~ ", prior_tag, covariate,
                       collapse = " \n ")
  hyper_struct <- paste0(hyper_name, " ~ ", prior_tag, covariate)
  model_bayes <- paste(model_measure, " \n ", fat_struct, " \n ", hyper_struct)

  list(model = model_bayes,
       model_measure = model_measure,
       measure = measure,
       factors = names(measure),
       hyper_name = hyper_name,
       covariate = covariate,
       items = tab$item,
       n_struct = length(measure) + 1L)
}


# --- 2. Assemble the analysis frame ------------------------------------------

# DimSEM z-scores continuous covariates internally (standardize_covariates =
# TRUE) while the replication script feeds the raw 0-1 range_0_1() version.
# The two differ by a constant factor, so coefficients are NOT on the same
# scale unless matched. standardize_covariate = TRUE reproduces DimSEM's
# metric; set FALSE to reproduce the original script exactly.
build_analysis_frame <- function(items_df, covariate_values,
                                 covariate_name = "cov",
                                 item_names,
                                 standardize_covariate = TRUE) {
  if (nrow(items_df) != length(covariate_values)) {
    stop("items_df and covariate_values have different row counts (",
         nrow(items_df), " vs ", length(covariate_values),
         "); align them before calling.")
  }
  x <- as.numeric(covariate_values)
  if (isTRUE(standardize_covariate)) {
    x <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }
  out <- items_df[, item_names, drop = FALSE]
  out[[covariate_name]] <- x
  # Mirrors `data <- data[!apply(is.na(data), 1, all), ]`: drop rows that are
  # entirely missing. Rows with SOME missingness are retained (see `missing=`).
  keep <- !apply(is.na(out), 1, all)
  out[keep, , drop = FALSE]
}


# --- 3. Fit -------------------------------------------------------------------

fit_blavaan_hierarchical <- function(syntax, data,
                                     run_lavaan_starts = TRUE,
                                     n.chains = 20, burnin = 600, sample = 400,
                                     cores = 20, adapt_delta = 0.9,
                                     max_treedepth = 20, mcmcfile = TRUE,
                                     missing = "ml", seed = 123) {
  if (!requireNamespace("blavaan", quietly = TRUE)) {
    stop("blavaan is required.")
  }

  # Vestigial DWLS fits from the original script (see header note). Wrapped in
  # tryCatch because DWLS on 50 ordinal-ish items can fail to converge without
  # affecting anything downstream.
  if (isTRUE(run_lavaan_starts)) {
    invisible(tryCatch({
      lavaan::sem(model = paste(syntax$model_measure, " \n ",
                                paste0(syntax$factors, " ~ ",
                                       syntax$covariate, collapse = " \n ")),
                  data = data, orthogonal = FALSE, missing = "pairwise",
                  estimator = "DWLS", optim.dx.tol = 0.01, fixed.x = FALSE,
                  optim.method = "GN", optim.gn.tol.x = 1e-04,
                  optim.gn.iter.max = 1000)
    }, error = function(e) {
      message("lavaan start-value fit failed (harmless): ",
              conditionMessage(e))
      NULL
    }))
  }

  new_prior <- blavaan::dpriors(nu = "normal(0, 2)", alpha = "normal(0, 2)",
                                lambda = "normal(0, 2)", beta = "normal(0, 2)")

  blavaan::bsem(
    model = syntax$model, data = data,
    fixed.x = FALSE, orthogonal = FALSE, missing = missing,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    bcontrol = list(cores = cores),
    dp = new_prior, mcmcfile = mcmcfile,
    n.chains = n.chains, burnin = burnin, sample = sample, seed = seed
  )
}


# --- 4. Extract the covariate decomposition on a comparable footing ----------

# Locate one parameter among blavaan's standardized-draw column names.
#
# Names are the lavaan partable form pasted together, e.g. "ideology=~Fat_1",
# "Fat_1~lib_con". Matching is done by EXACT string comparison rather than by
# regex: the operators "=~", "~~" and "~" contain characters that are regex
# metacharacters, and escaping them by hand is exactly the sort of thing that
# breaks (a bracket expression does not treat backslash as an escape in R's
# default TRE engine, so a hand-rolled class silently mis-parses).
#
# Three passes, each strictly less literal than the last:
#   1. exact match
#   2. whitespace-insensitive ("Fat_1 ~ lib_con" vs "Fat_1~lib_con")
#   3. punctuation-insensitive, which recovers names mangled by
#      make.names() -- e.g. a table round-tripped through read.csv(), where
#      "ideology=~Fat_1" becomes "ideology..Fat_1"
# An ambiguous match is an error, never a silent first-hit.
.dimsem_match_par <- function(nm, lhs, op, rhs) {
  target <- paste0(lhs, op, rhs)
  i <- which(nm == target)
  if (length(i) == 1L) return(i)

  no_ws <- function(x) gsub("[[:space:]]+", "", x)
  i <- which(no_ws(nm) == no_ws(target))
  if (length(i) == 1L) return(i)

  key <- function(x) gsub("[^[:alnum:]_]", "", x)
  i <- which(key(nm) == key(target))
  if (length(i) == 1L) return(i)

  if (length(i) > 1L) {
    stop("Ambiguous match for '", target, "': ",
         paste(nm[i], collapse = ", "), call. = FALSE)
  }
  near <- nm[grepl(lhs, nm, fixed = TRUE) | grepl(rhs, nm, fixed = TRUE)]
  stop("Could not locate '", target, "' among the standardized draws.\n",
       if (length(near)) paste0("Names mentioning '", lhs, "' or '", rhs,
                                "': ", paste(near, collapse = ", "))
       else paste0("First 20 available names: ",
                   paste(utils::head(nm, 20), collapse = ", ")),
       call. = FALSE)
}

# Inspect the available names directly when a lookup fails.
blavaan_draw_names <- function(fit) {
  colnames(blavaan::standardizedPosterior(fit))
}


# Returns, per facet: the hyper coefficient (alpha), the facet coefficients
# (beta_k), the LIKELIHOOD-IDENTIFIED totals t_k = beta_k + delta_k * alpha,
# and the totals re-split under DimSEM's sum-to-zero convention. Standardized
# draws are used so the two estimators' differing scale conventions (blavaan:
# marker item; DimSEM: unit disturbance variance) do not confound the
# comparison.
blavaan_covariate_decomposition <- function(fit, syntax, prob = 0.95) {
  if (!requireNamespace("blavaan", quietly = TRUE)) stop("blavaan required.")
  # Keep the matrix: as.data.frame() is name-safe here, but working with the
  # matrix avoids any downstream make.names() surprises.
  draws <- blavaan::standardizedPosterior(fit)
  nm <- colnames(draws)
  hyp <- syntax$hyper_name
  fac <- syntax$factors
  cov <- syntax$covariate

  col <- function(lhs, op, rhs) draws[, .dimsem_match_par(nm, lhs, op, rhs)]

  alpha <- col(hyp, "~", cov)
  beta <- vapply(fac, function(f) col(f, "~", cov), numeric(nrow(draws)))
  delta <- vapply(fac, function(f) col(hyp, "=~", f), numeric(nrow(draws)))
  beta <- matrix(beta, nrow = nrow(draws), dimnames = list(NULL, fac))
  delta <- matrix(delta, nrow = nrow(draws), dimnames = list(NULL, fac))

  totals <- beta + delta * alpha

  # DimSEM's convention, applied draw-wise to the same posterior: the
  # (loading-weighted) projection onto the sum-to-zero subspace.
  a_stz <- rowSums(delta * totals) / rowSums(delta^2)
  b_stz <- totals - delta * a_stz

  q <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  summarise <- function(x, label) {
    x <- as.matrix(x)
    data.frame(
      parameter = label,
      est = colMeans(x),
      se = apply(x, 2, stats::sd),
      ci.lower = apply(x, 2, stats::quantile, q[1]),
      ci.upper = apply(x, 2, stats::quantile, q[2]),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }

  out <- rbind(
    summarise(alpha, paste0(hyp, " ~ ", cov, " [as fitted]")),
    summarise(beta, paste0(fac, " ~ ", cov, " [as fitted]")),
    summarise(totals, paste0(fac, " ~ ", cov, " (total)")),
    summarise(a_stz, paste0(hyp, " ~ ", cov, " [sum-to-zero]")),
    summarise(b_stz, paste0(fac, " ~ ", cov, " [sum-to-zero]")),
    summarise(delta, paste0(hyp, " =~ ", fac))
  )

  attr(out, "ridge_diagnostics") <- data.frame(
    facet = fac,
    cor_alpha_beta = vapply(seq_along(fac), function(k)
      stats::cor(alpha, beta[, k]), numeric(1)),
    sd_beta = apply(beta, 2, stats::sd),
    sd_total = apply(totals, 2, stats::sd),
    ratio_sd = apply(beta, 2, stats::sd) / apply(totals, 2, stats::sd),
    row.names = NULL
  )
  attr(out, "sum_beta_as_fitted") <- sum(colMeans(beta))
  out
}


# --- 5. Side-by-side with the DimSEM fit -------------------------------------

compare_with_dimsem <- function(blav_decomp, dimsem_tab, covariate,
                                factor_prefix = "Fat_") {
  d <- dimsem_tab
  ds_alpha <- d[d$op == "~" & grepl("^G", d$lhs) & d$rhs == covariate, ]
  ds_beta <- d[d$op == "~" & grepl("^F[0-9]", d$lhs) & d$rhs == covariate, ]
  ds_tot <- d[d$op == ":=" & grepl("_total$", d$lhs), ]

  blav_key <- function(tag) {
    b <- blav_decomp[grepl(tag, blav_decomp$parameter, fixed = TRUE), ]
    b[order(b$parameter), ]
  }
  b_alpha <- blav_key("[sum-to-zero]")[1, ]
  b_beta <- blav_key("[sum-to-zero]")[-1, ]
  b_tot <- blav_key("(total)")

  cat("\n=== hyper-level coefficient (sum-to-zero convention) ===\n")
  cat(sprintf("  blavaan: %7.3f [%.3f, %.3f]\n", b_alpha$est,
              b_alpha$ci.lower, b_alpha$ci.upper))
  if (nrow(ds_alpha)) {
    cat(sprintf("  DimSEM : %7.3f [%.3f, %.3f]   rhat %.2f  ess %.0f\n",
                ds_alpha$est, ds_alpha$ci.lower, ds_alpha$ci.upper,
                ds_alpha$rhat, ds_alpha$ess_bulk))
  }

  cat("\n=== facet coefficients ===\n")
  n <- min(nrow(b_beta), nrow(ds_beta))
  for (i in seq_len(n)) {
    cat(sprintf("  %-6s blavaan %+7.3f [%+.3f, %+.3f] | DimSEM %+7.3f [%+.3f, %+.3f] rhat %.2f\n",
                ds_beta$lhs[i], b_beta$est[i], b_beta$ci.lower[i],
                b_beta$ci.upper[i], ds_beta$est[i], ds_beta$ci.lower[i],
                ds_beta$ci.upper[i], ds_beta$rhat[i]))
  }

  cat("\n=== total effects (the likelihood-identified quantity) ===\n")
  n <- min(nrow(b_tot), max(nrow(ds_tot), 1L))
  for (i in seq_len(nrow(b_tot))) {
    ds_line <- if (nrow(ds_tot) >= i) {
      sprintf("| DimSEM %+7.3f [%+.3f, %+.3f] rhat %.2f", ds_tot$est[i],
              ds_tot$ci.lower[i], ds_tot$ci.upper[i], ds_tot$rhat[i])
    } else {
      "| DimSEM (no *_total rows -- table predates beta_total)"
    }
    cat(sprintf("  %-8s blavaan %+7.3f [%+.3f, %+.3f] %s\n",
                paste0("t_", i), b_tot$est[i], b_tot$ci.lower[i],
                b_tot$ci.upper[i], ds_line))
  }

  cat("\n=== ridge diagnostics (blavaan posterior) ===\n")
  print(attr(blav_decomp, "ridge_diagnostics"), row.names = FALSE, digits = 3)
  cat("\nsum of facet coefficients as fitted:",
      round(attr(blav_decomp, "sum_beta_as_fitted"), 4),
      "(0 only under the sum-to-zero convention)\n")
  invisible(NULL)
}


# --- 6. Hyper-level variance: is the parent factor doing any work? -----------

# The live hypothesis for the ANES pathology is that a single compound-
# symmetric hyper variance cannot reconcile heterogeneous residual facet
# correlations, pushing tau^2 to its boundary. blavaan estimates the same
# quantity, so this is a direct cross-check that does not depend on DimSEM's
# sampler behaving well.
blavaan_hyper_strength <- function(fit, syntax) {
  draws <- blavaan::standardizedPosterior(fit)
  nm <- colnames(draws)
  delta <- vapply(syntax$factors,
                  function(f) draws[, .dimsem_match_par(nm, syntax$hyper_name,
                                                        "=~", f)],
                  numeric(nrow(draws)))
  delta <- matrix(delta, nrow = nrow(draws))
  # Standardized hyper loadings; the implied residual facet correlation for a
  # pair is delta_j * delta_k. Compound symmetry would make all of these equal.
  m <- ncol(delta)
  pairs <- utils::combn(m, 2)
  implied <- apply(pairs, 2, function(p) delta[, p[1]] * delta[, p[2]])
  data.frame(
    pair = apply(pairs, 2, function(p)
      paste(syntax$factors[p[1]], syntax$factors[p[2]], sep = "-")),
    implied_r = colMeans(implied),
    ci.lower = apply(implied, 2, stats::quantile, 0.025),
    ci.upper = apply(implied, 2, stats::quantile, 0.975),
    row.names = NULL
  )
}


# =============================================================================
# USAGE
# =============================================================================
# Assumes in memory: anes2000TS, all_items (50 policy items), bayes_2000_prop,
# and (optionally) lib_con_bayes_res_tab from the DimSEM run.
#
# range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
#                          (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
#
# syntax <- dimsem_to_blavaan_syntax(bayes_2000_prop, covariate = "lib_con")
# cat(syntax$model)                      # inspect before fitting
#
# dat <- build_analysis_frame(
#   items_df   = all_items,
#   covariate_values = range_0_1(anes2000TS$lib_con),
#   covariate_name   = "lib_con",
#   item_names = syntax$items,
#   standardize_covariate = TRUE)        # FALSE reproduces the original 0-1 scale
#
# blav_fit <- fit_blavaan_hierarchical(syntax, dat)   # ~hours at 20 chains
# saveRDS(blav_fit, "blav_lib_con_fit.rds")
#
# blavaan_draw_names(blav_fit)          # if a lookup ever fails, start here
#

ideo_test <- blavaan::summary(blav_fit, standardize = TRUE, fit.measures = TRUE)
decomp <- blavaan_covariate_decomposition(blav_fit, syntax)
print(decomp, row.names = FALSE, digits = 3)
write.csv(decomp, "blav_lib_con_decomposition.csv", row.names = FALSE)

res <- blavaan_parameterEstimates(blav_fit)

res_nostd <- blavaan_parameterEstimates(blav_fit, standardized = FALSE)

write.csv(res, "res.csv", row.names = F)

# compare_with_dimsem(decomp, lib_con_bayes_res_tab, covariate = "lib_con")
# print(blavaan_hyper_strength(blav_fit, syntax), row.names = FALSE, digits = 3)
#
# blavaan::blavInspect(blav_fit, "psrf")     # convergence
# blavaan::fitMeasures(blav_fit)
#
# WHAT TO LOOK FOR
#   1. Do the TOTALS agree between the two estimators? They are the quantity
#      the likelihood identifies, so they should match within MC error
#      regardless of parameterization. If they do, DimSEM's estimates are
#      right and only its sampling geometry is troubled. If they disagree,
#      the discrepancy is substantive and worth chasing.
#   2. Are the blavaan standardized hyper loadings small and the implied
#      pairwise facet correlations heterogeneous? That is the compound-
#      symmetry misspecification hypothesis, confirmed independently of
#      DimSEM's sampler.
#   3. Is corr(alpha, beta_k) strongly negative in the blavaan draws (as it
#      was in blav_symbolic / blav_party / blav_income)? That confirms the
#      ridge is present under blavaan too and is simply better hidden by the
#      prior anchor.
