# =============================================================================
# blavaan_parameterEstimates() -- a DimSEM_parameterEstimates() analogue
# =============================================================================
#
# Produces the same tidy layout DimSEM_parameterEstimates() returns, so the two
# estimators' output can be compared row-for-row (or rbind()-ed with an
# `engine` column):
#
#   lhs  op  rhs  est  se  ci.lower  ci.upper  rhat  ess_bulk  ess_tail
#
# Conventions mirrored from DimSEM_parameterEstimates():
#   * op vocabulary: "=~" loadings, "~" regressions, "~~" (co)variances,
#     "~1" intercepts, ":=" defined parameters
#   * residual/disturbance terms reported as VARIANCES
#   * parameters fixed by the model show se = 0 and NA diagnostics
#   * est = posterior mean, se = posterior SD, ci = central credible interval
#
# Added `:=` rows (add_defined = TRUE), matching DimSEM's naming so the tables
# line up:
#   <factor>_total   := <covariate>   the likelihood-identified total effect
#                                     t_k = beta_k + delta_k * alpha
#   <factor>_stz     := <covariate>   the same totals re-split under DimSEM's
#   <hyper>_stz      := <covariate>   sum-to-zero (loading-weighted) convention
#
# Standardized vs. unstandardized. blavaan identifies latent scale with a
# marker item; DimSEM fixes the disturbance variance to 1. Raw loadings and
# raw structural coefficients are therefore NOT comparable across the two.
# standardized = TRUE (the default) puts both on a common footing and is what
# the earlier blav_*.csv exports used.
#
# Requires: blavaan, lavaan, posterior.
# =============================================================================


# --- name matching (see blavaan_comparison.R for the rationale) --------------
# Exact comparison, never regex: "=~", "~~" and "~" are full of metacharacters
# and a hand-rolled character class mis-parses under R's default TRE engine.
.blav_match_par <- function(nm, target) {
  i <- which(nm == target)
  if (length(i) == 1L) return(i)
  no_ws <- function(x) gsub("[[:space:]]+", "", x)
  i <- which(no_ws(nm) == no_ws(target))
  if (length(i) == 1L) return(i)
  key <- function(x) gsub("[^[:alnum:]_]", "", x)
  i <- which(key(nm) == key(target))
  if (length(i) == 1L) return(i)
  NA_integer_                     # caller decides: fixed, aliased, or missing
}

.blav_par_name <- function(lhs, op, rhs) {
  # lavaan writes intercepts as lhs = "V1", op = "~1", rhs = ""
  if (identical(op, "~1")) paste0(lhs, "~1") else paste0(lhs, op, rhs)
}


# --- draws as an [iteration, chain, parameter] array -------------------------
#
# blavInspect(fit, "mcmc") preserves the chain structure; standardizedPosterior
# collapses it, stacking chains in order. We reshape the latter using the chain
# count and per-chain length taken from the former, and verify the total draw
# count matches before doing so -- diagnostics computed on a mis-shaped array
# would be silently wrong rather than obviously wrong.
.blav_draws_array <- function(fit, standardized = TRUE) {
  if (!requireNamespace("blavaan", quietly = TRUE)) stop("blavaan required.")
  mc <- blavaan::blavInspect(fit, "mcmc")
  n_chain <- length(mc)
  n_iter <- nrow(as.matrix(mc[[1]]))

  if (!standardized) {
    arr <- array(NA_real_, c(n_iter, n_chain, ncol(as.matrix(mc[[1]]))),
                 dimnames = list(NULL, NULL, colnames(as.matrix(mc[[1]]))))
    for (ch in seq_len(n_chain)) arr[, ch, ] <- as.matrix(mc[[ch]])
    return(arr)
  }

  sp <- as.matrix(blavaan::standardizedPosterior(fit))
  if (nrow(sp) != n_iter * n_chain) {
    stop("standardizedPosterior returned ", nrow(sp), " draws but the fit has ",
         n_chain, " chains x ", n_iter, " iterations (= ", n_chain * n_iter,
         "). Cannot recover the chain layout; use standardized = FALSE.",
         call. = FALSE)
  }
  arr <- array(NA_real_, c(n_iter, n_chain, ncol(sp)),
               dimnames = list(NULL, NULL, colnames(sp)))
  for (ch in seq_len(n_chain)) {
    arr[, ch, ] <- sp[((ch - 1L) * n_iter + 1L):(ch * n_iter), , drop = FALSE]
  }
  arr
}


# Always return an [iteration, chain] matrix, including when there is only
# one chain (where `arr[, , j]` would drop to a vector).
.blav_slice <- function(arr, j) {
  matrix(arr[, , j], nrow = dim(arr)[1], ncol = dim(arr)[2])
}


# --- summarize one parameter's draws ----------------------------------------
.blav_summarise <- function(x, prob = 0.95) {
  # x: [iteration, chain] matrix for a single parameter
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior is required for rhat / ESS.")
  }
  q <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  d <- posterior::as_draws_array(array(x, c(nrow(x), ncol(x), 1L)))
  v <- as.vector(x)
  data.frame(
    est = mean(v),
    se = stats::sd(v),
    ci.lower = unname(stats::quantile(v, q[1])),
    ci.upper = unname(stats::quantile(v, q[2])),
    rhat = posterior::rhat(d),
    ess_bulk = posterior::ess_bulk(d),
    ess_tail = posterior::ess_tail(d),
    stringsAsFactors = FALSE
  )
}


# --- main --------------------------------------------------------------------
blavaan_parameterEstimates <- function(fit, prob = 0.95, standardized = TRUE,
                                       add_defined = TRUE, hyper_name = NULL,
                                       covariates = NULL, verbose = TRUE) {
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan required.")
  arr <- .blav_draws_array(fit, standardized = standardized)
  nm <- dimnames(arr)[[3]]

  pt <- as.data.frame(lavaan::parTable(fit), stringsAsFactors = FALSE)
  pt <- pt[pt$op %in% c("=~", "~", "~~", "~1"), , drop = FALSE]

  rows <- vector("list", nrow(pt))
  missed <- character(0)
  for (i in seq_len(nrow(pt))) {
    target <- .blav_par_name(pt$lhs[i], pt$op[i], pt$rhs[i])
    j <- .blav_match_par(nm, target)

    if (is.na(j) && !is.null(pt$label) && nzchar(pt$label[i])) {
      # Equality-constrained parameters (e.g. ideology =~ a*Fat_k) may be
      # stored once under the shared label rather than per row.
      j <- .blav_match_par(nm, pt$label[i])
    }
    if (is.na(j) && !is.null(pt$free) && pt$free[i] > 0L) {
      # Same free index as an earlier row -> reuse that row's column.
      sib <- which(pt$free == pt$free[i])
      sib <- sib[sib != i]
      for (k in sib) {
        j <- .blav_match_par(nm, .blav_par_name(pt$lhs[k], pt$op[k], pt$rhs[k]))
        if (!is.na(j)) break
      }
    }

    if (!is.na(j)) {
      s <- .blav_summarise(.blav_slice(arr, j), prob = prob)
    } else if (!is.null(pt$free) && pt$free[i] == 0L) {
      # Fixed by the model: mirror DimSEM (se = 0, no diagnostics).
      val <- if (!is.null(pt$est)) pt$est[i] else pt$ustart[i]
      s <- data.frame(est = val, se = 0, ci.lower = val, ci.upper = val,
                      rhat = NA_real_, ess_bulk = NA_real_,
                      ess_tail = NA_real_)
    } else {
      missed <- c(missed, target)
      s <- data.frame(est = NA_real_, se = NA_real_, ci.lower = NA_real_,
                      ci.upper = NA_real_, rhat = NA_real_,
                      ess_bulk = NA_real_, ess_tail = NA_real_)
    }
    rows[[i]] <- cbind(
      data.frame(lhs = pt$lhs[i], op = pt$op[i], rhs = pt$rhs[i],
                 stringsAsFactors = FALSE), s)
  }
  out <- do.call(rbind, rows)

  if (length(missed) > 0 && isTRUE(verbose)) {
    warning(length(missed), " free parameter(s) had no matching draw column ",
            "and are reported as NA: ",
            paste(utils::head(missed, 5), collapse = ", "),
            if (length(missed) > 5) ", ..." else "",
            "\nInspect the available names with ",
            "colnames(blavaan::standardizedPosterior(fit)).", call. = FALSE)
  }

  if (isTRUE(add_defined)) {
    dfn <- .blav_defined_rows(arr, pt, prob = prob, hyper_name = hyper_name,
                              covariates = covariates)
    if (!is.null(dfn)) out <- rbind(out, dfn)
  }

  attr(out, "standardized") <- standardized
  attr(out, "n_chains") <- dim(arr)[2]
  attr(out, "n_draws") <- dim(arr)[1] * dim(arr)[2]
  rownames(out) <- NULL
  out
}


# --- defined parameters: total effects and the sum-to-zero re-split ----------
#
# The hyper factor is detected structurally: a latent that loads on other
# latents (i.e. appears as lhs of "=~" with rhs entries that are themselves
# lhs of "=~"). Covariates are the observed variables on the rhs of "~".
.blav_defined_rows <- function(arr, pt, prob = 0.95, hyper_name = NULL,
                               covariates = NULL) {
  nm <- dimnames(arr)[[3]]
  lat <- unique(pt$lhs[pt$op == "=~"])
  if (is.null(hyper_name)) {
    is_hyper <- vapply(lat, function(L) {
      kids <- pt$rhs[pt$op == "=~" & pt$lhs == L]
      length(kids) > 0 && all(kids %in% lat)
    }, logical(1))
    hyper_name <- lat[is_hyper]
  }
  if (length(hyper_name) != 1L) return(NULL)     # no unique hyper factor

  kids <- pt$rhs[pt$op == "=~" & pt$lhs == hyper_name]
  if (length(kids) < 2L) return(NULL)

  if (is.null(covariates)) {
    covariates <- unique(pt$rhs[pt$op == "~" & pt$lhs == hyper_name])
  }
  covariates <- covariates[nzchar(covariates)]
  if (length(covariates) == 0L) return(NULL)

  get <- function(target) {
    j <- .blav_match_par(nm, target)
    if (is.na(j)) NULL else .blav_slice(arr, j)
  }
  delta <- lapply(kids, function(k) get(paste0(hyper_name, "=~", k)))
  if (any(vapply(delta, is.null, logical(1)))) return(NULL)

  res <- list()
  for (cv in covariates) {
    alpha <- get(paste0(hyper_name, "~", cv))
    beta <- lapply(kids, function(k) get(paste0(k, "~", cv)))
    if (is.null(alpha) || any(vapply(beta, is.null, logical(1)))) next

    totals <- lapply(seq_along(kids), function(k) beta[[k]] + delta[[k]] * alpha)
    # Loading-weighted projection onto the sum-to-zero subspace: alpha becomes
    # sum(delta * t) / sum(delta^2), i.e. DimSEM's convention. Under equality-
    # constrained (unity) loadings this is the plain mean of the totals.
    num <- Reduce(`+`, lapply(seq_along(kids), function(k)
      delta[[k]] * totals[[k]]))
    den <- Reduce(`+`, lapply(delta, function(d) d^2))
    a_stz <- num / den
    b_stz <- lapply(seq_along(kids), function(k)
      totals[[k]] - delta[[k]] * a_stz)

    for (k in seq_along(kids)) {
      res[[length(res) + 1L]] <- cbind(
        data.frame(lhs = paste0(kids[k], "_total"), op = ":=", rhs = cv,
                   stringsAsFactors = FALSE),
        .blav_summarise(totals[[k]], prob))
    }
    res[[length(res) + 1L]] <- cbind(
      data.frame(lhs = paste0(hyper_name, "_stz"), op = ":=", rhs = cv,
                 stringsAsFactors = FALSE),
      .blav_summarise(a_stz, prob))
    for (k in seq_along(kids)) {
      res[[length(res) + 1L]] <- cbind(
        data.frame(lhs = paste0(kids[k], "_stz"), op = ":=", rhs = cv,
                   stringsAsFactors = FALSE),
        .blav_summarise(b_stz[[k]], prob))
    }
  }
  if (length(res) == 0L) NULL else do.call(rbind, res)
}


# --- convenience: stack a blavaan and a DimSEM table for comparison ----------
compare_parameter_tables <- function(blav_tab, dimsem_tab,
                                     ops = c("~", ":=")) {
  b <- blav_tab[blav_tab$op %in% ops, , drop = FALSE]
  d <- dimsem_tab[dimsem_tab$op %in% ops, , drop = FALSE]
  b$engine <- "blavaan"
  d$engine <- "DimSEM"
  keep <- c("engine", "lhs", "op", "rhs", "est", "se", "ci.lower", "ci.upper",
            "rhat", "ess_bulk")
  d <- d[, intersect(keep, names(d)), drop = FALSE]
  b <- b[, intersect(keep, names(b)), drop = FALSE]
  out <- rbind(b, d)
  out[order(out$lhs, out$rhs, out$engine), , drop = FALSE]
}


# =============================================================================
# USAGE
# =============================================================================
# source("blavaan_parameterEstimates.R")
# blav_fit <- readRDS("blav_lib_con_fit.rds")
#
# tab <- blavaan_parameterEstimates(blav_fit)          # standardized
# print(head(tab, 20), row.names = FALSE, digits = 3)
# write.csv(tab, "blav_lib_con_parameterEstimates.csv", row.names = FALSE)
#
# # unstandardized (blavaan's own marker-item metric)
# tab_raw <- blavaan_parameterEstimates(blav_fit, standardized = FALSE)
#
# # side by side with the DimSEM table
# compare_parameter_tables(tab, lib_con_bayes_res_tab)
#
# # worst-mixing rows, same triage as for the DimSEM output
# head(tab[order(tab$ess_bulk), ], 10)
#
# NOTE ON COMPARABILITY
#   The `~` rows are on blavaan's own (prior-anchored) parameterization, so
#   they are NOT directly comparable to a DimSEM run fitted with
#   covariate_identification = "sum_to_zero". Compare the `<factor>_total`
#   rows -- identified by the likelihood under either convention -- or the
#   `_stz` rows, which apply DimSEM's convention to the blavaan posterior.
