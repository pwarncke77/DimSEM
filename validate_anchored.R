# =============================================================================
# Validating the orientation- and init-aware regime
# =============================================================================
#
# Three fixes have accumulated, and only a factorial design separates their
# contributions:
#
#   anchor_signs   constrains each facet's reference loading to the
#                  orientation implied by the proposal CFA, removing the
#                  spurious mirror optima that trap chains at tau^2 ~ 0
#   cfa_inits      starting values from the proposal CFA (ON by default
#                  whenever the proposal carries one, so switching it OFF is
#                  what isolates it)
#   marginal       integrates the latent disturbances out, removing the
#                  tau <-> latent-cloud funnel
#
# What we already know, from hand-supplied positive inits on toy_simple_cov
# (true rho = 0.30):
#
#   conditional + random init    rho 0.146   ESS 1.3    37.4s
#   conditional + positive init  rho 0.142   ESS 1.3     --
#   marginal    + random init    rho 0.020   ESS 12.2   24.4s
#   marginal    + positive init  rho 0.282   ESS 3093    5.0s
#
# Neither fix worked alone. Section 1 asks whether the ANCHOR reproduces the
# bottom-right cell WITHOUT hand-supplied inits -- i.e. whether it carries
# the orientation by itself, which is what matters at scale, since nobody can
# hand-supply inits for a real battery.
#
# Section 2 then runs ANES against the blavaan reference numbers.
# =============================================================================

devtools::document()
stopifnot(requireNamespace("posterior", quietly = TRUE))


# --- one-row summary of a fit ------------------------------------------------
summarise_fit <- function(fit, label) {
  tb <- DimSEM_parameterEstimates(fit)
  d  <- tb[!is.na(tb$rhat), ]
  rho <- tb[grepl("^rho", tb$lhs), ]
  a   <- tb$est[tb$op == "~"  & grepl("^G", tb$lhs)]
  tt  <- tb$est[tb$op == ":=" & grepl("_total$", tb$lhs)]
  tm  <- tryCatch(fit$fit$time()$chains$total, error = function(e) NA_real_)
  data.frame(
    config      = label,
    rho         = if (nrow(rho)) rho$est[1] else NA_real_,
    rho_lo      = if (nrow(rho)) rho$ci.lower[1] else NA_real_,
    rho_hi      = if (nrow(rho)) rho$ci.upper[1] else NA_real_,
    rho_ess     = if (nrow(rho)) rho$ess_bulk[1] else NA_real_,
    alpha       = if (length(a)) a[1] else NA_real_,
    mean_total  = if (length(tt)) mean(tt) else NA_real_,
    min_ess     = min(d$ess_bulk), max_rhat = max(d$rhat),
    secs        = if (all(is.na(tm))) NA_real_ else max(tm),
    chain_ratio = if (all(is.na(tm))) NA_real_ else max(tm) / min(tm),
    row.names = NULL, stringsAsFactors = FALSE)
}

fit_cfg <- function(prop, data, cov, hl, anchor, cfa_init, marginal,
                    bayes_extra = list(), seed = 123) {
  DimSEM_estimate(
    prop, data = data, covariates = cov, covariate_targets = "all",
    engine = "bayes", hyper_loadings = hl,
    covariate_identification = "sum_to_zero",
    bayes_args = c(list(anchor_signs = anchor, cfa_inits = cfa_init,
                        marginal = marginal), bayes_extra),
    seed = seed)
}


# =============================================================================
# SECTION 1 -- toy_simple_cov, factorial over the three switches
# =============================================================================
# Unity loadings, because that is where the sign trap lives. Ground truth in
# the unity metric: rho = 0.30, alpha = 0.45, totals (0.60, 0.30, 0.45).
# Likert coarsening attenuates, so ~0.28 / ~0.43 is the realistic target.

data(toy_simple_cov)
prop_toy <- DimSEM_propose(toy_simple_cov[, letters[1:12]], seed = 123)
stopifnot(!is.null(prop_toy$cfa), prop_toy$selected$n_dim >= 3L)

toy_args <- list(chains = 2L, parallel_chains = 2L,
                 iter_warmup = 500L, iter_sampling = 2000L,
                 adapt_delta = 0.9, likelihood = "normal")

grid <- expand.grid(anchor = c(FALSE, TRUE), cfa_init = c(FALSE, TRUE),
                    marginal = c(FALSE, TRUE))
res <- list()
for (i in seq_len(nrow(grid))) {
  lab <- sprintf("anchor=%-5s init=%-6s marg=%-5s",
                 grid$anchor[i], if (grid$cfa_init[i]) "cfa" else "random",
                 grid$marginal[i])
  message("\n>>> ", lab)
  f <- try(suppressWarnings(fit_cfg(prop_toy, toy_simple_cov, "cov", "unity",
                                    grid$anchor[i], grid$cfa_init[i],
                                    grid$marginal[i], toy_args)), silent = TRUE)
  if (inherits(f, "try-error")) {
    message("   FAILED: ", conditionMessage(attr(f, "condition"))); next
  }
  res[[lab]] <- summarise_fit(f, lab)
  if (grid$anchor[i] && grid$cfa_init[i] && grid$marginal[i]) fit_best <- f
}
toy_tab <- do.call(rbind, res)

cat("\n============ SECTION 1: toy factorial (truth rho = 0.30) ============\n")
print(toy_tab, row.names = FALSE, digits = 3)
cat("\nref_sign chosen for the toy:",
    paste(fit_best$stan$data$ref_sign, collapse = " "), "\n")
cat("\nHOW TO READ THIS\n")
cat(" * The decisive row is anchor=TRUE init=random marg=TRUE. If its rho is\n")
cat("   near 0.28 with ESS in the thousands, the ANCHOR carries the\n")
cat("   orientation on its own and the CFA inits are a speed-up, not a\n")
cat("   crutch. That is the result that generalizes.\n")
cat(" * Compare anchor=TRUE marg=FALSE against anchor=TRUE marg=TRUE: this\n")
cat("   is now the honest test of whether marginalization is still needed\n")
cat("   once the sign trap is gone.\n")
cat(" * chain_ratio is the max/min chain runtime. A ratio near 1 means the\n")
cat("   chains are starting in comparable places; the ~3x spread seen\n")
cat("   earlier is the signature of some chains fighting a bad mode.\n")


# =============================================================================
# SECTION 2 -- ANES, 32 items, party_scale
# =============================================================================
stopifnot(exists("both_items"), exists("anes2000TS"))
range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
                         (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
dat_ds <- both_items
dat_ds$party_scale <- range_0_1(anes2000TS$party)
dat_x_complete <- dat_ds[!is.na(dat_ds$party_scale), ]

prop32 <- DimSEM_propose(both_items, seed = 123)
fac_names <- paste0("F", sort(unique(prop32$selected$partition[
  !is.na(prop32$selected$partition)])))
cat("\npermutation test found:",
    paste(prop32$hyper$permutation$blocks, collapse = " | "), "\n")
prop32$hyper$permutation$blocks  <- list(G1 = fac_names)
prop32$hyper$permutation$orphans <- character(0)
prop32$hyper$permutation$verdict <- "single"

cat("\nCFA inter-factor correlations (what the orientation is read from):\n")
print(round(prop32$cfa$factor_cor, 3))

# Priors matched to the blavaan reference fit. hyper_scale stays at its
# default "tau": the rho reparameterization was a workaround for the funnel,
# and its beta(2, 8) prior would now pull rho toward 0.2 for no reason.
blav_priors <- list(
  lambda = "normal(0, 2)", nu = "normal(0, 2)",
  sigma  = "gamma(1, 0.5)", tau_g = "gamma(1, 0.5)",
  delta  = "normal(0, 2)", beta = "normal(0, 0.5)")

anes_args <- list(chains = 16L, parallel_chains = 16L,
                  iter_warmup = 800L, iter_sampling = 400L,
                  adapt_delta = 0.9, max_treedepth = 20L,
                  likelihood = "normal", priors = blav_priors)

# fit32 <- fit_cfg(prop32, dat_x_complete, "party_scale", "unity",
#                  anchor = TRUE, cfa_init = TRUE, marginal = TRUE,
#                  bayes_extra = c(anes_args,
#                                  list(standardize_covariates = FALSE)))
# saveRDS(fit32, "dimsem_32_party_unity_anchored.rds")
fit32 <- readRDS("dimsem_32_party_unity_anchored.rds")

tab32 <- DimSEM_parameterEstimates(fit32, standardize = "std.all")


tab32 <- DimSEM_parameterEstimates(fit32)
write.csv(tab32, "dimsem_32_party_unity_anchored.csv", row.names = FALSE)


?DimSEM_parameterEstimates()

library(blavaan)

?standardizedPosterior()


cat("\n============ SECTION 2: ANES structural block ============\n")
print(tab32[tab32$op %in% c("~", ":="), ], row.names = FALSE, digits = 4)
cat("\nref_sign:", paste(fit32$stan$data$ref_sign, collapse = " "),
    " (any -1 means a facet the CFA orients opposite to the rest --\n",
    " under unity that signals a misspecified block, not a convention)\n")

cat("\nESS by block:\n")
for (o in c("=~", "~", "~~", "~1", ":=")) {
  s <- tab32[tab32$op == o & !is.na(tab32$ess_bulk), ]
  if (nrow(s)) cat(sprintf("  %-3s n=%3d  ess %7.0f - %8.0f   max rhat %.4f\n",
                           o, nrow(s), min(s$ess_bulk), max(s$ess_bulk),
                           max(s$rhat)))
}
print(summarise_fit(fit32, "ANES anchored+cfa+marginal"),
      row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------------
# TARGETS from the blavaan reference fit (same model, same metric):
#
#   rho_G1  ~ 0.319  [0.283, 0.355]   NOT the 0.101-0.121 the earlier
#                                      unity runs reported
#   totals  ~ 1.688 / 2.238 / 0.915 / 1.130   (match by size, not by label)
#   alpha   ~ mean(totals) ~ 1.493    (sum-to-zero; blavaan's own 1.194 is
#                                      its prior-anchored split and is NOT
#                                      the comparison target)
#   ESS     in the thousands, rhat ~ 1.00
#
# One legitimate small difference: under anchoring the K reference loadings
# carry a half-normal rather than a normal prior, while blavaan left all 32
# free. That affects 4 of 32 loadings and should be immaterial -- but it is a
# real difference, not an identity.
# ---------------------------------------------------------------------------
rho <- tab32[grepl("^rho", tab32$lhs), ]
tt  <- tab32$est[tab32$op == ":=" & grepl("_total$", tab32$lhs)]
a   <- tab32$est[tab32$op == "~"  & grepl("^G", tab32$lhs)]
cat(sprintf("\nrho_G1 = %.4f [%.4f, %.4f]  ess %.0f  rhat %.4f\n",
            rho$est[1], rho$ci.lower[1], rho$ci.upper[1],
            rho$ess_bulk[1], rho$rhat[1]))
cat("blavaan reference: 0.3186 [0.2827, 0.3547], ess 3576, rhat 1.000\n")
cat(sprintf("totals : %s\nblavaan: 1.688 2.238 0.915 1.130\n",
            paste(sprintf("%.3f", sort(tt, decreasing = TRUE)), collapse = " ")))
cat(sprintf("alpha = %.4f   mean(totals) = %.4f   (equal under sum-to-zero)\n",
            a, mean(tt)))

cat("\nVERDICT\n")
if (!is.na(rho$ess_bulk[1]) && rho$ess_bulk[1] > 400 &&
    abs(rho$est[1] - 0.3186) < 0.06) {
  cat(" rho matches the blavaan reference at high ESS. The unity path is\n")
  cat(" fixed; earlier unity results were biased and should be discarded.\n")
} else if (!is.na(rho$ess_bulk[1]) && rho$ess_bulk[1] > 400) {
  cat(" ESS is healthy but rho disagrees with blavaan. That is a MODEL\n")
  cat(" disagreement, not a sampling one -- stop and investigate before\n")
  cat(" using anything downstream.\n")
} else {
  cat(" ESS still low. Report which block is slow; do not interpret rho.\n")
}


# =============================================================================
# SECTION 3 -- prior sensitivity (run only if Section 2 passed)
# =============================================================================
# With the funnel and the sign trap both gone, the likelihood should dominate,
# so a prior change that previously moved rho substantially ought to be close
# to inert. If rho still tracks the prior, tau^2 is not as well identified as
# the blavaan interval suggests.
#
p2 <- modifyList(blav_priors, list(tau_g = "gamma(1, 2)"))
fit32b <- fit_cfg(prop32, dat_x_complete, "party_scale", "unity",
                  TRUE, TRUE, TRUE,
                  bayes_extra = c(modifyList(anes_args, list(priors = p2)),
                                  list(standardize_covariates = FALSE)))
rbind(summarise_fit(fit32,  "tau_g ~ gamma(1, 0.5)"),
      summarise_fit(fit32b, "tau_g ~ gamma(1, 2)"))



# > # =============================================================================
# > # SECTION 2 -- ANES, 32 items, party_scale
#   > # =============================================================================
# > stopifnot(exists("both_items"), exists("anes2000TS"))
# > range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
#   +                          (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
# > dat_ds <- both_items
# > dat_ds$party_scale <- range_0_1(anes2000TS$party)
# > dat_x_complete <- dat_ds[!is.na(dat_ds$party_scale), ]
# >
#   > prop32 <- DimSEM_propose(both_items, seed = 123)
# Running EGA...
# Auditing item-community memberships (permutation)...
# Fitting CFA (lavaan) on the EGA partition...
# Evaluating hyper-structure: permutation-validated higher-order Louvain...
# Permutation test: 2000 replicates on 31 cores (doSNOW).
# > fac_names <- paste0("F", sort(unique(prop32$selected$partition[
#   +   !is.na(prop32$selected$partition)])))
# > cat("\npermutation test found:",
#       +     paste(prop32$hyper$permutation$blocks, collapse = " | "), "\n")
#
# permutation test found: c("F1", "F2", "F3", "F4")
# > prop32$hyper$permutation$blocks  <- list(G1 = fac_names)
# > prop32$hyper$permutation$orphans <- character(0)
# > prop32$hyper$permutation$verdict <- "single"
# >
#   > cat("\nCFA inter-factor correlations (what the orientation is read from):\n")
#
# CFA inter-factor correlations (what the orientation is read from):
#   > print(round(prop32$cfa$factor_cor, 3))
# F1    F2    F3    F4
# F1 1.000
# F2 0.723 1.000
# F3 0.540 0.228 1.000
# F4 0.273 0.426 0.388 1.000
# >
#   > # Priors matched to the blavaan reference fit. hyper_scale stays at its
#   > # default "tau": the rho reparameterization was a workaround for the funnel,
#   > # and its beta(2, 8) prior would now pull rho toward 0.2 for no reason.
#   > blav_priors <- list(
#     +   lambda = "normal(0, 2)", nu = "normal(0, 2)",
#     +   sigma  = "gamma(1, 0.5)", tau_g = "gamma(1, 0.5)",
#     +   delta  = "normal(0, 2)", beta = "normal(0, 0.5)")
# >
#   > anes_args <- list(chains = 16L, parallel_chains = 16L,
#                       +                   iter_warmup = 800L, iter_sampling = 400L,
#                       +                   adapt_delta = 0.9, max_treedepth = 20L,
#                       +                   likelihood = "normal", priors = blav_priors)
# >
#   > fit32 <- fit_cfg(prop32, dat_x_complete, "party_scale", "unity",
#                      +                  anchor = TRUE, cfa_init = TRUE, marginal = TRUE,
#                      +                  bayes_extra = c(anes_args,
#                                                         +                                  list(standardize_covariates = FALSE)))
# |                                                                                                                                                                       |   0%
# Warning: `data` has 1776 rows but the proposal was derived from 1807 rows; make sure this is the same data set.
# |===================================================================                                                                                                    |  40%
# Standardized covariate(s): party_scale.
# |====================================================================================================                                                                   |  60%
# Sign anchoring: each facet's reference loading is constrained to the orientation implied by the proposal CFA. Under unity hyper loadings a lone facet reflection is not a symmetry, so the mirror configurations are spurious optima that trap chains at a near-zero hyper variance.
# Likelihood: marginalized over the latent disturbances (one Cholesky per missingness pattern). Latent factor scores are not sampled under this program.
# Covariate identification: exact sum-to-zero constraint on child-dimension coefficients within hyper-block(s) G1 (per covariate). Each hyper-level coefficient is thereby identified as the unweighted mean of its children's total covariate effects; each child coefficient as that dimension's deviation from the mean.
# Compiling Stan program...
# Starting values: derived from the proposal CFA.
# Init values were only set for a subset of parameters.
# Missing init values for the following parameters:
#  - chain 1: beta_free, beta_g
#  - chain 2: beta_free, beta_g
#  - chain 3: beta_free, beta_g
#  - chain 4: beta_free, beta_g
#  - chain 5: beta_free, beta_g
#  - chain 6: beta_free, beta_g
#  - chain 7: beta_free, beta_g
#  - chain 8: beta_free, beta_g
#  - chain 9: beta_free, beta_g
#  - chain 10: beta_free, beta_g
#  - chain 11: beta_free, beta_g
#  - chain 12: beta_free, beta_g
#  - chain 13: beta_free, beta_g
#  - chain 14: beta_free, beta_g
#  - chain 15: beta_free, beta_g
#  - chain 16: beta_free, beta_g
#
# To disable this message use options(cmdstanr_warn_inits = FALSE).
# Running MCMC with 16 parallel chains...
#
# > saveRDS(fit32, "dimsem_32_party_unity_anchored.rds")
# > tab32 <- DimSEM_parameterEstimates(fit32)
# > write.csv(tab32, "dimsem_32_party_unity_anchored.csv", row.names = FALSE)
# >
#   > cat("\n============ SECTION 2: ANES structural block ============\n")
#
# ============ SECTION 2: ANES structural block ============
#   > print(tab32[tab32$op %in% c("~", ":="), ], row.names = FALSE, digits = 4)
# lhs op         rhs      est      se ci.lower ci.upper   rhat ess_bulk ess_tail
# F1  ~ party_scale  0.06916 0.02862  0.01489  0.12523 0.9999    10063     5659
# F2  ~ party_scale  0.26706 0.03125  0.20634  0.32826 1.0001     8734     5144
# F3  ~ party_scale -0.20756 0.03374 -0.27361 -0.13917 0.9999     9617     5382
# F4  ~ party_scale -0.12866 0.02881 -0.18408 -0.07359 1.0001    10644     5635
# G1  ~ party_scale  0.51754 0.02428  0.47027  0.56596 1.0000     7617     4770
# F1_total := party_scale  0.58670 0.03553  0.51726  0.65739 1.0000     9156     5298
# F2_total := party_scale  0.78460 0.04053  0.70689  0.86480 1.0003     7352     4798
# F3_total := party_scale  0.30998 0.04331  0.22504  0.39774 0.9999     8387     4854
# F4_total := party_scale  0.38888 0.03665  0.31682  0.45922 0.9999    11910     5134
# rho_G1 :=              0.32119 0.01890  0.28489  0.35841 1.0002     4765     5375
# > cat("\nref_sign:", paste(fit32$stan$data$ref_sign, collapse = " "),
#       +     " (any -1 means a facet the CFA orients opposite to the rest --\n",
#       +     " under unity that signals a misspecified block, not a convention)\n")
#
# ref_sign: 1 1 1 1  (any -1 means a facet the CFA orients opposite to the rest --
#                       under unity that signals a misspecified block, not a convention)
# >
#   > cat("\nESS by block:\n")
#
# ESS by block:
#   > for (o in c("=~", "~", "~~", "~1", ":=")) {
#     +   s <- tab32[tab32$op == o & !is.na(tab32$ess_bulk), ]
#     +   if (nrow(s)) cat(sprintf("  %-3s n=%3d  ess %7.0f - %8.0f   max rhat %.4f\n",
#                                  +                            o, nrow(s), min(s$ess_bulk), max(s$ess_bulk),
#                                  +                            max(s$rhat)))
#     + }
# =~  n= 32  ess    6425 -    14925   max rhat 1.0010
# ~   n=  5  ess    7617 -    10644   max rhat 1.0001
# ~~  n= 33  ess    4765 -    18397   max rhat 1.0016
# ~1  n= 32  ess    6998 -    14847   max rhat 1.0008
# :=  n=  5  ess    4765 -    11910   max rhat 1.0003
# > print(summarise_fit(fit32, "ANES anchored+cfa+marginal"),
#         +       row.names = FALSE, digits = 4)
# config    rho rho_lo rho_hi rho_ess  alpha mean_total min_ess max_rhat  secs chain_ratio
# ANES anchored+cfa+marginal 0.3212 0.2849 0.3584    4765 0.5175     0.5175    4765    1.002 484.9       1.223
# >
#   > # ---------------------------------------------------------------------------
# > # TARGETS from the blavaan reference fit (same model, same metric):
#   > #
#   > #   rho_G1  ~ 0.319  [0.283, 0.355]   NOT the 0.101-0.121 the earlier
#   > #                                      unity runs reported
#   > #   totals  ~ 1.688 / 2.238 / 0.915 / 1.130   (match by size, not by label)
#   > #   alpha   ~ mean(totals) ~ 1.493    (sum-to-zero; blavaan's own 1.194 is
#   > #                                      its prior-anchored split and is NOT
#   > #                                      the comparison target)
#   > #   ESS     in the thousands, rhat ~ 1.00
#   > #
#   > # One legitimate small difference: under anchoring the K reference loadings
#   > # carry a half-normal rather than a normal prior, while blavaan left all 32
#   > # free. That affects 4 of 32 loadings and should be immaterial -- but it is a
#   > # real difference, not an identity.
#   > # ---------------------------------------------------------------------------
# > rho <- tab32[grepl("^rho", tab32$lhs), ]
# > tt  <- tab32$est[tab32$op == ":=" & grepl("_total$", tab32$lhs)]
# > a   <- tab32$est[tab32$op == "~"  & grepl("^G", tab32$lhs)]
# > cat(sprintf("\nrho_G1 = %.4f [%.4f, %.4f]  ess %.0f  rhat %.4f\n",
#               +             rho$est[1], rho$ci.lower[1], rho$ci.upper[1],
#               +             rho$ess_bulk[1], rho$rhat[1]))
#
# rho_G1 = 0.3212 [0.2849, 0.3584]  ess 4765  rhat 1.0002
# > cat("blavaan reference: 0.3186 [0.2827, 0.3547], ess 3576, rhat 1.000\n")
# blavaan reference: 0.3186 [0.2827, 0.3547], ess 3576, rhat 1.000
# > cat(sprintf("totals : %s\nblavaan: 1.688 2.238 0.915 1.130\n",
#               +             paste(sprintf("%.3f", sort(tt, decreasing = TRUE)), collapse = " ")))
# totals : 0.785 0.587 0.389 0.310
# blavaan: 1.688 2.238 0.915 1.130
# > cat(sprintf("alpha = %.4f   mean(totals) = %.4f   (equal under sum-to-zero)\n",
#               +             a, mean(tt)))
# alpha = 0.5175   mean(totals) = 0.5175   (equal under sum-to-zero)
# >
#   > cat("\nVERDICT\n")
#
# VERDICT
# > if (!is.na(rho$ess_bulk[1]) && rho$ess_bulk[1] > 400 &&
#       +     abs(rho$est[1] - 0.3186) < 0.06) {
#   +   cat(" rho matches the blavaan reference at high ESS. The unity path is\n")
#   +   cat(" fixed; earlier unity results were biased and should be discarded.\n")
#   + } else if (!is.na(rho$ess_bulk[1]) && rho$ess_bulk[1] > 400) {
#     +   cat(" ESS is healthy but rho disagrees with blavaan. That is a MODEL\n")
#     +   cat(" disagreement, not a sampling one -- stop and investigate before\n")
#     +   cat(" using anything downstream.\n")
#     + } else {
#       +   cat(" ESS still low. Report which block is slow; do not interpret rho.\n")
#       + }
# rho matches the blavaan reference at high ESS. The unity path is
# fixed; earlier unity results were biased and should be discarded.
# >
#
#
# > # p2 <- modifyList(blav_priors, list(tau_g = "gamma(1, 2)"))
#   > # fit32b <- fit_cfg(prop32, dat_x_complete, "party_scale", "unity",
#   > #                   TRUE, TRUE, TRUE,
#   > #                   bayes_extra = c(modifyList(anes_args, list(priors = p2)),
#   > #                                   list(standardize_covariates = FALSE)))
#   > # rbind(summarise_fit(fit32,  "tau_g ~ gamma(1, 0.5)"),
#   > #       summarise_fit(fit32b, "tau_g ~ gamma(1, 2)"))
#   > p2 <- modifyList(blav_priors, list(tau_g = "gamma(1, 2)"))
# > fit32b <- fit_cfg(prop32, dat_x_complete, "party_scale", "unity",
#                     +                   TRUE, TRUE, TRUE,
#                     +                   bayes_extra = c(modifyList(anes_args, list(priors = p2)),
#                                                         +                                   list(standardize_covariates = FALSE)))
# |                                                                                                                                                                       |   0%
# Warning: `data` has 1776 rows but the proposal was derived from 1807 rows; make sure this is the same data set.
#
# rbind(summarise_fit(fit32,  "tau_g ~ gamma(1, 0.5)"),
#       +       summarise_fit(fit32b, "tau_g ~ gamma(1, 2)"))
# config       rho    rho_lo    rho_hi  rho_ess     alpha mean_total  min_ess max_rhat    secs chain_ratio
# 1 tau_g ~ gamma(1, 0.5) 0.3211942 0.2848909 0.3584092 4765.486 0.5175403  0.5175403 4765.486 1.001576 484.869    1.222783
# 2   tau_g ~ gamma(1, 2) 0.3198990 0.2829681 0.3563903 4827.670 0.5169519  0.5169519 4827.669 1.001663 482.845    1.199317
# > print(summarise_fit(fit32, "ANES anchored+cfa+marginal"),
#         +       row.names = FALSE, digits = 4)
# config    rho rho_lo rho_hi rho_ess  alpha mean_total min_ess max_rhat  secs chain_ratio
# ANES anchored+cfa+marginal 0.3212 0.2849 0.3584    4765 0.5175     0.5175    4765    1.002 484.9       1.223
