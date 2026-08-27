# =============================================================================
# Validating the marginalized program
# =============================================================================
#
# Three sections, cheapest first:
#   1. EQUIVALENCE  -- same toy data, same seed, conditional vs marginal.
#                      They fit the SAME model, so the posteriors must agree
#                      within Monte Carlo error. This is the direct test of
#                      whether the marginalization is implemented correctly,
#                      and it costs minutes.
#   2. PATTERNS     -- the toy is complete data (one missingness pattern), so
#                      section 1 never exercises the pattern machinery. Punch
#                      holes in it and re-run.
#   3. FULL         -- ANES, with blavaan's numbers as the target.
#
# Changes from the previous script, all deliberate:
#   * bayes_args$marginal = TRUE                (otherwise nothing changes)
#   * hyper_scale back to the default "tau"     (the rho workaround is no
#                                                longer needed, and its
#                                                beta(2,8) prior would pull
#                                                rho away from the likelihood)
#   * tau_g = "gamma(1, 0.5)"                   (matches blavaan's default, so
#                                                the comparison is like-for-like;
#                                                gamma(1,2) is the sensitivity
#                                                check, run it afterwards)
# =============================================================================

devtools::document()
devtools::load_all()

check_weighted_constraint <- function(fit, label = "") {
  dr <- posterior::as_draws_matrix(fit$fit$draws(c("beta_f_cor", "delta")))
  nm <- colnames(dr)
  b  <- dr[, grep("^beta_f_cor\\[", nm), drop = FALSE]
  d  <- dr[, grep("^delta\\[", nm), drop = FALSE]
  s  <- rowSums(b * d)
  cat(sprintf("%s max |sum(delta*beta)| = %.3e  (< 1e-7 at default sig_figs)\n",
              label, max(abs(s))))
  invisible(s)
}

# Compare two fits on the quantities that must agree.
compare_fits <- function(fa, fb, la = "conditional", lb = "marginal") {
  ta <- DimSEM_parameterEstimates(fa)
  tb <- DimSEM_parameterEstimates(fb)
  key <- function(x) paste(x$lhs, x$op, x$rhs)
  m <- merge(data.frame(k = key(ta), est_a = ta$est, se_a = ta$se,
                        ess_a = ta$ess_bulk, rhat_a = ta$rhat),
             data.frame(k = key(tb), est_b = tb$est, se_b = tb$se,
                        ess_b = tb$ess_bulk, rhat_b = tb$rhat), by = "k")
  m <- m[!is.na(m$se_a) & m$se_a > 0, ]
  # z-scale the discrepancy by the pooled MC error: |diff| should be small
  # relative to how precisely each fit locates the parameter.
  mcse <- sqrt(m$se_a^2 / pmax(m$ess_a, 1) + m$se_b^2 / pmax(m$ess_b, 1))
  m$z <- (m$est_a - m$est_b) / mcse
  cat(sprintf("\n--- %s vs %s: %d shared parameters ---\n", la, lb, nrow(m)))
  cat(sprintf("  max |est difference|      : %.4f\n",
              max(abs(m$est_a - m$est_b))))
  cat(sprintf("  max |difference| / MCSE   : %.1f   (a few is fine; >10 is not)\n",
              max(abs(m$z))))
  cat(sprintf("  ESS ratio marginal/cond   : median %.1f, min %.1f, max %.1f\n",
              median(m$ess_b / m$ess_a), min(m$ess_b / m$ess_a),
              max(m$ess_b / m$ess_a)))
  print(head(m[order(-abs(m$z)), c("k", "est_a", "est_b", "z")], 8),
        row.names = FALSE, digits = 4)
  invisible(m)
}


# =============================================================================
# SECTION 1 -- EQUIVALENCE on toy_simple_cov
# =============================================================================
data(toy_simple_cov)
items_toy <- letters[1:12]
prop_toy <- DimSEM_propose(toy_simple_cov[, items_toy], seed = 123)
stopifnot(prop_toy$selected$n_dim >= 3L)

toy_args <- list(chains = 2L, parallel_chains = 2L,
                 iter_warmup = 500L, iter_sampling = 2000L,
                 adapt_delta = 0.9, likelihood = "normal")

fit_cond <- DimSEM_estimate(
  prop_toy, data = toy_simple_cov, covariates = "cov",
  covariate_targets = "all", engine = "bayes",
  hyper_loadings = "free", covariate_identification = "sum_to_zero",
  bayes_args = toy_args, seed = 123)

fit_marg <- DimSEM_estimate(
  prop_toy, data = toy_simple_cov, covariates = "cov",
  covariate_targets = "all", engine = "bayes",
  hyper_loadings = "free", covariate_identification = "sum_to_zero",
  bayes_args = c(toy_args, list(marginal = TRUE)), seed = 123)

check_weighted_constraint(fit_marg, "marginal toy:")
cmp <- compare_fits(fit_cond, fit_marg)

# Ground truth (free metric): delta ~ 0.656 for all three, alpha ~ 0.686,
# beta ~ (+0.15, -0.15, 0.00), totals ~ (0.60, 0.30, 0.45). Match facets by
# the size of their total, never by label.
tab_m <- DimSEM_parameterEstimates(fit_marg)
cat("\n--- marginal toy: structural block ---\n")
print(tab_m[tab_m$op %in% c("~", ":=") |
            (tab_m$op == "=~" & grepl("^G", tab_m$lhs)), ],
      row.names = FALSE, digits = 3)
stopifnot(all(tab_m$rhat[!is.na(tab_m$rhat)] < 1.05))
cat("\nSECTION 1 PASSED.\n")


# =============================================================================
# SECTION 2 -- the pattern machinery
# =============================================================================
# toy_simple_cov is complete, so section 1 ran with n_pat = 1 and never
# tested the missingness path -- which is most of the new data builder.
set.seed(99)
toy_miss <- toy_simple_cov
for (j in items_toy) {
  toy_miss[sample(nrow(toy_miss), 40), j] <- NA
}
cat("\ncomplete cases after punching holes:",
    sum(complete.cases(toy_miss[, items_toy])), "of", nrow(toy_miss), "\n")

fit_miss <- DimSEM_estimate(
  prop_toy, data = toy_miss, covariates = "cov",
  covariate_targets = "all", engine = "bayes",
  hyper_loadings = "free", covariate_identification = "sum_to_zero",
  bayes_args = c(toy_args, list(marginal = TRUE)), seed = 123)

cat("missingness patterns used:", fit_miss$stan$data$n_pat, "\n")
cat("respondents retained     :", fit_miss$stan$data$N, "\n")
stopifnot(fit_miss$stan$data$n_pat > 1L)
tab_miss <- DimSEM_parameterEstimates(fit_miss)
stopifnot(all(tab_miss$rhat[!is.na(tab_miss$rhat)] < 1.05))
cat("SECTION 2 PASSED: multiple patterns, all rhat < 1.05.\n")


# =============================================================================
# SECTION 3 -- FULL: 32 items, party_scale, unity loadings
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

# Priors matched to the blavaan reference fit so the comparison is
# like-for-like. gamma(1, 0.5) on the hyper SD is blavaan's own default.
blav_priors <- list(
  lambda = "normal(0, 2)", nu = "normal(0, 2)",
  sigma  = "gamma(1, 0.5)", tau_g = "gamma(1, 0.5)",
  delta  = "normal(0, 2)", beta = "normal(0, 0.5)")

fit_marg32 <- DimSEM_estimate(
  prop32, data = dat_x_complete,
  covariates = "party_scale", covariate_targets = "all",
  engine = "bayes",
  hyper_loadings = "unity",
  covariate_identification = "sum_to_zero",
  standardize_covariates = FALSE,
  bayes_args = list(
    chains = 16L, parallel_chains = 16L,
    iter_warmup = 800L, iter_sampling = 400L,   # short is fine if ESS is high
    adapt_delta = 0.9, max_treedepth = 20L,
    likelihood = "normal", priors = blav_priors,
    marginal = TRUE),
  seed = 123)

saveRDS(fit_marg32, "dimsem_32_party_unity_marginal.rds")
tab32 <- DimSEM_parameterEstimates(fit_marg32)
write.csv(tab32, "dimsem_32_party_unity_marginal.csv", row.names = FALSE)

cat("\n--- FULL (marginal): structural block ---\n")
print(tab32[tab32$op %in% c("~", ":="), ], row.names = FALSE, digits = 4)

cat("\nESS by block:\n")
for (o in c("=~", "~", "~~", "~1", ":=")) {
  s <- tab32[tab32$op == o & !is.na(tab32$ess_bulk), ]
  if (nrow(s)) cat(sprintf("  %-3s n=%3d  ess %7.0f - %8.0f   max rhat %.4f\n",
                           o, nrow(s), min(s$ess_bulk), max(s$ess_bulk),
                           max(s$rhat)))
}

# ---------------------------------------------------------------------------
# TARGETS from the blavaan reference fit (same model, same metric, and a
# marginalized likelihood, which is why it sampled cleanly):
#
#   rho_G1   ~ 0.319  [0.283, 0.355]      NOT the 0.101 the conditional
#                                          program reported
#   totals   ~ 1.688 / 2.238 / 0.915 / 1.130   (match by size, not by label)
#   alpha    ~ mean(totals) ~ 1.493       (sum-to-zero; blavaan's own 1.194 is
#                                          its prior-anchored split, so it is
#                                          NOT the comparison target)
#   ESS      in the thousands, rhat ~ 1.00
#
# If rho lands near 0.32 with high ESS, marginalization is correct and the
# bias is gone. If it lands near 0.10 WITH high ESS, the implementation
# disagrees with blavaan on the model, not just on the sampling -- stop and
# investigate before using anything.
# ---------------------------------------------------------------------------
rho <- tab32[grepl("^rho", tab32$lhs), ]
if (nrow(rho)) {
  cat(sprintf("\nrho_G1 = %.4f [%.4f, %.4f]   ess %.0f   rhat %.4f\n",
              rho$est[1], rho$ci.lower[1], rho$ci.upper[1],
              rho$ess_bulk[1], rho$rhat[1]))
  cat("blavaan reference: 0.3186 [0.2827, 0.3547], ess 3576, rhat 1.000\n")
}
tt <- tab32$est[tab32$op == ":=" & grepl("_total$", tab32$lhs)]
a  <- tab32$est[tab32$op == "~" & grepl("^G", tab32$lhs)]
cat(sprintf("\ntotals: %s\nblavaan: 1.688 2.238 0.915 1.130\n",
            paste(sprintf("%.3f", sort(tt, decreasing = TRUE)), collapse = " ")))
cat(sprintf("alpha = %.4f   mean(totals) = %.4f   (should match under sum-to-zero)\n",
            a, mean(tt)))

# Sensitivity: rerun with tau_g = "gamma(1, 2)" and confirm rho barely moves.
# With the funnel gone the likelihood should dominate, so a prior change that
# previously mattered should now be close to inert.


# saveRDS(fit_marg32, "dimsem_32_party_unity_marginal.rds")
# > tab32 <- DimSEM_parameterEstimates(fit_marg32)
# > write.csv(tab32, "dimsem_32_party_unity_marginal.csv", row.names = FALSE)
# >
#   > cat("\n--- FULL (marginal): structural block ---\n")
#
# --- FULL (marginal): structural block ---
#   > print(tab32[tab32$op %in% c("~", ":="), ], row.names = FALSE, digits = 4)
# lhs op         rhs     est      se   ci.lower ci.upper  rhat ess_bulk ess_tail
# F1  ~ party_scale  0.1907 0.07714  4.196e-02   0.3415 1.000  3021.17   4538.3
# F2  ~ party_scale  0.6889 0.09575  5.092e-01   0.8864 1.005   106.03    508.8
# F3  ~ party_scale -0.5680 0.09290 -7.514e-01  -0.3839 1.005   532.54   1131.8
# F4  ~ party_scale -0.3116 0.08761 -4.880e-01  -0.1464 1.000   115.99   1325.4
# G1  ~ party_scale  1.2991 0.13950  1.120e+00   1.5875 1.012    33.75    271.2
# F1_total := party_scale  1.4897 0.16150  1.246e+00   1.8262 1.012    37.46    272.8
# F2_total := party_scale  1.9880 0.20063  1.696e+00   2.4104 1.011    36.53    272.9
# F3_total := party_scale  0.7311 0.16047  4.663e-01   1.0811 1.003    50.08    247.6
# F4_total := party_scale  0.9875 0.13109  7.755e-01   1.2703 1.014    49.76    283.7
# rho_G1 :=              0.1207 0.13813  3.840e-06   0.3450 1.021    26.21    237.4
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
# =~  n= 32  ess      32 -      144   max rhat 1.0214
# ~   n=  5  ess      34 -     3021   max rhat 1.0116
# ~~  n= 33  ess      26 -    10631   max rhat 1.0215
# ~1  n= 32  ess     126 -     4810   max rhat 1.0030
# :=  n=  5  ess      26 -       50   max rhat 1.0215
# >
#   > # ---------------------------------------------------------------------------
# > # TARGETS from the blavaan reference fit (same model, same metric, and a
#   > # marginalized likelihood, which is why it sampled cleanly):
#   > #
#   > #   rho_G1   ~ 0.319  [0.283, 0.355]      NOT the 0.101 the conditional
#   > #                                          program reported
#   > #   totals   ~ 1.688 / 2.238 / 0.915 / 1.130   (match by size, not by label)
#   > #   alpha    ~ mean(totals) ~ 1.493       (sum-to-zero; blavaan's own 1.194 is
#   > #                                          its prior-anchored split, so it is
#   > #                                          NOT the comparison target)
#   > #   ESS      in the thousands, rhat ~ 1.00
#   > #
#   > # If rho lands near 0.32 with high ESS, marginalization is correct and the
#   > # bias is gone. If it lands near 0.10 WITH high ESS, the implementation
#   > # disagrees with blavaan on the model, not just on the sampling -- stop and
#   > # investigate before using anything.
#   > # ---------------------------------------------------------------------------
# > rho <- tab32[grepl("^rho", tab32$lhs), ]
# > if (nrow(rho)) {
#   +   cat(sprintf("\nrho_G1 = %.4f [%.4f, %.4f]   ess %.0f   rhat %.4f\n",
#                   +               rho$est[1], rho$ci.lower[1], rho$ci.upper[1],
#                   +               rho$ess_bulk[1], rho$rhat[1]))
#   +   cat("blavaan reference: 0.3186 [0.2827, 0.3547], ess 3576, rhat 1.000\n")
#   + }
#
# rho_G1 = 0.1207 [0.0000, 0.3450]   ess 26   rhat 1.0215
# blavaan reference: 0.3186 [0.2827, 0.3547], ess 3576, rhat 1.000
# > tt <- tab32$est[tab32$op == ":=" & grepl("_total$", tab32$lhs)]
# > a  <- tab32$est[tab32$op == "~" & grepl("^G", tab32$lhs)]
# > cat(sprintf("\ntotals: %s\nblavaan: 1.688 2.238 0.915 1.130\n",
#               +             paste(sprintf("%.3f", sort(tt, decreasing = TRUE)), collapse = " ")))
#
# totals: 1.988 1.490 0.987 0.731
# blavaan: 1.688 2.238 0.915 1.130
# > cat(sprintf("alpha = %.4f   mean(totals) = %.4f   (should match under sum-to-zero)\n",
#               +             a, mean(tt)))
# alpha = 1.2991   mean(totals) = 1.2991   (should match under sum-to-zero)
# >
#   > # Sensitivity: rerun with tau_g = "gamma(1, 2)" and confirm rho barely moves.
#   > # With the funnel gone the likelihood should dominate, so a prior change that
#   > # previously mattered should now be close to inert
