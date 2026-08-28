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


fit_cond_u <- DimSEM_estimate(prop_toy, data = toy_simple_cov, covariates = "cov",
                              covariate_targets = "all", engine = "bayes",
                              hyper_loadings = "unity", covariate_identification = "sum_to_zero",
                              bayes_args = toy_args, seed = 123)

saveRDS(fit_cond_u, "fit_cond_u.rds")

fit_cond_u_tab <- DimSEM_parameterEstimates(fit_cond_u)
write.csv(fit_cond_u_tab, "fit_cond_u_tab.csv", row.names = FALSE)


fit_marg_u <- DimSEM_estimate(prop_toy, data = toy_simple_cov, covariates = "cov",
                              covariate_targets = "all", engine = "bayes",
                              hyper_loadings = "unity", covariate_identification = "sum_to_zero",
                              bayes_args = c(toy_args, list(marginal = TRUE, hyper_scale = "rho")), seed = 123)

saveRDS(fit_marg_u, "fit_marg_u.rds")

fit_marg_u_tab <- DimSEM_parameterEstimates(fit_marg_u)
write.csv(fit_marg_u_tab, "fit_marg_u_tab.csv", row.names = FALSE)


## refit the marginal arm with the SAME priors as the conditional arm
fit_marg_u <- DimSEM_estimate(prop_toy, data = toy_simple_cov, covariates = "cov",
                              covariate_targets = "all", engine = "bayes",
                              hyper_loadings = "unity", covariate_identification = "sum_to_zero",
                              bayes_args = c(toy_args, list(marginal = TRUE)), seed = 123)

## 1. do the two programs agree? (this is the correctness test)
cmp <- compare_fits(fit_cond_u, fit_marg_u, "conditional", "marginal")

## 2. constraint check that works under BOTH loading modes
con_check <- function(fit, lab = "") {
  dr <- posterior::as_draws_matrix(fit$fit$draws("beta_f_cor"))
  b  <- dr[, grep("^beta_f_cor\\[", colnames(dr)), drop = FALSE]
  d  <- tryCatch({
    dd <- posterior::as_draws_matrix(fit$fit$draws("delta"))
    dd[, grep("^delta\\[", colnames(dd)), drop = FALSE]
  }, error = function(e) matrix(1, nrow(b), ncol(b)))   # unity: delta == 1
  cat(sprintf("%s max |sum(delta*beta)| = %.3e\n", lab, max(abs(rowSums(b * d)))))
}
con_check(fit_cond_u, "conditional:")
con_check(fit_marg_u, "marginal:   ")

## 3. did the pattern machinery engage, and what did it cost?
cat("n_pat:", fit_marg_u$stan$data$n_pat,
    " N:", fit_marg_u$stan$data$N, "\n")
cat("elapsed  conditional:", round(fit_cond_u$timing$elapsed, 1),
    " marginal:", round(fit_marg_u$timing$elapsed, 1), "\n")

## 4. structural blocks side by side, against known truth
print(fit_cond_u_tab[fit_cond_u_tab$op %in% c("~", ":=") |
                       grepl("^rho", fit_cond_u_tab$lhs), ], row.names = FALSE, digits = 4)
print(DimSEM_parameterEstimates(fit_marg_u)[
  {t <- DimSEM_parameterEstimates(fit_marg_u);
  t$op %in% c("~", ":=") | grepl("^rho", t$lhs)}, ],
  row.names = FALSE, digits = 4)



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

cat(fit_marg_u$stan$model_code)

## 2. what the data actually say, to adjudicate 0.097 vs 0.012
Y <- toy_simple_cov[, letters[1:12]]
R <- cor(Y, use = "pairwise")
cat("within-facet mean r:", round(mean(R[1:4,1:4][upper.tri(diag(4))]), 4), "\n")
cat("cross-facet  mean r:", round(mean(R[1:4,5:8]), 4), "\n")



# Done in 24.4 s.
# >
#   > ## 1. do the two programs agree? (this is the correctness test)
#   > cmp <- compare_fits(fit_cond_u, fit_marg_u, "conditional", "marginal")
#
# --- conditional vs marginal: 45 shared parameters ---
#   max |est difference|      : 0.1976
# max |difference| / MCSE   : 4.6   (a few is fine; >10 is not)
# ESS ratio marginal/cond   : median 17.6, min 1.8, max 4152.1
# k    est_a    est_b      z
# i ~~ i  0.43239  0.43889 -4.578
# k ~~ k  0.57494  0.56975  4.232
# c ~~ c  0.52641  0.53219 -4.087
# e ~~ e  0.49735  0.50769 -3.194
# f ~~ f  0.49591  0.48515  2.674
# b ~~ b  0.46761  0.47049 -2.516
# F3 ~ cov -0.01522 -0.01158 -2.417
# b ~1   3.01847  3.02093 -1.616
# >
#   > ## 2. constraint check that works under BOTH loading modes
#   > con_check <- function(fit, lab = "") {
#     +   dr <- posterior::as_draws_matrix(fit$fit$draws("beta_f_cor"))
#     +   b  <- dr[, grep("^beta_f_cor\\[", colnames(dr)), drop = FALSE]
#     +   d  <- tryCatch({
#       +     dd <- posterior::as_draws_matrix(fit$fit$draws("delta"))
#       +     dd[, grep("^delta\\[", colnames(dd)), drop = FALSE]
#       +   }, error = function(e) matrix(1, nrow(b), ncol(b)))   # unity: delta == 1
#     +   cat(sprintf("%s max |sum(delta*beta)| = %.3e\n", lab, max(abs(rowSums(b * d)))))
#     + }
# > con_check(fit_cond_u, "conditional:")
# conditional: max |sum(delta*beta)| = 1.000e-08
# > con_check(fit_marg_u, "marginal:   ")
# marginal:    max |sum(delta*beta)| = 1.000e-08
# >
#   > ## 3. did the pattern machinery engage, and what did it cost?
#   > cat("n_pat:", fit_marg_u$stan$data$n_pat,
#         +     " N:", fit_marg_u$stan$data$N, "\n")
# n_pat: 1  N: 400
# > cat("elapsed  conditional:", round(fit_cond_u$timing$elapsed, 1),
#       +     " marginal:", round(fit_marg_u$timing$elapsed, 1), "\n")
# elapsed  conditional: 37.4  marginal: 24.4
# >
#   > ## 4. structural blocks side by side, against known truth
#   > print(fit_cond_u_tab[fit_cond_u_tab$op %in% c("~", ":=") |
#                            +                        grepl("^rho", fit_cond_u_tab$lhs), ], row.names = FALSE, digits = 4)
# lhs op rhs      est      se   ci.lower ci.upper  rhat ess_bulk ess_tail
# F1  ~ cov  0.18201 0.05260  8.112e-02  0.28816 1.088    8.955   105.00
# F2  ~ cov -0.16679 0.05160 -2.676e-01 -0.06677 1.081   10.974   134.70
# F3  ~ cov -0.01522 0.05037 -1.133e-01  0.08580 1.002 1408.718  2309.78
# G1  ~ cov  0.40055 0.05343  3.090e-01  0.51228 1.455    1.993    16.72
# F1_total := cov  0.58256 0.08279  4.353e-01  0.75769 1.406    2.123    16.40
# F2_total := cov  0.23376 0.06613  1.102e-01  0.36778 1.077   11.748   134.46
# F3_total := cov  0.38532 0.07282  2.499e-01  0.53372 1.209    3.706    19.69
# rho_G1 :=      0.14553 0.14116  2.104e-05  0.35192 2.123    1.317    18.46
# > print(DimSEM_parameterEstimates(fit_marg_u)[
#   +   {t <- DimSEM_parameterEstimates(fit_marg_u);
#   +   t$op %in% c("~", ":=") | grepl("^rho", t$lhs)}, ],
#   +   row.names = FALSE, digits = 4)
# lhs op rhs       est       se   ci.lower ci.upper   rhat ess_bulk ess_tail
# F1  ~ cov  0.164204 0.048908  7.055e-02  0.26260 1.0001   5671.7     3071
# F2  ~ cov -0.152625 0.048343 -2.480e-01 -0.05845 1.0000   5905.0     3046
# F3  ~ cov -0.011579 0.050176 -1.109e-01  0.08515 1.0008   5366.4     2890
# G1  ~ cov  0.367681 0.034508  2.993e-01  0.43637 1.0001   5952.7     3180
# F1_total := cov  0.531886 0.062129  4.140e-01  0.65666 1.0000   5214.0     2852
# F2_total := cov  0.215056 0.057197  1.039e-01  0.32392 0.9999   6818.0     3078
# F3_total := cov  0.356102 0.060725  2.371e-01  0.47761 1.0004   5698.7     3017
# rho_G1 :=      0.006313 0.008535  9.070e-06  0.03038 1.0147    483.1     2068
#


## A. Are the chains in different sign modes? (30 seconds, no refit)
m  <- posterior::as_draws_matrix(fit_marg_u$fit$draws("lambda_un"))
ch <- posterior::as_draws_df(fit_marg_u$fit$draws("lambda_un"))$.chain
fac <- rep(1:3, each = 4)
for (k in 1:3)
  cat("facet", k, "mean loading sign by chain:",
      round(tapply(rowMeans(m[, fac == k, drop = FALSE]), ch, function(x) mean(sign(x))), 2), "\n")

## B. Start every loading positive and see if rho moves to ~0.35
init_pos <- function() list(lambda_un = rep(0.6, 12), tau_g = array(0.6, 1))
fit_fix <- DimSEM_estimate(prop_toy, data = toy_simple_cov, covariates = "cov",
                           covariate_targets = "all", engine = "bayes",
                           hyper_loadings = "unity", covariate_identification = "sum_to_zero",
                           bayes_args = c(toy_args, list(marginal = TRUE, init = init_pos)), seed = 123)
DimSEM_parameterEstimates(fit_fix)[c("lhs","op","rhs","est","rhat","ess_bulk")] |>
  subset(grepl("^rho|^G1", lhs))




# > ## A. Are the chains in different sign modes? (30 seconds, no refit)
#   > m  <- posterior::as_draws_matrix(fit_marg_u$fit$draws("lambda_un"))
# > ch <- posterior::as_draws_df(fit_marg_u$fit$draws("lambda_un"))$.chain
# > fac <- rep(1:3, each = 4)
# > for (k in 1:3)
#   +   cat("facet", k, "mean loading sign by chain:",
#           +       round(tapply(rowMeans(m[, fac == k, drop = FALSE]), ch, function(x) mean(sign(x))), 2), "\n")
# facet 1 mean loading sign by chain: 1 1
# facet 2 mean loading sign by chain: -1 1
# facet 3 mean loading sign by chain: 1 -1
# >
#   > ## B. Start every loading positive and see if rho moves to ~0.35
#   > init_pos <- function() list(lambda_un = rep(0.6, 12), tau_g = array(0.6, 1))
# > fit_fix <- DimSEM_estimate(prop_toy, data = toy_simple_cov, covariates = "cov",
#                              +                            covariate_targets = "all", engine = "bayes",
#                              +                            hyper_loadings = "unity", covariate_identification = "sum_to_zero",
#                              +                            bayes_args = c(toy_args, list(marginal = TRUE, init = init_pos)), seed = 123)
#
# Both chains finished successfully.
# Mean chain execution time: 4.4 seconds.
# Total execution time: 5.0 seconds.
#
# |=======================================================================================================================================================================| 100%
# Done in 5.7 s.
# > DimSEM_parameterEstimates(fit_fix)[c("lhs","op","rhs","est","rhat","ess_bulk")] |>
#   +   subset(grepl("^rho|^G1", lhs))
# lhs op rhs       est      rhat ess_bulk
# 13     G1 =~  F1 1.0000000        NA       NA
# 14     G1 =~  F2 1.0000000        NA       NA
# 15     G1 =~  F3 1.0000000        NA       NA
# 43     G1 ~~  G1 0.3975517 0.9998905 3093.204
# 47     G1  ~ cov 0.4335106 1.0014184 5107.892
# 51 rho_G1 :=     0.2820611 0.9998905 3093.204

## 1. Is marginalization even necessary, or was the sign trap the whole story?


