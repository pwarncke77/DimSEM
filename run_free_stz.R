# =============================================================================
# delta-weighted sum-to-zero under FREE loadings: trial run, then the full run
# =============================================================================
#
# SECTION 1 is a cheap smoke test on toy_simple_cov, where the true
#   decomposition is known analytically. Run it first: it compiles the new
#   Stan program, exercises the projection, and checks the answer against
#   ground truth in a few minutes rather than a few hours.
#
# SECTION 2 is the real thing: the 32-item battery with party_scale, matched
#   to the settings of the earlier free+prior (A) and unity+sum_to_zero (B)
#   runs so the three are directly comparable.
#
# SECTION 3 verifies the constraint draw-wise and lines the three runs up.
#
# No `force_recompile` needed anywhere: the generated program changed (new
# prj_g data plus the projection block), and cmdstanr keys its cache on the
# code, so it will rebuild on its own.
# =============================================================================

# --- shared check: does the constraint actually hold, per draw? ---------------
# The reported quantities are sign-corrected, and the constraint survives the
# correction (a per-factor reflection flips delta_k and beta_k together, so
# each product is unchanged). Checking on the corrected draws is therefore
# both valid and the more useful test, since those are what gets reported.
check_weighted_constraint <- function(fit, label = "") {
  dr <- posterior::as_draws_matrix(fit$fit$draws(c("beta_f_cor", "delta")))
  nm <- colnames(dr)
  b  <- dr[, grep("^beta_f_cor\\[", nm), drop = FALSE]
  d  <- dr[, grep("^delta\\[", nm), drop = FALSE]
  s  <- rowSums(b * d)
  cat(sprintf("%s  max |sum(delta*beta)| across draws: %.3e   (expect ~1e-12)\n",
              label, max(abs(s))))
  invisible(s)
}


# =============================================================================
# SECTION 1 -- TRIAL: toy_simple_cov, known ground truth
# =============================================================================
#
# toy_simple_cov was generated as
#     eta_k = t_k*cov + g + z_k,   var(g) = 0.43, var(z_k) = 1, cov ~ N(0,1)
# with total effects t = (0.60, 0.30, 0.45). Re-expressed in DimSEM's free
# parameterization (var(G) disturbance = 1, var(z_k) = 1) that is
#
#     delta_k = sqrt(0.43) = 0.656   (identical across facets by construction)
#     alpha   = sum(delta*t)/sum(delta^2) = 0.686
#     beta    = t - delta*alpha            = (+0.15, -0.15, 0.00)
#
# Note alpha is 0.686 here, not the 0.45 stored in attr(., "truth"): that
# figure is in the UNITY metric, where the hyper factor carries variance
# tau^2 rather than 1. The totals and the betas are the same in both.
#
# The items are 5-point Likert, so expect mild attenuation under the normal
# likelihood. Judge the trial on ordering, on the constraint holding, and on
# coverage -- not on exact equality.

data(toy_simple_cov)
items_toy <- letters[1:12]

prop_toy <- DimSEM_propose(toy_simple_cov[, items_toy], seed = 123)
cat("\ntrial partition:", prop_toy$selected$n_dim, "dimension(s);",
    "hyper:", paste(names(prop_toy$hyper$permutation$blocks), collapse = ", "),
    "\n")
if (prop_toy$selected$n_dim < 3L) {
  stop("Fewer than three facets proposed; free hyper loadings need at least ",
       "three children. Inspect prop_toy before continuing.")
}

fit_toy <- DimSEM_estimate(
  prop_toy, data = toy_simple_cov,
  covariates = "cov", covariate_targets = "all",
  engine = "bayes",
  hyper_loadings = "free",
  covariate_identification = "sum_to_zero",
  seed = 123,
  bayes_args = list(chains = 2L, parallel_chains = 2L,
                    iter_warmup = 500L, iter_sampling = 2000L,
                    adapt_delta = 0.9, likelihood = "normal"))

tab_toy <- DimSEM_parameterEstimates(fit_toy)
cat("\n--- TRIAL: structural block ---\n")
print(tab_toy[tab_toy$op %in% c("~", ":=") |
              (tab_toy$op == "=~" & grepl("^G", tab_toy$lhs)), ],
      row.names = FALSE, digits = 3)

check_weighted_constraint(fit_toy, "TRIAL:")

write.csv(tab_toy, "tab_toy.csv")

cat("\nexpected under the free parameterization:\n")
cat("  delta_k  ~ 0.656 (all three similar)\n")
cat("  alpha    ~ 0.686\n")
cat("  beta     ~ (+0.15, -0.15, 0.00), delta-weighted sum 0\n")
cat("  totals   ~ (0.60, 0.30, 0.45)   [facet labels are arbitrary --\n")
cat("             match by which total is largest/smallest, not by number]\n")
cat("  diagnostics: rhat ~ 1.00 throughout, NO rho_g row (free loadings\n")
cat("             have no hyper variance parameter, which is the point)\n")

stopifnot(all(tab_toy$rhat[!is.na(tab_toy$rhat)] < 1.05))
cat("\nTRIAL PASSED: all rhat < 1.05.\n")


# =============================================================================
# SECTION 2 -- FULL: 32 items, party_scale
# =============================================================================
stopifnot(exists("both_items"), exists("anes2000TS"))
range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
                         (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

dat_ds <- both_items
dat_ds$party_scale <- range_0_1(anes2000TS$party)
stopifnot(nrow(dat_ds) == nrow(anes2000TS))

prop32 <- DimSEM_propose(both_items, seed = 123)
cat("\nfull-run partition:", prop32$selected$n_dim, "dimension(s)\n")
print(table(paste0("F", prop32$selected$partition)))

# Impose one hyper factor over all facets, as runs A and B did, so the three
# are comparable. (Reported first, so a change in the permutation test's own
# verdict is visible rather than silently overwritten.)
fac_names <- paste0("F", sort(unique(prop32$selected$partition[
  !is.na(prop32$selected$partition)])))
cat("permutation test found:",
    paste(prop32$hyper$permutation$blocks, collapse = " | "), "\n")
prop32$hyper$permutation$blocks  <- list(G1 = fac_names)
prop32$hyper$permutation$orphans <- character(0)
prop32$hyper$permutation$verdict <- "single"
cat("imposed: G1 = {", paste(fac_names, collapse = ", "), "}\n")

if (length(fac_names) < 3L) {
  stop("Free hyper loadings need at least three children; got ",
       length(fac_names), ".")
}

# Settings identical to runs A and B.
blav_priors <- list(
  lambda = "normal(0, 2)", nu = "normal(0, 2)",
  sigma  = "gamma(1, 0.5)", tau_g = "gamma(1, 2)",
  delta  = "normal(0, 2)", beta = "normal(0, 0.5)"
)

dat_x_complete <- dat_ds[!is.na(dat_ds$party_scale),]

# fit_wstz <- DimSEM_estimate(
#   prop32, data = dat_x_complete,
#   covariates = "party_scale", covariate_targets = "all",
#   engine = "bayes",
#   hyper_loadings = "unity",
#   covariate_identification = "sum_to_zero",
#   standardize_covariates = FALSE,
#   bayes_args = list(
#     chains = 16L, parallel_chains = 16L,
#     iter_warmup = 800L, iter_sampling = 1600L,
#     adapt_delta = 0.95, max_treedepth = 20L,
#     likelihood = "normal", priors = blav_priors
#   ),seed = 123)
#
# saveRDS(fit_wstz, "dimsem_32_party_unity_stz.rds")
fit_wstz <- readRDS("dimsem_32_party_unity_stz.rds")
tab_wstz <- DimSEM_parameterEstimates(fit_wstz)

sd_x <- sd(range_0_1(anes2000TS$party), na.rm = TRUE)
tab_stnd <- DimSEM_standardize_surviving(fit_wstz, sd_x)
# tab_stnd <- DimSEM_standardized(fit_wstz) ## Error in object$covariate_scaling$standardized[[cv]] : subscript out of bounds

write.csv(tab_wstz, "dimsem_32_party_unity_stz.csv", row.names = FALSE)

# =============================================================================
# SECTION 3 -- verification and three-way comparison
# =============================================================================

cat("\n--- FULL: structural block ---\n")
print(tab_wstz[tab_wstz$op %in% c("~", ":=") |
               (tab_wstz$op == "=~" & grepl("^G", tab_wstz$lhs)), ],
      row.names = FALSE, digits = 3)

check_weighted_constraint(fit_wstz, "FULL: ")

# (a) alpha should equal the delta-weighted projection of the totals.
d  <- tab_wstz$est[tab_wstz$op == "=~" & grepl("^G", tab_wstz$lhs)]
tt <- tab_wstz$est[tab_wstz$op == ":=" & grepl("_total$", tab_wstz$lhs)]
a  <- tab_wstz$est[tab_wstz$op == "~"  & grepl("^G", tab_wstz$lhs)]
cat(sprintf("\nalpha = %.4f   sum(delta*t)/sum(delta^2) = %.4f\n",
            a, sum(d * tt) / sum(d^2)))
cat("(small gaps are the usual mean-of-a-ratio vs ratio-of-means slack)\n")

# (b) Mixing: the whole point of moving the constraint to free loadings.
cat("\nESS summary (bulk):\n")
for (o in c("=~", "~", ":=")) {
  s <- tab_wstz[tab_wstz$op == o & !is.na(tab_wstz$ess_bulk), ]
  if (nrow(s)) cat(sprintf("  %-3s n=%3d  ess %6.0f - %6.0f   max rhat %.3f\n",
                           o, nrow(s), min(s$ess_bulk), max(s$ess_bulk),
                           max(s$rhat)))
}
cat("no rho_g row should exist here -- free loadings have no hyper variance\n")
cat("parameter, so the boundary that throttled run B is simply absent.\n")

# (c) Totals against runs A and B. These are the likelihood-identified
#     quantities and must agree across all three identifications.
if (file.exists("dimsem_32_party_free_prior.rds") &&
    file.exists("dimsem_32_party_unity_stz.rds"))
  tA <- DimSEM_parameterEstimates(readRDS("dimsem_32_party_free_prior.rds"))
  tB <- DimSEM_parameterEstimates(readRDS("dimsem_32_party_unity_stz.rds"))
  gt <- function(x) x$est[x$op == ":=" & grepl("_total$", x$lhs)]
  cmp <- data.frame(facet = fac_names,
                    A_free_prior = gt(tA), B_unity_stz = gt(tB),
                    C_free_stz = tt)
  cat("\n--- totals across the three runs (raw metric) ---\n")
  print(cmp, row.names = FALSE, digits = 4)
  cat("\nNOTE: raw totals are NOT directly comparable between free and unity\n")
  cat("runs (the facet scale differs). Compare the STANDARDIZED totals:\n")
  sd_x <- sd(dat_ds$party_scale, na.rm = TRUE)
  sA <- DimSEM_standardize_surviving(readRDS("dimsem_32_party_free_prior.rds"), sd_x)
  sB <- DimSEM_standardize_surviving(readRDS("dimsem_32_party_unity_stz.rds"), sd_x)
  sC <- DimSEM_standardize_surviving(fit_wstz, sd_x)
  pick <- function(z) z$est[grepl("std TOTAL", z$parameter)]
  print(data.frame(facet = fac_names, A = pick(sA), B = pick(sB),
                   C = pick(sC)), row.names = FALSE, digits = 4)
  cat("\nThese three columns should agree closely. Where they do, the\n")
  cat("substantive result is invariant to the identification convention;\n")
  cat("the alpha/beta SPLIT is where the three runs differ by design.\n")


# (d) Ridge report: should now show near-zero ridge shares.
if (exists("DimSEM_ridgeReport")) {
  cat("\n--- ridge report (expect shares ~0 and no strong alpha/beta trade) ---\n")
  print(DimSEM_ridgeReport(fit_wstz), row.names = FALSE, digits = 3)
}


# =============================================================================
# WHAT TO EXPECT IN THE FULL RUN
# =============================================================================
# delta       ~ (2.83, 0.87, 0.47, 0.20), as recovered from run A
# alpha       ~ 1.58  -- run A's PROJECTION value, not its prior-anchored 1.42
# betas       delta-weighted sum exactly 0
# totals      match runs A and B once standardized
# mixing      like run A (ess in the thousands), not like run B (ess ~30);
#             there is no rho_g and no tau^2 boundary to get stuck against
# se(beta)/se(total) ~ 1, as in run B -- the ridge is gone
#
# IF MIXING IS STILL POOR, the first suspects are (i) delta_1 being large and
# weakly determined, which also gave run A's F1_total an ess of only 248, and
# (ii) 400 sampling iterations being thin for 16 chains. Raising
# iter_sampling to 1000 is the cheap first response before anything
# structural.

  # =============================================================================
  # BONUS Section -- FULL: 32 items, income_scale
  # =============================================================================
  stopifnot(exists("both_items"), exists("anes2000TS"))
  range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

  dat_ds <- both_items
  dat_ds$income_scale <- range_0_1(anes2000TS$income)
  stopifnot(nrow(dat_ds) == nrow(anes2000TS))

  prop32 <- DimSEM_propose(both_items, seed = 123)
  cat("\nfull-run partition:", prop32$selected$n_dim, "dimension(s)\n")
  print(table(paste0("F", prop32$selected$partition)))

  # Impose one hyper factor over all facets, as runs A and B did, so the three
  # are comparable. (Reported first, so a change in the permutation test's own
  # verdict is visible rather than silently overwritten.)
  fac_names <- paste0("F", sort(unique(prop32$selected$partition[
    !is.na(prop32$selected$partition)])))
  cat("permutation test found:",
      paste(prop32$hyper$permutation$blocks, collapse = " | "), "\n")
  prop32$hyper$permutation$blocks  <- list(G1 = fac_names)
  prop32$hyper$permutation$orphans <- character(0)
  prop32$hyper$permutation$verdict <- "single"
  cat("imposed: G1 = {", paste(fac_names, collapse = ", "), "}\n")

  if (length(fac_names) < 3L) {
    stop("Free hyper loadings need at least three children; got ",
         length(fac_names), ".")
  }

  # Settings identical to runs A and B.
  blav_priors <- list(
    lambda = "normal(0, 2)", nu = "normal(0, 2)",
    sigma  = "gamma(1, 0.5)", tau_g = "gamma(1, 0.5)",
    delta  = "normal(0, 2)", beta = "normal(0, 0.5)"
  )

  # fit_wstz <- DimSEM_estimate(
  #   prop32, data = dat_ds,
  #   covariates = "income_scale", covariate_targets = "all",
  #   engine = "bayes",
  #   hyper_loadings = "free",
  #   covariate_identification = "sum_to_zero",
  #   standardize_covariates = FALSE,
  #   bayes_args = list(
  #     chains = 16L, parallel_chains = 16L,
  #     iter_warmup = 800L, iter_sampling = 400L,
  #     adapt_delta = 0.9, max_treedepth = 20L,
  #     likelihood = "normal", priors = blav_priors),
  #   seed = 123
  # )
  #
  # saveRDS(fit_wstz, "dimsem_32_income_free_stz.rds")
  fit_wstz <- readRDS("dimsem_32_income_free_stz.rds")
  tab_wstz <- DimSEM_parameterEstimates(fit_wstz)

  sd_x <- sd(range_0_1(anes2000TS$income), na.rm = TRUE)
  tab_stnd <- DimSEM_standardize_surviving(fit_wstz, sd_x)
  # tab_stnd <- DimSEM_standardized(fit_wstz) ## Error in object$covariate_scaling$standardized[[cv]] : subscript out of bounds

  write.csv(tab_wstz, "dimsem_32_income_free_stz.csv", row.names = FALSE)

