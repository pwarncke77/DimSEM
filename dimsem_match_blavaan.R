# =============================================================================
# DimSEM run matched to the blavaan replication setup (32 items, party_scale)
# =============================================================================
#
# Goal: estimate the same hierarchical model as
#
#   bsem(model_bayes, data, fixed.x = FALSE, orthogonal = FALSE,
#        control = list(adapt_delta = 0.9, max_treedepth = 20),
#        bcontrol = list(cores = 16), dp = new_prior,
#        n.chains = 16, burnin = 800, sample = 400)
#
# but under the sum-to-zero identification blavaan cannot represent.
#
# -----------------------------------------------------------------------------
# which DimSEM mode corresponds to your blavaan model
# -----------------------------------------------------------------------------
# I previously told you `ideology =~ a*Fat_k` was a reparameterization of
# DimSEM's hyper_loadings = "unity". That was wrong, and the distinction
# matters here.
#
#   blavaan     hyper loadings fixed at 1, facet disturbances psi_k FREE
#               => cov(eta_j, eta_k) = psi_G  (constant)
#                  var(eta_k)        = psi_G + psi_k  (varies)
#                  correlation       = u_j * u_k,  u_k = sqrt(psi_G/(psi_G+psi_k))
#               i.e. a RANK-1 correlation structure with facet-specific loadings.
#
#   DimSEM "free"   delta_k free, var(G) = 1, var(z_k) = 1
#               => correlation = v_j * v_k,  v_k = delta_k/sqrt(delta_k^2+1)
#               Same family. This is the exact analogue of your blavaan model
#               (the two are reparameterizations of one another).
#
#   DimSEM "unity"  delta_k = 1 AND var(z_k) = 1
#               => correlation = tau^2/(tau^2+1) for EVERY pair: compound
#               symmetry. Strictly nested inside the above.
#
# Your fit says compound symmetry is not innocuous here: the standardized
# hyper loadings are 0.798, 0.554, 0.825, 0.702, 0.559, 0.470, implying facet
# correlations from about 0.22 to 0.66. Unity would force them all equal.
#
# CONSEQUENCE. DimSEM currently implements the sum-to-zero constraint only
# under "unity" (the constraint is exact there; under free loadings the
# symmetry-consistent version is the delta-weighted one, which is not yet
# implemented). So you cannot presently have BOTH the exact structural
# analogue AND the hard constraint in one fit. This script therefore runs two
# models and treats their comparison as the result:
#
#   RUN A  free  + prior        <- structural analogue of blavaan, prior-identified
#   RUN B  unity + sum_to_zero  <- the constrained fit you asked for
#
# Run A should closely reproduce your blavaan numbers (a fidelity check on the
# engine). Run B answers the identification question. Where they agree, the
# conclusion is robust to both choices; where they differ, the difference is
# attributable to compound symmetry, not to the identification convention.
# =============================================================================

## --- 0. Data ------------------------------------------------------------------
# Same frame the blavaan run used: the 32 items plus the 0-1 rescaled covariate.
stopifnot(exists("both_items"), exists("anes2000TS"))
range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
                         (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

dat_ds <- both_items
dat_ds$party_scale <- range_0_1(anes2000TS$party)
stopifnot(nrow(dat_ds) == nrow(anes2000TS))


## --- 1. Proposal with the SAME partition as the blavaan run --------------------
# The walktrap-first pipeline reproduces the published 6-community solution
# exactly, but verify rather than assume: a silent partition change would
# invalidate the whole comparison.

prop32 <- DimSEM_propose(both_items, seed = 123)

blav_measure <- list(   # from your blavaan summary, in its own factor order
  Fat_1 = c("V000680_r","V000676_r","V000679_r"),
  Fat_2 = c("V000620_r","V000614_r","V000550_r","V000587_r","V000681_r","V000693_r"),
  Fat_3 = c("V000685_r","V000683_r","V000682_r","V000731_r","V000776_r","V000677_r"),
  Fat_4 = c("V000687_r","V000674a_r","V000752_r","V000645_r","V000747_r",
            "V000802_r","V000806_r"),
  Fat_5 = c("V000760_r","V000748_r","V000694_r","V001481_r","V000702_r",
            "V000705_r","V000727_r"),
  Fat_6 = c("V000510_r","V000686_r","V000745_r")
)

part <- prop32$selected$partition
ds_measure <- split(names(part), paste0("F", part))

# Match DimSEM's factors to blavaan's by CONTENT (community labels are
# arbitrary in both, so never map by number).
map <- vapply(ds_measure, function(items) {
  hit <- vapply(blav_measure, function(b) setequal(items, b), logical(1))
  if (sum(hit) == 1L) names(blav_measure)[hit] else NA_character_
}, character(1))

cat("\n--- partition check ---\n")
print(data.frame(DimSEM = names(map), blavaan = unname(map),
                 n = lengths(ds_measure[names(map)])), row.names = FALSE)


## --- 2. Impose the same single hyper factor -----------------------------------
# The blavaan model imposes one hyper factor over all six facets by
# construction; it does not test for it. To compare like with like, impose the
# same structure here rather than letting the permutation test decide.

fac_names <- paste0("F", sort(unique(part[!is.na(part)])))
cat("\nhyper block found by the permutation test: ",
    paste(prop32$hyper$permutation$blocks, collapse = " | "), "\n")

prop32$hyper$permutation$blocks  <- list(G1 = fac_names)   # force, to match blavaan
prop32$hyper$permutation$orphans <- character(0)
prop32$hyper$permutation$verdict <- "single"
cat("imposed: G1 = {", paste(fac_names, collapse = ", "), "}\n")


## --- 3. Settings matched to the bsem() call -----------------------------------
#
# Sampler: identical (16 chains, 800 warmup, 400 sampling, adapt_delta 0.9,
#          max_treedepth 20).
# Priors:  matched where the parameter means the same thing in both engines --
#          item intercepts (nu) and item residual SDs (sigma) live on the raw
#          item scale in both, so blavaan's normal(0, 2) and gamma(1, .5)[sd]
#          transfer directly. lambda is normal(0, 2) in both by default.
#          beta is normal(0, 0.5) in both.
#          NOTE: blavaan's gamma(1,.5) priors on the FACET disturbances have no
#          DimSEM counterpart -- those variances are fixed at 1 for scale
#          identification, with the item loadings absorbing the scale instead.
#          Exact prior equivalence is therefore impossible; this is as close as
#          the two parameterizations allow.
# Covariate: NOT z-scored, so party_scale stays on the same 0-1 scale blavaan
#          used and the normal(0, 0.5) prior carries the same meaning.

blav_priors <- list(
  lambda = "normal(0, 2)",
  nu     = "normal(0, 2)",       # blavaan dp: nu = normal(0, 2)
  sigma  = "gamma(1, 0.5)",      # blavaan dp default: gamma(1,.5)[sd]
  tau_g  = "gamma(1, 0.5)",
  delta  = "normal(0, 2)",
  beta   = "normal(0, 0.5)"      # the prior your dp change set
)

bayes_common <- list(
  chains          = 16L,
  parallel_chains = 16L,
  iter_warmup     = 800L,
  iter_sampling   = 400L,
  adapt_delta     = 0.9,
  max_treedepth   = 20L,
  likelihood      = "normal",    # blavaan treats the items as continuous
  priors          = blav_priors
)


## --- RUN A: structural analogue of blavaan (free loadings, prior-identified) ---
#
# fit_free <- DimSEM_estimate(
#   prop32, data = dat_ds,
#   covariates = "party_scale", covariate_targets = "all",
#   engine = "bayes",
#   hyper_loadings = "free",
#   covariate_identification = "prior",
#   standardize_covariates = FALSE,
#   bayes_args = bayes_common,
#   seed = 123
# )
# saveRDS(fit_free, "dimsem_32_party_free_prior.rds")

fit_free <- readRDS("dimsem_32_party_free_prior.rds")
tab_free <- DimSEM_parameterEstimates(fit_free)
write.csv(tab_free, "dimsem_32_party_free_prior.csv", row.names = FALSE)
tab_free

fit_free %>% DimSEM_standardized()

## --- RUN B: the constrained fit (unity loadings, sum-to-zero) ------------------

# fit_stz <- DimSEM_estimate(
#   prop32, data = dat_ds,
#   covariates = "party_scale", covariate_targets = "all",
#   engine = "bayes",
#   hyper_loadings = "unity",
#   covariate_identification = "sum_to_zero",
#   standardize_covariates = FALSE,
#   bayes_args = bayes_common,
#   seed = 123
# )
# saveRDS(fit_stz, "dimsem_32_party_unity_stz.rds")

fit_stz <- readRDS("dimsem_32_party_unity_stz.rds")
tab_stz <- DimSEM_parameterEstimates(fit_stz)
write.csv(tab_stz, "dimsem_32_party_unity_stz.csv", row.names = FALSE)
tab_stz

fit_stz_post <- DimSEM_posterior(fit_stz)

## --- 4. What to look at -------------------------------------------------------

cat("\n=== RUN B: sum-to-zero identified coefficients ===\n")
print(tab_stz[tab_stz$op %in% c("~", ":="), ], row.names = FALSE, digits = 3)

# (a) The constraint holds exactly.
b <- tab_stz$est[tab_stz$op == "~" & grepl("^F[0-9]", tab_stz$lhs)]
cat("\nsum of facet coefficients:", format(sum(b), digits = 3), " (expect ~0)\n")

# (b) alpha should equal the mean of the totals.
tot <- tab_stz$est[tab_stz$op == ":=" & grepl("_total$", tab_stz$lhs)]
a   <- tab_stz$est[tab_stz$op == "~" & grepl("^G", tab_stz$lhs)]
cat("alpha =", round(a, 4), " mean(totals) =", round(mean(tot), 4), "\n")

# (c) Is compound symmetry tenable? Compare Run B's single rho_g against
#     Run A's implied pairwise correlations delta_j*delta_k/sqrt(...).
d <- tab_free$est[tab_free$op == "=~" & grepl("^G", tab_free$lhs)]
v <- d / sqrt(d^2 + 1)
imp <- outer(v, v)[upper.tri(diag(length(v)))]
cat("\nRun A implied facet correlations:", sprintf("%.3f", range(imp)), "(range)\n")
cat("Run B single rho_g             :",
    round(tab_stz$est[tab_stz$lhs == "rho_G1"], 3), "\n")
cat("-> a wide Run A range means unity is imposing a restriction blavaan did not.\n")

# (d) SCALE CAVEAT for cross-engine comparison. blavaan identifies the facet
#     scale with a marker item; DimSEM fixes the disturbance variance to 1.
#     RAW coefficients are therefore NOT comparable between engines. Compare
#     either the standardized quantities or the ratio structure. A DimSEM total
#     on the standardized scale is t_k * sd(x) / sd(eta_k), with
#     var(eta_k) = t_k^2 * var(x) + tau^2 + 1 under unity.
sd_x <- stats::sd(dat_ds$party_scale, na.rm = TRUE)
tau2 <- tab_stz$est[tab_stz$op == "~~" & tab_stz$lhs == "G1"]
std_tot <- tot * sd_x / sqrt(tot^2 * sd_x^2 + tau2 + 1)
cat("\nRun B totals, standardized:", sprintf("%+.3f", std_tot), "\n")
cat("blavaan totals, standardized: +0.354 +0.573 +0.481 +0.441 +0.281 +0.172\n")
cat("(match the ORDER and relative spacing; exact equality is not expected\n")
cat(" because unity constrains the facet correlation structure.)\n")



# =============================================================================
# Standardizing what survived in a reloaded DimSEM fit
# =============================================================================
#
# Both saved fits kept `beta_f_cor`, `beta_g_cor` and `beta_total` (and
# `rho_g` for the unity fit), because those were read by $summary() during
# fitting. That turns out to be exactly enough: every scale factor needed for
# standardization is recoverable from them.
#
# DimSEM's parameterization (see the generated Stan program):
#
#   free loadings :  g_score = X*alpha + z_g        -> var(G) = alpha^2 var(X) + 1
#                    delta_k is a free parameter, tau_g is FIXED at 1
#   unity loadings:  g_score = X*alpha + tau_g*z_g  -> var(G) = alpha^2 var(X) + tau^2
#                    delta_k == 1, tau_g free, rho_g = tau^2/(tau^2+1)
#
#   eta_k = X*beta_k + delta_k*G + z_k,  var(z_k) = 1
#   => var(eta_k) = t_k^2 var(X) + delta_k^2 * tau_g^2 + 1,  t_k = beta_k + delta_k*alpha
#
# The two missing ingredients are recovered rather than refitted:
#
#   delta_k (free fit)  = (beta_total - beta_f_cor) / beta_g_cor, per draw.
#                         This is exact: the generated quantities block defines
#                         beta_total[q,k] = beta_f_cor[q,k] + delta[k]*beta_g_cor[q,g],
#                         so the identity holds draw by draw. It is safe only
#                         while alpha stays away from 0 -- checked below.
#   tau^2   (unity fit) = rho_g / (1 - rho_g), per draw.
#
# Standardized quantities, all computed PER DRAW:
#   beta_k^std  = beta_k  * sd(X) / sd(eta_k)
#   alpha^std   = alpha   * sd(X) / sd(G)
#   t_k^std     = t_k     * sd(X) / sd(eta_k)
#   delta_k^std = delta_k * sd(G) / sd(eta_k)
#
# NOT recoverable without the item-level draws: standardized item loadings,
# intercepts, residual variances. Those need a refit with the patched package.
# =============================================================================

DimSEM_standardize_surviving <- function(object, sd_x, prob = 0.95,
                                         alpha_zero_tol = 0.05) {
  stopifnot(inherits(object, "DimSEM_estimate"))
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior is required.", call. = FALSE)
  }
  unity <- identical(object$param_map$hyper_loadings, "unity")
  facs  <- object$param_map$factors
  m     <- length(facs)

  gm <- function(v) {
    d <- tryCatch(object$fit$draws(variables = v), error = function(e) NULL)
    if (is.null(d)) stop("`", v, "` is not present in this fit.", call. = FALSE)
    posterior::as_draws_matrix(d)
  }

  beta  <- gm("beta_f_cor")                 # draws x K  (single covariate)
  alpha <- as.numeric(gm("beta_g_cor")[, 1])
  tot   <- gm("beta_total")
  stopifnot(ncol(beta) == m, ncol(tot) == m)

  if (unity) {
    rho   <- as.numeric(gm("rho_g")[, 1])
    tau2  <- rho / (1 - rho)
    delta <- matrix(1, length(alpha), m)
  } else {
    # alpha in the denominator: verify it is nowhere near zero before trusting
    # the recovered delta.
    frac <- mean(abs(alpha) < alpha_zero_tol)
    if (frac > 0) {
      warning(sprintf(paste0("%.2f%% of draws have |alpha| < %g; the recovered ",
                             "delta is unstable in those draws. Interpret with ",
                             "care."), 100 * frac, alpha_zero_tol),
              call. = FALSE)
    }
    tau2  <- rep(1, length(alpha))
    delta <- sweep(tot - beta, 1, alpha, "/")
  }

  vx    <- sd_x^2
  v_G   <- alpha^2 * vx + tau2
  v_eta <- tot^2 * vx + delta^2 * tau2 + 1

  b_std <- beta * sd_x / sqrt(v_eta)
  t_std <- tot  * sd_x / sqrt(v_eta)
  d_std <- delta * sqrt(v_G) / sqrt(v_eta)
  a_std <- alpha * sd_x / sqrt(v_G)

  qs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  smry <- function(x, lab) data.frame(
    parameter = lab, est = mean(x), se = stats::sd(x),
    ci.lower = unname(stats::quantile(x, qs[1])),
    ci.upper = unname(stats::quantile(x, qs[2])),
    ess_bulk = posterior::ess_bulk(posterior::as_draws_array(
      array(x, c(length(x), 1L, 1L)))),
    row.names = NULL, stringsAsFactors = FALSE)

  out <- rbind(
    smry(a_std, "G1 ~ cov (std)"),
    do.call(rbind, lapply(seq_len(m), function(k)
      smry(b_std[, k], paste0(facs[k], " ~ cov (std)")))),
    do.call(rbind, lapply(seq_len(m), function(k)
      smry(t_std[, k], paste0(facs[k], " ~ cov (std TOTAL)")))),
    do.call(rbind, lapply(seq_len(m), function(k)
      smry(d_std[, k], paste0("G1 =~ ", facs[k], " (std)"))))
  )
  if (unity) {
    out <- rbind(out, smry(sqrt(rho), "implied facet loading sqrt(rho_g)"))
  }

  # Ridge diagnostics on the RAW draws (scale-free conclusions).
  d_bar <- colMeans(delta)
  s_a <- s_b <- 0.5                       # the priors this fit actually used
  sd_c <- (1 / s_a^2 + sum(d_bar^2) / s_b^2)^(-0.5)
  ridge <- data.frame(
    parameter   = c(facs, "G1"),
    se_raw      = c(apply(beta, 2, stats::sd), stats::sd(alpha)),
    se_total    = c(apply(tot, 2, stats::sd), NA_real_),
    ridge_sd    = c(d_bar * sd_c, sd_c),
    ridge_share = pmin(1, (c(d_bar * sd_c, sd_c) /
                             c(apply(beta, 2, stats::sd), stats::sd(alpha)))^2),
    cor_alpha_beta = c(vapply(seq_len(m), function(k)
      stats::cor(alpha, beta[, k]), numeric(1)), NA_real_),
    row.names = NULL)

  attr(out, "ridge") <- ridge
  attr(out, "delta_posterior_mean") <- d_bar
  attr(out, "hyper_loadings") <- object$param_map$hyper_loadings
  out
}


# =============================================================================
# USAGE
# =============================================================================
# source("DimSEM_standardize_surviving.R")
#
sd_x <- sd(range_0_1(anes2000TS$party), na.rm = TRUE)   # ~0.345
sd_x

std_free <- DimSEM_standardize_surviving(fit_free, sd_x)
std_stz  <- DimSEM_standardize_surviving(fit_stz,  sd_x)

std_free
std_stz

# print(std_free, row.names = FALSE, digits = 3)
# print(attr(std_free, "ridge"), row.names = FALSE, digits = 3)
#
# print(std_stz, row.names = FALSE, digits = 3)
# print(attr(std_stz, "ridge"), row.names = FALSE, digits = 3)
#
# WHAT TO LOOK FOR
#   ridge_share ~ 1  -> that coefficient's interval is the prior's width.
#                       Expect high values in the free/prior fit (my
#                       arithmetic from your summary table put F1 at ~98%,
#                       alpha at ~74%) and ~0 in the sum-to-zero fit.
#   cor_alpha_beta strongly negative -> alpha and beta trade off; the ridge.
#   se_raw vs se_total -> the totals are what the likelihood pins down.
#
# CAVEAT for the unity fit: rho_g sits near its boundary there (0.085, lower
# CI 5.7e-06) with ESS ~26, so tau^2 -- and therefore every standardized
# quantity that divides by sd(eta_k) -- carries that poor mixing. Read the
# standardized unity numbers as indicative, not precise.
