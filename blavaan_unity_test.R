# =============================================================================
# Does marginalizing the latents fix the unity model? A blavaan existence proof
# =============================================================================
#
# WHY THIS TEST
#
# DimSEM keeps the hyper-factor disturbance u_i as an explicit parameter (one
# per respondent, 1807 of them) and writes
#
#     eta_ik = X_i * t_k + tau * u_i + z_ik
#
# so tau multiplies a 1807-dimensional latent cloud. Any move in tau requires
# all 1807 u_i to rescale coherently, which a single global step size cannot
# do -- and that coupling is present at every value of tau, not just near
# zero. That is the last remaining explanation for the bulk-ESS ceiling
# (~20-35) after longer chains and the rho reparameterization both failed to
# lift it.
#
# blavaan's target = "stan" integrates the latents out analytically and
# evaluates a marginal multivariate normal over the items. If the SAME model
# mixes cleanly there, the diagnosis is confirmed and marginalizing inside
# DimSEM is worth building. If blavaan struggles too, the diagnosis is wrong
# and several days of implementation are saved.
#
# WHAT MAKES THIS "THE SAME MODEL"
#
# DimSEM's unity mode fixes the hyper loadings at 1 AND the facet disturbance
# variances at 1, letting the item loadings carry the scale. Every pairwise
# facet correlation is then tau^2/(tau^2+1) -- compound symmetry.
#
# The obvious blavaan translation is to keep its marker-variable convention
# and constrain the facet disturbance variances EQUAL via a shared label.
# That was tried and blavaan silently ignored it (npar unchanged at 104 with
# and without the constraint), which is the third constraint form it has
# refused. So this script does not use equality constraints at all. It
# reproduces DimSEM's parameterization directly, with fixed values only:
#
#     Fat_k =~ NA*<first item> + <rest>   free ALL loadings (no marker)
#     Fat_k ~~ 1*Fat_k                    scale set by the disturbance
#     ideology =~ 1*Fat_1 + 1*Fat_2 + ... hyper loadings fixed at 1
#     ideology ~~ ideology                free, equals tau^2
#
# This is not merely equivalent to DimSEM's unity model up to a rescaling --
# it is the same parameterization in the same metric, so the raw parameters
# are directly comparable between the two engines, tau^2 included.
#
# NOTE: no sum-to-zero constraint is needed here. Under unity all hyper
# loadings equal 1, so alpha = mean(totals) exactly, and
# blavaan_parameterEstimates() derives the constrained decomposition
# draw-wise as its `_stz` rows. This test is about MIXING, not identification.
#
# Requires: blavaan, lavaan, posterior, plus blavaan_parameterEstimates.R.
# =============================================================================

library(blavaan)
source("blavaan_parameterEstimates.R")

stopifnot(exists("both_items"), exists("anes2000TS"), exists("prop32"))
range_0_1 <- function(x) (x - min(x, na.rm = TRUE)) /
                         (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))


## --- 1. Same partition, same data as the DimSEM run ---------------------------

part <- prop32$selected$partition
part <- part[!is.na(part)]
measure <- split(names(part), paste0("Fat_", part))
fac <- names(measure)
cat("facets:", length(fac), "| items per facet:", lengths(measure), "\n")
stopifnot(length(fac) >= 3L)

dat_b <- both_items[, unlist(measure), drop = FALSE]
dat_b$party_scale <- range_0_1(anes2000TS$party)
dat_b <- dat_b[!is.na(dat_b$party_scale), ]          # matches dat_x_complete
cat("N used:", nrow(dat_b), "\n")

# Free every loading: NA* on the first item overrides lavaan's auto.fix.first,
# so the facet scale comes from the fixed disturbance variance instead.
fat_measure <- paste(
  vapply(fac, function(f) {
    it <- measure[[f]]
    paste0(f, " =~ NA*", it[1], if (length(it) > 1)
           paste0(" + ", paste(it[-1], collapse = " + ")) else "")
  }, character(1)),
  collapse = " \n ")

# Facet disturbance variances FIXED at 1 (DimSEM: var(z_k) = 1).
fixed_psi <- paste(paste0(fac, " ~~ 1*", fac), collapse = " \n ")

# Hyper loadings FIXED at 1 (DimSEM: delta_k = 1); ideology's own residual
# variance stays free and is tau^2.
hyper_measure <- paste0("ideology =~ ", paste0("1*", fac, collapse = " + "))
hyper_var <- "ideology ~~ ideology"

struct <- paste(paste0(fac, " ~ party_scale"), collapse = " \n ")
hyper_struct <- "ideology ~ party_scale"

model_unity <- paste(fat_measure, " \n ", fixed_psi, " \n ", hyper_measure,
                     " \n ", hyper_var, " \n ", struct, " \n ", hyper_struct)
cat("\n---- model syntax ----\n"); cat(model_unity); cat("\n\n")

new_prior <- dpriors(nu = "normal(0, 2)", alpha = "normal(0, 2)",
                     lambda = "normal(0, 2)", beta = "normal(0, 0.5)")


## --- 2. Verify the restrictions actually took ---------------------------------
# npar is not a reliable witness here (it did not budge for the equality
# version even though the constraint was requested). Check the parameter
# table directly: fixed values are unambiguous there.

pt_check <- lavaan::lavaanify(model_unity, auto.var = TRUE, auto.fix.first = TRUE,
                              auto.cov.lv.x = TRUE, meanstructure = TRUE)
psi_rows <- pt_check[pt_check$op == "~~" & pt_check$lhs %in% fac &
                     pt_check$lhs == pt_check$rhs, ]
lam_rows <- pt_check[pt_check$op == "=~" & pt_check$lhs %in% fac, ]
hyp_rows <- pt_check[pt_check$op == "=~" & pt_check$lhs == "ideology", ]

cat("facet disturbance variances (want free == 0, ustart == 1):\n")
print(psi_rows[, c("lhs", "op", "rhs", "free", "ustart")], row.names = FALSE)
cat("\nhyper loadings (want free == 0, ustart == 1):\n")
print(hyp_rows[, c("lhs", "op", "rhs", "free", "ustart")], row.names = FALSE)
cat("\nitem loadings fixed by auto.fix.first (want NONE):",
    sum(lam_rows$free == 0), "of", nrow(lam_rows), "\n")

stopifnot(
  "facet disturbance variances are not fixed at 1" =
    all(psi_rows$free == 0) && all(psi_rows$ustart == 1),
  "hyper loadings are not fixed at 1" =
    all(hyp_rows$free == 0) && all(hyp_rows$ustart == 1),
  "some item loading is still fixed -- NA* did not free the marker" =
    all(lam_rows$free > 0)
)
cat("\nrestrictions verified: this is DimSEM's unity model, same metric.\n")

## --- 3. The real run ----------------------------------------------------------
# Sampler settings deliberately LIGHTER than the DimSEM run they are being
# compared against (16 chains x 1600 there). If blavaan reaches acceptable ESS
# on fewer draws, the contrast is only sharper.

# blav_unity <- bsem(model = model_unity, data = dat_b,
#                    fixed.x = FALSE, orthogonal = FALSE,
#                    control = list(adapt_delta = 0.9, max_treedepth = 20),
#                    bcontrol = list(cores = 16),
#                    dp = new_prior, mcmcfile = TRUE,
#                    n.chains = 16, burnin = 800, sample = 400, seed = 123)
#
# saveRDS(blav_unity, "blav_32_party_unity.rds")

blav_unity <- readRDS("blav_32_party_unity.rds")

## --- 4. The comparison that decides it ----------------------------------------

tab_b <- blavaan_parameterEstimates(blav_unity, standardized = FALSE)

cat("\n=== blavaan (marginalized) ESS by block ===\n")
for (o in c("=~", "~", "~~", "~1", ":=")) {
  s <- tab_b[tab_b$op == o & !is.na(tab_b$ess_bulk), ]
  if (nrow(s)) cat(sprintf("  %-3s n=%3d  ess_bulk %7.0f - %8.0f   max rhat %.4f\n",
                           o, nrow(s), min(s$ess_bulk), max(s$ess_bulk),
                           max(s$rhat)))
}

# The single decisive number: the hyper-factor variance. In DimSEM this is
# rho_G1 / G1 ~~ G1, stuck at bulk ESS 17-24 with rhat up to 1.033.
psi_g <- tab_b[tab_b$op == "~~" & tab_b$lhs == "ideology", ]
cat("\nhyper-factor variance (blavaan):\n")
print(psi_g, row.names = FALSE, digits = 4)
cat("\nDimSEM equivalent: bulk ESS 23.5, rhat 1.033\n")

cat("\n--- VERDICT ---\n")
if (nrow(psi_g) && !is.na(psi_g$ess_bulk[1])) {
  if (psi_g$ess_bulk[1] > 400) {
    cat("blavaan samples the SAME model cleanly. The latent-coupling\n",
        "diagnosis is confirmed and marginalizing inside DimSEM is\n",
        "justified.\n")
  } else if (psi_g$ess_bulk[1] > 100) {
    cat("Substantially better than DimSEM but short of the >400 guideline.\n",
        "Marginalization helps but may not fully resolve it -- worth\n",
        "implementing, with tempered expectations.\n")
  } else {
    cat("blavaan struggles too. The problem is NOT the parameterization,\n",
        "so marginalization would not fix it. Do not implement; the\n",
        "difficulty is the compound-symmetry restriction itself.\n")
  }
}

# Do the two engines agree on the estimand? In this parameterization the
# implied facet correlation is simply tau^2/(tau^2+1) with tau^2 = the
# ideology residual variance -- the same formula DimSEM reports as rho_g, and
# on the same scale, since the two now share a parameterization.
dr <- blavInspect(blav_unity, "mcmc")
nm <- colnames(as.matrix(dr[[1]]))
i_g <- grep("^ideology~~ideology$", nm)
if (length(i_g) == 1) {
  m <- do.call(rbind, lapply(dr, as.matrix))
  tau2 <- m[, i_g]
  rho_b <- tau2 / (tau2 + 1)
  cat(sprintf("\nimplied facet correlation, blavaan: %.4f [%.4f, %.4f]\n",
              mean(rho_b), quantile(rho_b, .025), quantile(rho_b, .975)))
  cat("DimSEM (rho parameterization)     : 0.1009 [0.0006, 0.3404]\n")
  cat("These should agree closely. If they do, the two engines are fitting\n")
  cat("the same model in the same metric and only the sampling geometry\n")
  cat("differs -- which is exactly the claim being tested.\n")
} else {
  cat("\nCould not locate ideology~~ideology; inspect colnames(as.matrix(dr[[1]])).\n")
}


# > saveRDS(blav_unity, "blav_32_party_unity.rds")
# >
#   >
#   > ## --- 4. The comparison that decides it ----------------------------------------
# >
#   > tab_b <- blavaan_parameterEstimates(blav_unity, standardized = FALSE)
# >
#   > cat("\n=== blavaan (marginalized) ESS by block ===\n")
#
# === blavaan (marginalized) ESS by block ===
#   > for (o in c("=~", "~", "~~", "~1", ":=")) {
#     +   s <- tab_b[tab_b$op == o & !is.na(tab_b$ess_bulk), ]
#     +   if (nrow(s)) cat(sprintf("  %-3s n=%3d  ess_bulk %7.0f - %8.0f   max rhat %.4f\n",
#                                  +                            o, nrow(s), min(s$ess_bulk), max(s$ess_bulk),
#                                  +                            max(s$rhat)))
#     + }
# =~  n= 32  ess_bulk    5859 -    10698   max rhat 1.0007
# ~   n=  5  ess_bulk    2659 -     2904   max rhat 1.0003
# ~~  n= 34  ess_bulk    3576 -    17094   max rhat 1.0010
# ~1  n= 33  ess_bulk    3485 -    16665   max rhat 1.0007
# >
#   > # The single decisive number: the hyper-factor variance. In DimSEM this is
#   > # rho_G1 / G1 ~~ G1, stuck at bulk ESS 17-24 with rhat up to 1.033.
#   > psi_g <- tab_b[tab_b$op == "~~" & tab_b$lhs == "ideology", ]
# > cat("\nhyper-factor variance (blavaan):\n")
#
# hyper-factor variance (blavaan):
#   > print(psi_g, row.names = FALSE, digits = 4)
# lhs op      rhs    est     se ci.lower ci.upper rhat ess_bulk ess_tail
# ideology ~~ ideology 0.4687 0.0399   0.3941   0.5497    1     3576     4786
# > cat("\nDimSEM equivalent: bulk ESS 23.5, rhat 1.033\n")
#
# DimSEM equivalent: bulk ESS 23.5, rhat 1.033
# >
#   > cat("\n--- VERDICT ---\n")
#
# --- VERDICT ---
#   > if (nrow(psi_g) && !is.na(psi_g$ess_bulk[1])) {
#     +   if (psi_g$ess_bulk[1] > 400) {
#       +     cat("blavaan samples the SAME model cleanly. The latent-coupling\n",
#                 +         "diagnosis is confirmed and marginalizing inside DimSEM is\n",
#                 +         "justified.\n")
#       +   } else if (psi_g$ess_bulk[1] > 100) {
#         +     cat("Substantially better than DimSEM but short of the >400 guideline.\n",
#                   +         "Marginalization helps but may not fully resolve it -- worth\n",
#                   +         "implementing, with tempered expectations.\n")
#         +   } else {
#           +     cat("blavaan struggles too. The problem is NOT the parameterization,\n",
#                     +         "so marginalization would not fix it. Do not implement; the\n",
#                     +         "difficulty is the compound-symmetry restriction itself.\n")
#           +   }
#     + }
# blavaan samples the SAME model cleanly. The latent-coupling
# diagnosis is confirmed and marginalizing inside DimSEM is
# justified.
# >
#   > # Do the two engines agree on the estimand? In this parameterization the
#   > # implied facet correlation is simply tau^2/(tau^2+1) with tau^2 = the
#   > # ideology residual variance -- the same formula DimSEM reports as rho_g, and
#   > # on the same scale, since the two now share a parameterization.
#   > dr <- blavInspect(blav_unity, "mcmc")
# > nm <- colnames(as.matrix(dr[[1]]))
# > i_g <- grep("^ideology~~ideology$", nm)
# > if (length(i_g) == 1) {
#   +   m <- do.call(rbind, lapply(dr, as.matrix))
#   +   tau2 <- m[, i_g]
#   +   rho_b <- tau2 / (tau2 + 1)
#   +   cat(sprintf("\nimplied facet correlation, blavaan: %.4f [%.4f, %.4f]\n",
#                   +               mean(rho_b), quantile(rho_b, .025), quantile(rho_b, .975)))
#   +   cat("DimSEM (rho parameterization)     : 0.1009 [0.0006, 0.3404]\n")
#   +   cat("These should agree closely. If they do, the two engines are fitting\n")
#   +   cat("the same model in the same metric and only the sampling geometry\n")
#   +   cat("differs -- which is exactly the claim being tested.\n")
#   + } else {
#     +   cat("\nCould not locate ideology~~ideology; inspect colnames(as.matrix(dr[[1]])).\n")
#     + }
#
# implied facet correlation, blavaan: 0.3186 [0.2827, 0.3547]
# DimSEM (rho parameterization)     : 0.1009 [0.0006, 0.3404]
# These should agree closely. If they do, the two engines are fitting
# the same model in the same metric and only the sampling geometry
# differs -- which is exactly the claim being tested.

tab_b <- blavaan_parameterEstimates(blav_unity, standardized = TRUE)
tab_b[tab_b$op == ":=" & grepl("_total$", tab_b$lhs), ]
write.csv(tab_b, "tab_b.csv", row.names = FALSE)
