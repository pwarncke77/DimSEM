#' Estimate the DimSem-proposed structural equation model
#'
#' `DimSEM_estimate()` takes a proposal object produced by
#' [DimSEM_propose()] and estimates the corresponding (hierarchical)
#' structural equation model: first-order latent dimensions defined by the
#' proposed item partition, optional second-order hyper-factors (unit
#' loadings, following the family-resemblance convention), and external
#' covariates regressed onto the latent variables.
#'
#' @param proposal A `"DimSEM_proposal"` object from [DimSEM_propose()].
#' @param data The data set the proposal was derived from (the proposal
#'   object does not store raw data). Must contain all analyzed item
#'   columns and any requested covariates. Row count is checked against
#'   the proposal.
#' @param covariates Optional character vector of covariate column names in
#'   `data` to regress latent variables on.
#' @param covariate_targets Which latent variables receive the covariate
#'   regressions. `"all"` (default) regresses every first-order dimension
#'   and every hyper-factor on all covariates simultaneously. Alternatively,
#'   a named list mapping latent names to character vectors of covariates, e.g.
#'   `list(G1 = c("age", "educ"), F2 = "age")`. When a hyper-factor and all
#'   of its child dimensions are targeted at once (as under `"all"`), the
#'   Bayesian engine with `hyper_loadings = "unity"` identifies the
#'   decomposition through an exact sum-to-zero constraint on the child
#'   coefficients within each such block; see Details.
#' @param model_syntax Optional complete lavaan-style model syntax. When
#'   supplied, it overides the syntax generated from the proposal entirely
#'   (measurement, hyper, and structural parts); `covariates` and
#'   `include_hyper` are then ignored for syntax construction (covariate
#'   columns are still passed to the estimator).
#' @param hyper_loadings How hyper-factor loadings are treated. `"unity"` (default)
#'   fixes all hyper loadings to 1 with a freely estimated hyper-factor variance,
#'   following Warncke & Azevedo (n.d.): every sub-dimension resembles the family
#'   prototype equally, implying a uniform (compound-symmetric) latent correlation
#'   within each hyper-block.
#'   `"free"`estimates each first-order dimension's loading on its
#'   hyper-factor freely, with the hyper-factor variance fixed to 1 for
#'   identification; this permits differential family resemblance across
#'   sub-dimensions.  Free loadings require at least three sub-dimensions per
#'   hyper-factor for identification; two-child blocks automatically fall
#'   back to unity loadings with a warning.
#' @param covariate_identification How the covariate decomposition is
#'   identified when a hyper-factor and all of its children are targeted
#'   simultaneously (Bayesian engine only). `"sum_to_zero"` (default)
#'   imposes the exact within-block constraint described in Details --
#'   unweighted under `hyper_loadings = "unity"`, loading-weighted under
#'   `"free"`;
#'   `"prior"` imposes no constraint and lets the coefficient priors anchor
#'   the split, reproducing the conventional blavaan-style
#'   parameterization. Under `"prior"` the child coefficients do not sum to
#'   zero and the hyper-level coefficient is a prior-shrunken mean of the
#'   children's total effects, so the estimand depends on the prior scale
#'   ratio and on the number of detected facets; a prior sensitivity check
#'   should accompany results reported this way. The constraint is the
#'   limiting case of `"prior"` as the hyper-coefficient prior becomes
#'   diffuse relative to the child priors.
#' @param standardize_covariates Logical; if `TRUE` (default), continuous
#'   and ordered-categorical covariates with more than two distinct values
#'   are z-standardized before estimation, while binary covariates are kept
#'   in their raw metric (so their coefficients remain interpretable as
#'   group contrasts). The classification and scaling applied to each
#'   covariate is recorded in the `covariate_scaling` element of the
#'   returned object.
#' @param include_hyper Logical; include the proposed hyper-factor level
#'   (default `TRUE`). Hyper blocks are taken from the permutation strategy
#'   when the permutation test validates them. Hyper loadings are
#'   fixed to 1 with a freely estimated hyper-factor variance, following
#'   the family-resemblance convention; with standardized first-order
#'   disturbances this implies the compound-symmetry structure discussed in
#'   the package documentation.
#' @param engine `"ml"` (lavaan; default) or `"bayes"` (custom Stan
#'   backend, see Details).
#' @param estimator Character; any maximum-likelihood estimator supported
#'   by lavaan (`"ML"`, `"MLR"` (default), `"MLM"`, `"MLMV"`, `"MLMVS"`,
#'   `"MLF"`). Non-ML estimators (e.g. `"WLSMV"`) are accepted and passed
#'   through, but note that FIML missing handling requires the ML family.
#' @param missing Character; any missing-data method supported by lavaan:
#'   `"ml"` (FIML; default), `"ml.x"`, `"two.stage"`,
#'   `"robust.two.stage"`, `"pairwise"`, `"listwise"`, or
#'   `"doubly.robust"`. Passed to lavaan unmodified.
#' @param lavaan_args Named list of further arguments for [lavaan::sem()]
#'   (merged over `list(std.lv = TRUE)`).
#' @param bayes_args Named list controlling the Bayesian backend. Key
#'   entries: `backend` (`"cmdstanr"`, default, or `"dry_run"` to return
#'   the generated Stan program and data without compiling),
#'   `likelihood` (`"normal"`, default, or `"ordinal"` for an
#'   ordered-logistic item model), `chains` (4), `parallel_chains`
#'   (defaults to `chains`), `threads_per_chain` (within-chain
#'   parallelism via `reduce_sum`; default 1), `opencl_ids` (e.g.
#'   `c(0, 0)` to enable GPU acceleration via OpenCL), `iter_warmup`
#'   (1000), `iter_sampling` (1000), `adapt_delta` (.9), `max_treedepth`,
#'   `hyper_scale`, `hyper_param`, `marginal` (see Details),
#'   `init`, `step_size`, `refresh`,
#'   `rho_floor` (.05; see Details), and `priors` (a named list overriding
#'   any prior in the generated program). Regression coefficients default
#'   to `normal(0, 0.5)` for standardized covariates, following Warncke &
#'   Azevedo (n.d.), and `normal(0, 1)` for raw binary predictors (one raw
#'   unit is about two predictor SDs for a balanced binary variable, so the
#'   doubled scale preserves the same per-SD regularization).
#' @param force_recompile Logical; if `TRUE`, any cached Stan executable
#'   for the generated program is deleted before compilation, forcing a
#'   full rebuild (Bayesian engine only; ignored by the ML engine). This is
#'   the reset button for stale-binary situations that cmdstanr's cache
#'   cannot detect on its own: the cache key is the generated Stan code, so
#'   a change that leaves the code identical but alters the required
#'   compile flags -- most notably toggling `threads_per_chain` on the CPU
#'   variant (which requires `STAN_THREADS`) -- would otherwise silently
#'   reuse a binary built with the old flags. Also useful after CmdStan
#'   upgrades or interrupted/corrupted builds. Where the installed
#'   cmdstanr version supports it, the flag is additionally passed to
#'   [cmdstanr::cmdstan_model()] as `force_recompile`.
#' @param progress Logical; display stage progress and (for the ML engine
#'   on larger data) a pilot-based estimate of time to completion. The
#'   Bayesian backend additionally streams cmdstan's native per-chain
#'   iteration progress.
#' @param seed Optional integer seed (passed to the estimator).
#' @param verbose Logical; print progress messages.
#'
#' @details
#' The Bayesian engine translates the assembled model into a Stan program
#' via DimSem's own syntax-generation helpers rather than calling blavaan.
#' The generated program uses a non-centered parameterization, models only
#' the observed item entries (full-information under MAR; no case
#' deletion), supports within-chain parallelism through `reduce_sum`
#' (enable with `threads_per_chain > 1`), and can be compiled with OpenCL
#' support for GPU acceleration (`opencl_ids`). Portions of the model
#' structure follow conventions established in blavaan (Merkle & Rosseel,
#' 2018); see the package citation file.
#'
#' \strong{Identification of simultaneous covariate effects.} When a
#' covariate is regressed on a hyper-factor and on all of its child
#' dimensions at once, the child-level mean structure informs only the
#' total effects (child coefficient plus loading-weighted hyper path);
#' shifting a constant out of every child coefficient and into the
#' hyper-level coefficient is exactly likelihood-invariant, so the
#' decomposition is a flat ridge that only the priors would resolve --
#' with prior-dependent estimands, and (where the residual shared
#' variance is also small) severely degraded MCMC mixing as the
#' practical symptoms. Under `hyper_loadings = "unity"` and the default
#' `covariate_identification = "sum_to_zero"`, the Bayesian
#' engine therefore imposes an exact sum-to-zero constraint on the child
#' coefficients within every fully targeted hyper-block (implemented
#' through a marginal-variance-preserving orthonormal contrast basis).
#' The hyper-level coefficient is then identified as the unweighted mean
#' of the children's total covariate effects -- extending the unweighted
#' family-resemblance convention of the unity loading structure to the
#' structural layer -- and each child coefficient as that dimension's
#' deviation from the mean. Child coefficients within a constrained block
#' are consequently contrasts: they sum to zero by construction, and a
#' negative value indicates a below-average, not necessarily an opposing,
#' covariate effect. Under `hyper_loadings = "free"` no analogous fixed
#' constraint is loading-weighted, \eqn{\sum_k \delta_k \beta_k = 0}. Since
#' the loadings are sampled, it cannot be pre-built as a fixed basis and is
#' applied instead as a runtime projection inside the Stan program, giving
#' \eqn{\alpha = \sum_k \delta_k t_k / \sum_k \delta_k^2} -- the
#' loading-weighted component of the children's total effects -- with each
#' child coefficient the residual from it.
#'
#' Which loading mode to pair the constraint with is a modelling decision,
#' not a formality. `"unity"` fixes both the hyper loadings and the facet
#' disturbance variances, which forces every pairwise facet correlation to
#' the same value; when the facets relate to the hyper factor at genuinely
#' different strengths, the single hyper variance is driven toward its lower
#' boundary and the boundary-adjacent geometry degrades mixing for that
#' variance and for everything downstream of it (the `rho_g` diagnostic
#' below flags this). `"free"` expresses the same heterogeneity through
#' unconstrained loadings, which have no boundary, so the constraint is
#' available without the compound-symmetry restriction.
#' Setting `covariate_identification = "prior"` deliberately selects that
#' soft, prior-anchored identification for the unity case too. It leaves
#' the coefficients unconstrained and reproduces the conventional
#' blavaan-style parameterization, in which the ridge settles at the
#' penalized projection
#' \eqn{\alpha^* = \sum_k \delta_k t_k / (s_\beta^2/s_\alpha^2 +
#' \sum_k \delta_k^2)}: the hyper-level coefficient is a shrunken
#' (loading-weighted) mean of the children's total effects \eqn{t_k},
#' and the child coefficients retain a common residue of the shared
#' effect rather than summing to zero. This is a legitimate Bayesian
#' identification strategy provided its two dependencies are reported --
#' the split is a function of the prior scale ratio, and the shrinkage
#' factor drifts toward one as the number of detected facets grows, so
#' the hyper-level estimand is not directly comparable across item pools
#' of different size. Under either strategy the total effects are what
#' the likelihood identifies. Under
#' the ML engine the same redundancy manifests as inflated standard
#' errors, as documented. With covariates present, the reported
#' within-block correlation `rho_g` is conditional on the covariates, not
#' the marginal correlation.
#'
#' \strong{Marginalized likelihood.} `bayes_args$marginal = TRUE` switches
#' to a program that integrates the latent disturbances out analytically
#' instead of sampling them. In the default program the hyper-factor scale
#' multiplies an \eqn{N}-dimensional cloud of latent disturbances, so every
#' move in that scale requires all \eqn{N} of them to rescale coherently --
#' something a single global step size cannot do. On a 1807-respondent ANES
#' battery this capped effective sample size for the hyper variance at
#' roughly 20-35 and, as a `blavaan` cross-check established, biased it
#' downward, because the chains never left the low-variance region. Longer
#' chains, a bounded reparameterization and centering all failed to lift it.
#'
#' Since everything in the normal-likelihood model is Gaussian, the latents
#' can be integrated out exactly. Writing \eqn{t_k} for the total covariate
#' effect on facet \eqn{k}, the marginal moments per respondent are
#' \eqn{E[y_j] = \nu_j + \lambda_j X t_{k(j)}} and a covariance that is
#' diagonal plus rank \eqn{K + G}: one column per facet carrying its own
#' loadings, one per hyper-block carrying the loadings weighted by the
#' parent. Evaluated with the Woodbury identity and the matrix determinant
#' lemma, at \eqn{O(P_{obs} R^2)} per missingness pattern rather than
#' \eqn{O(P_{obs}^3)}, with one Cholesky factorization per pattern shared by
#' its respondents. No latent cloud survives and the hyper scale becomes an
#' ordinary smooth parameter inside a covariance matrix.
#'
#' Restrictions: normal likelihood only (ordinal latents cannot be
#' integrated out in closed form), mutually exclusive with the GPU program
#' (the work is a Cholesky per pattern, not one large vectorized density),
#' and latent factor scores are not sampled. Identification, the sum-to-zero
#' machinery and the sign relabeling are all unchanged -- marginalization
#' alters only how the likelihood is evaluated, not the model.
#'
#' \strong{Hyper-disturbance parameterization.} Under `hyper_loadings =
#' "unity"` the hyper factor carries a free variance, and when the data pull
#' that variance toward zero the sampler must work against a boundary: the
#' latent disturbance it multiplies becomes unidentified there, and effective
#' sample size for the variance -- and for everything whose scale depends on
#' it -- collapses. Two `bayes_args` switches address this, and they are
#' independent:
#' \itemize{
#'   \item `hyper_scale`: `"tau"` (default) samples the hyper SD directly
#'     against a lower bound of zero. `"rho"` instead samples the implied
#'     within-block correlation on \eqn{(0, 1)} and derives
#'     \eqn{\tau = \sqrt{\rho/(1-\rho)}}. The sampler then works on
#'     \eqn{\mathrm{logit}(\rho)}, a bounded and better-scaled coordinate,
#'     and the prior (`bayes_args$priors$rho_g`, default `beta(2, 8)`, whose
#'     density vanishes at both ends) keeps mass off the degenerate point at
#'     which the hyper factor disappears. Excluding that point is consistent
#'     with a block that the permutation test in [DimSEM_propose()] has
#'     already validated as non-independent, but it is a deliberate prior
#'     choice and should be reported with a sensitivity check.
#'   \item `hyper_param`: `"noncentered"` (default) samples the disturbances
#'     as \eqn{\tau \times} standard normal; `"centered"` samples them on
#'     their own scale. Neither dominates: centering is more efficient when
#'     the likelihood is informative about the disturbance relative to its
#'     prior, non-centering when it is not. Since the prior precision is
#'     \eqn{1/\tau^2}, the comparison depends on where \eqn{\tau} sits --
#'     so the two switches interact, and centering typically pays off only
#'     once `hyper_scale = "rho"` has moved \eqn{\rho} away from zero. Test
#'     `"rho"` alone first, then add centering.
#' }
#' Both switches concern the hyper-factor disturbance only. Under
#' `hyper_loadings = "free"` the hyper variance is fixed at 1, the two
#' parameterizations are the same program, and both settings are ignored with
#' a warning. Within constrained blocks the reported
#' decomposition is additionally made orbit-invariant: when the residual
#' shared variance is negligible, per-factor reflections become
#' near-symmetries of the posterior, so the sign relabeling re-derives
#' the reported coefficients per draw from sign-corrected total effects
#' (reducing exactly to the ordinary correction whenever a block's signs
#' are uniform). This applies under both identification strategies, with
#' the appropriate group action in each: a lone child reflection leaves the
#' hyper-level coefficient invariant while negating that child's total, so
#' under `"prior"` the hyper coefficient's own prior-anchored posterior is
#' left untouched rather than being sharpened by the correction.
#' If the posterior median of `rho_g` for a
#' covariate-targeted block falls below `bayes_args$rho_floor` (default
#' 0.05), a diagnostic warning notes that the hyper level is empirically
#' vacuous conditional on the covariates -- typically because a covariate
#' absorbs the block's common variance -- and recommends interpreting the
#' total effects. The shipped `toy_simple_cov` dataset provides a
#' latent-level covariate with known ground truth for this decomposition
#' (see `attr(toy_simple_cov, "truth")`).
#'
#' @return An object of class `"DimSEM_estimate"`: a list with elements
#'   `engine`, `model_syntax`, `fit` (lavaan object or cmdstanr fit),
#'   `converged`, `factor_cor`, `structural` (covariate coefficient
#'   table), `fit_measures`, `stan` (generated Stan program and data list;
#'   Bayesian engine only), `timing`, and `call`.
#'
#' @export
DimSEM_estimate <- function(proposal,
                            data,
                            covariates = NULL,
                            covariate_targets = "all",
                            model_syntax = NULL,
                            include_hyper = TRUE,
                            hyper_loadings = c("unity", "free"),
                            covariate_identification = c("sum_to_zero",
                                                         "prior"),
                            standardize_covariates = TRUE,
                            engine = c("ml", "bayes"),
                            estimator = "MLR",
                            missing = "ml",
                            lavaan_args = list(),
                            bayes_args = list(),
                            force_recompile = FALSE,
                            progress = TRUE,
                            seed = NULL,
                            verbose = TRUE) {
  engine <- match.arg(engine)
  hyper_loadings <- match.arg(hyper_loadings)
  covariate_identification <- match.arg(covariate_identification)
  t_start <- Sys.time()

  pb <- .dimsem_stage_progress(
    stages = c("Validating inputs", "Assembling model syntax",
               "Preparing data", "Estimating model", "Post-processing"),
    enabled = isTRUE(progress)
  )

  # --- Stage 1: validation ------------------------------------------------
  spec <- .dimsem_extract_proposal(proposal, data,
                                   covariates = covariates,
                                   include_hyper = include_hyper)
  pb$tick()

  # --- Stage 2: model syntax ----------------------------------------------
  if (!is.null(model_syntax)) {
    syntax <- paste(model_syntax, collapse = "\n")
    syntax_source <- "user"
  } else {
    syntax <- .dimsem_assemble_syntax(
      spec,
      covariates = covariates,
      covariate_targets = covariate_targets,
      hyper_loadings = hyper_loadings
    )
    syntax_source <- "proposal"
  }
  pb$tick()

  # --- Stage 3: data ------------------------------------------------------
  covariate_info <- NULL
  if (!is.null(covariates) && length(covariates) > 0) {
    scaled <- .dimsem_scale_covariates(data, covariates,
                                       standardize = standardize_covariates)
    data <- scaled$data
    covariate_info <- scaled$info
    if (isTRUE(verbose)) {
      std_cv <- covariate_info$covariate[covariate_info$standardized]
      bin_cv <- covariate_info$covariate[covariate_info$type == "binary"]
      if (length(std_cv) > 0) {
        message("Standardized covariate(s): ",
                paste(std_cv, collapse = ", "), ".")
      }
      if (length(bin_cv) > 0) {
        message("Binary covariate(s) kept in raw metric: ",
                paste(bin_cv, collapse = ", "), ".")
      }
    }
  }
  est_data <- data[, unique(c(spec$item_names, covariates)), drop = FALSE]
  pb$tick()

  # --- Stage 4: estimation ------------------------------------------------
  if (identical(engine, "ml")) {
    fit_out <- .dimsem_estimate_ml(
      syntax = syntax, data = est_data, estimator = estimator,
      missing = missing, lavaan_args = lavaan_args,
      progress = progress, seed = seed, verbose = verbose
    )
  } else {
    fit_out <- .dimsem_estimate_bayes(
      spec = spec, syntax = syntax, data = est_data,
      covariates = covariates, covariate_targets = covariate_targets,
      covariate_info = covariate_info, hyper_loadings = hyper_loadings,
      covariate_identification = covariate_identification,
      bayes_args = bayes_args, force_recompile = force_recompile,
      seed = seed, verbose = verbose
    )
  }
  pb$tick()

  # --- Stage 5: post-processing -------------------------------------------
  out <- c(
    list(
      engine = engine,
      model_syntax = syntax,
      syntax_source = syntax_source,
      proposal_source = spec$source,
      hyper_blocks = spec$hyper_blocks,
      hyper_loadings = hyper_loadings,
      covariate_identification = covariate_identification,
      covariates = covariates,
      covariate_scaling = covariate_info
    ),
    fit_out,
    list(
      timing = list(started = t_start, finished = Sys.time(),
                    elapsed = difftime(Sys.time(), t_start, units = "secs")),
      call = match.call()
    )
  )
  class(out) <- "DimSEM_estimate"
  pb$tick()
  pb$close()
  out
}

#' @export
print.DimSEM_estimate <- function(x, ...) {
  cat("DimSem estimate (engine:", x$engine, ")\n")
  cat("  Syntax source:", x$syntax_source,
      "| proposal partition source:", x$proposal_source, "\n")
  if (length(x$hyper_blocks) > 0) {
    cat("  Hyper-factors:",
        paste(vapply(seq_along(x$hyper_blocks), function(i) {
          paste0(names(x$hyper_blocks)[i], " = {",
                 paste(x$hyper_blocks[[i]], collapse = ", "), "}")
        }, character(1)), collapse = "; "), "\n")
  }
  if (length(x$hyper_blocks) > 0) {
    cat("  Hyper loadings:",
        if (identical(x$hyper_loadings, "unity")) {
          "fixed to unity (Warncke & Azevedo, n.d.)"
        } else {
          "freely estimated"
        }, "\n")
  }
  if (!is.null(x$covariates)) {
    cat("  Covariates:", paste(x$covariates, collapse = ", "), "\n")
    if (!is.null(x$covariate_scaling)) {
      std_cv <- x$covariate_scaling$covariate[x$covariate_scaling$standardized]
      bin_cv <- x$covariate_scaling$covariate[
        x$covariate_scaling$type == "binary"]
      if (length(std_cv) > 0) {
        cat("    standardized:", paste(std_cv, collapse = ", "), "\n")
      }
      if (length(bin_cv) > 0) {
        cat("    binary (raw metric):", paste(bin_cv, collapse = ", "), "\n")
      }
    }
  }
  cat("  Converged:", isTRUE(x$converged), "\n")
  if (!is.null(x$fit_measures)) {
    fm <- x$fit_measures
    if (isTRUE(attr(fm, "scaled"))) cat("  Robust (scaled) fit statistics:\n")
    cat(sprintf(
      "  chisq(%d) = %.2f, p = %.3f | CFI = %.3f, TLI = %.3f, RMSEA = %.3f, SRMR = %.3f\n",
      as.integer(fm[["df"]]), fm[["chisq"]], fm[["pvalue"]],
      fm[["cfi"]], fm[["tli"]], fm[["rmsea"]], fm[["srmr"]]))
  }
  if (!is.null(x$factor_cor) && ncol(x$factor_cor) > 1) {
    cat("  Latent correlations:\n")
    print(round(x$factor_cor, 3))
  }
  if (!is.null(x$structural) && nrow(x$structural) > 0) {
    cat("  Structural (covariate) coefficients:\n")
    print(x$structural, row.names = FALSE, digits = 3)
  }
  if (!is.null(x$stan) && is.null(x$fit)) {
    cat("  [dry run] Stan program generated (",
        length(strsplit(x$stan$model_code, "\n")[[1]]),
        "lines ); compile with cmdstanr locally.\n")
  }
  cat("  Elapsed:", format(round(as.numeric(x$timing$elapsed), 1)), "s\n")
  invisible(x)
}
