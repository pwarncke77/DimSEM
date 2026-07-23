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
#'   `list(G1 = c("age", "educ"), F2 = "age")`.
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
#'   (1000), `iter_sampling` (1000), `adapt_delta` (.9), `refresh`.
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
