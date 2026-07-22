#' Propose a dimensionality-sensitive SEM schematic
#'
#' `DimSEM_propose()` runs one or more dimensionality-detection procedures and
#' returns a proposed item-to-dimension partition that can later be converted
#' into a hierarchical SEM model. By default, the proposal is based on
#' Exploratory Graph Analysis (EGA), which directly returns item communities. If
#' Parallel Analysis (PA) or the Empirical Kaiser Criterion (EKC) is
#' requested, the function estimates the number of dimensions and then partitions the regularized item
#' partial-correlation network into exactly that many communities using
#' weighted walktrap community detection ([igraph::cluster_walktrap()] +
#' [igraph::cut_at()]) by default; see `assignment_args` for alternatives.
#'
#' @param data A data frame or matrix with rows as cases/respondents. If
#'   `items` is left `NULL`, every column is treated as a manifest item and
#'   must be numeric. If `items` is supplied, `data` may be a raw data set
#'   containing additional (even non-numeric) columns; only the selected item
#'   columns are checked and analyzed.
#' @param items Optional character vector of column names (or integer vector
#'   of column positions) selecting the subset of columns in `data` that
#'   represent manifest items. Defaults to `NULL`, in which case all columns
#'   of `data` are passed downstream.
#' @param methods Character vector determining methods to run. Any combination of `"EGA"`,
#'   `"PA"`, `"EKC"`, or `"all"` is accepted, with `"EGA"` used by default.
#' @param item_names Optional character vector of display names for the
#'   analyzed items. If supplied, must have length equal to the number of
#'   analyzed items (i.e., `length(items)` when `items` is used, otherwise
#'   `ncol(data)`). Defaults to the selected column names or generated names
#'   of the form `item_1`, `item_2`, ...
#' @param ega_args Optional named list of arguments passed to [EGAnet::EGA()].
#'   Defaults are merged with `list(corr = "auto", model = "glasso",
#'   algorithm = "louvain", consensus.method = "lowest_tefi",
#'   uni.method = "LE", plot.EGA = FALSE)`. By default, DimSEM additionally
#'   enforces a minimum community size of three items via Louvain
#'   resolution annealing: whenever the EGA solution contains a community
#'   of fewer than `min_community_size` (default 3) items that is connected
#'   to the rest of the network, the Louvain `resolution` parameter is
#'   lowered from `resolution_start` (default 1.2) in decrements of
#'   `resolution_step` (default 0.01), coarsening the partition, until all
#'   communities reach the minimum size or `resolution_min` (default 0.05)
#'   is hit. Sub-threshold communities that form their own disconnected
#'   network component (isolated items or isolated doublets) can never be
#'   merged by coarsening and are instead converted to unassigned orphans
#'   (NA), as are any communities still below the threshold at the
#'   resolution floor; orphaned items are excluded from the proposed
#'   measurement syntax and SEM schematic. Rationale: a two-item community
#'   cannot identify a latent dimension (the three-indicator rule; Gorsuch,
#'   1997), consistent with the two-node-community handling in Taxonomic
#'   Graph Analysis (Samo et al., 2025) and with the hyper-level
#'   `min_block_size` rule. All four annealing controls are intercepted
#'   convenience arguments and are not passed to [EGAnet::EGA()]. Set
#'   `min_community_size = 1` to disable the rule, or supply a fixed
#'   `resolution` (or a non-Louvain `algorithm`, e.g. `"walktrap"`) to
#'   disable the search while retaining orphan conversion.
#' @param pa_args Optional named list of arguments passed to
#'   [psych::fa.parallel()]. Set `n_cores` inside this list to temporarily set
#'   `options(mc.cores = n)` during the PA call. Defaults are merged with
#'   `list(fm = "minres", fa = "fa", n.iter = 100, quant = .95, cor = "cor",
#'   use = "pairwise", plot = FALSE)`. A second convenience argument,
#'   `suppress_heywood` (default `TRUE`), muffles the "(ultra-)Heywood case"
#'   warnings that [psych::fa.parallel()] emits while deliberately
#'   overextracting factors from observed and resampled correlation matrices;
#'   these throwaway EFA solutions do not affect the eigenvalue-based
#'   retention decision. All other warnings are passed through. Set
#'   `suppress_heywood = FALSE` to surface them.
#' @param ekc_args Optional named list of arguments for the Empirical Kaiser
#'   Criterion (Braeken & van Assen, 2017). Currently supports `use`, passed to
#'   [stats::cor()] (default `"pairwise.complete.obs"`). EKC compares sample
#'   eigenvalues of the item correlation matrix against serially updated
#'   Marchenko-Pastur reference eigenvalues (floored at the classical Kaiser
#'   value of 1) and retains dimensions until the first eigenvalue falls below
#'   its reference. It is deterministic, closed-form, and well suited to
#'   Likert-type item data. A raw retention count of 0 (e.g., pure-noise data)
#'   is reported as such in `results$EKC$n_dim_raw` but clamped to 1 for
#'   partition construction.
#' @param partition_source Character scalar. Which method should define the
#'   final item partition? `"auto"` prefers EGA, then PA, then EKC. `"EGA"`
#'   uses the EGA community vector directly. `"PA"` or `"EKC"` uses the estimated
#'   number of dimensions and assigns items by partitioning the regularized
#'   partial-correlation network into exactly that many communities (weighted
#'   walktrap by default; see `assignment_args`).
#' @param assignment_args Optional named list controlling how PA/EKC
#'   dimension counts are turned into item partitions. Supports `algorithm`,
#'   one of `"walktrap"` (default; [igraph::cluster_walktrap()] with edge
#'   weights, cut at exactly the requested number of communities via
#'   [igraph::cut_at()]), `"fast_greedy"` ([igraph::cluster_fast_greedy()]
#'   with weights, also cut via `cut_at()`), or `"fluid"`
#'   ([igraph::cluster_fluid_communities()]; ignores edge weights and is
#'   retained only for backward compatibility); `glasso_args` for arguments passed to
#'   [qgraph::EBICglasso()] when an EGA network is unavailable; `cor = "auto"`
#'   or `"pearson"` for the correlation matrix used by the glasso fallback;
#'   `use` for `stats::cor()` when `cor = "pearson"`; `weight_transform`
#'   (`"absolute"` or `"positive"`) for turning the regularized signed matrix
#'   into graph topology; and `connect_disconnected`, a logical that allows tiny
#'   bridge edges to be added when the regularized network is disconnected.
#' @param fit_cfa Logical. If `TRUE` (default), fit a confirmatory factor
#'   analysis with [lavaan::cfa()] using the proposed item partition as the
#'   measurement model and store it in the `cfa` element of the returned
#'   object, together with the estimated inter-factor correlation matrix and
#'   standard fit measures. Requires the `lavaan` package.
#' @param cfa_source Character scalar. Which item partition should define the
#'   CFA measurement model? `"selected"` (default) uses the partition already
#'   chosen via `partition_source` -- under auto-selection this is the EGA
#'   solution whenever EGA was run and produced item communities. `"EGA"`,
#'   `"PA"`, or `"EKC"` instead force the CFA onto that specific method's
#'   partition, which may differ from the reported `selected` partition.
#' @param cfa_args Optional named list of arguments passed to [lavaan::cfa()].
#'   Defaults are merged with `list(std.lv = TRUE, ordered = TRUE, missing = "pairwise")`:
#'   `std.lv = TRUE` fixes latent variances to 1 so that all loadings are
#'   estimated and latent covariances are correlations; `ordered = TRUE`
#'   treats all items as ordered-categorical (Likert-type or binary), making
#'   lavaan estimate the model from polychoric/tetrachoric correlations with
#'   a WLSMV-type estimator, and the reported fit measures are then the
#'   robust (scaled) versions. Pass `list(ordered = FALSE)` to treat items as
#'   continuous instead. Items left unassigned by the partition (NA
#'   dimension) are excluded from the measurement model.
#' @param propose_hyper Logical. If `TRUE` (default), evaluate whether the
#'   proposed level-1 dimensions warrant one or more second-order
#'   ("hyper") factors, using network-psychometric criteria only. Results
#'   are stored in the `hyper` element of the returned object.
#' @param hyper_args Optional named list controlling the hyper-structure
#'   evaluation. `method` selects the strategy: `"tga"` runs a two-level
#'   adaptation of Taxonomic Graph Analysis (Samo et al., 2025): network
#'   loadings and scores are computed from the level-1 network and partition
#'   (Christensen et al., 2025; Golino et al., 2022), a level-2 EBICglasso
#'   network is estimated on the score correlations, lower-order Louvain
#'   community detection with most-common consensus identifies level-2
#'   communities (Blondel et al., 2008; Jimenez et al., 2023; Lancichinetti
#'   & Fortunato, 2012), and the unidim index (Revelle & Condon, 2025)
#'   statistically gates the single-hyper-factor verdict. `"permutation"`
#'   runs a permutation-validated higher-order Louvain: the null hypothesis
#'   that level-1 communities are mutually independent conditional on their
#'   internal structure is simulated exactly by permuting data rows
#'   independently within each community's column block; the
#'   between-community edge mass of the re-estimated regularized network
#'   (the between-community component of the weighted-modularity numerator)
#'   is the test statistic; BH-FDR-validated community pairs form a
#'   meta-graph on which Louvain consensus identifies hyper-blocks, with
#'   unconnected dimensions remaining stand-alone. `"both"` (default) runs
#'   the two strategies and reports their agreement. Further options:
#'   `n_perm` permutations (default: adaptive,
#'   `max(500, ceiling(2 * n_pairs / alpha))`, so that a single true
#'   between-community connection can survive FDR correction);
#'   `alpha` (default .05) FDR level
#'   and omnibus gate; `min_block_size` (default 3) the minimum number of
#'   level-1 dimensions required to support a hyper-factor -- a hyper-factor
#'   over two dimensions is an exact, zero-df reparameterization of their
#'   correlation and hence unfalsifiable, so validated groupings below this
#'   size are reported as correlated level-1 dimensions instead;
#'   `member_rule` (`"majority"`, the default, or `"all"`) the
#'   family-resemblance membership criterion in the permutation strategy: a
#'   dimension is retained in a candidate hyper-community only if its
#'   MARGINAL between-dimension dependence with the required share of
#'   co-members is individually permutation-validated, with weakly attached
#'   "stray" dimensions removed by backward elimination and reported in
#'   `$hyper$permutation$evicted`; `u_hi`/`u_lo` (defaults .90/.50) unidim verdict
#'   thresholds per Revelle & Condon's simulations, with the interval
#'   between them reported as inconclusive; `consensus_reps` (default 500)
#'   Louvain consensus repetitions; `resolution` (default 1) Louvain
#'   resolution; `glasso_args` passed to [qgraph::EBICglasso()] at level 2.
#' @param parallel_hyper Logical. If `TRUE` (default), evaluate the
#'   permutation replicates of the hyper-structure test in parallel on a
#'   `doSNOW` cluster. The permutation loop is the dominant cost of
#'   `DimSem_propose()` (each replicate re-estimates an [qgraph::EBICglasso()]
#'   network, and the default `n_perm` grows with the number of level-1
#'   dimension pairs), and replicates are mutually independent, so the
#'   speed-up is close to linear in the number of cores. Requires the
#'   `foreach`, `doSNOW`, and `parallel` packages; if any is unavailable, or
#'   if only one core is usable, the test falls back to serial evaluation
#'   with a warning. RESULTS ARE UNAFFECTED: all permutation indices are
#'   drawn serially from the master RNG stream before any work is
#'   dispatched, and workers perform no random number generation, so the
#'   output for a given `seed` is bit-identical to serial evaluation
#'   regardless of `parallel_hyper` or the number of cores. The TGA strategy
#'   is not parallelized (it fits a single level-2 network and is cheap by
#'   comparison).
#' @param parallel_hyper_ncores Integer or `NULL`. Number of worker
#'   processes to use when `parallel_hyper = TRUE`. Defaults to `NULL`,
#'   which resolves to `parallel::detectCores() - 1`, leaving one core for
#'   background tasks. The value is clamped to the number of available
#'   cores and never exceeds the number of permutation replicates.
#' @param make_plots Logical. If `TRUE`, create a community/network plot and a
#'   simple directed SEM schematic.
#' @param plot_args Optional named list controlling plotting. Currently supports
#'   `network_layout`, `sem_layout`, `node_size`, and `label_size`.
#' @param seed Optional integer seed applied before stochastic procedures.
#' @param verbose Logical. If `TRUE`, print method-progress messages.
#'
#' @return An object of class `"DimSEM_proposal"`, a list with the following
#'   elements:
#'   \describe{
#'     \item{results}{Raw and summarized method-specific results for EGA, PA, and/or EKC. PA/EKC results include their graph-partition item allocations.}
#'     \item{selected}{The selected partition source, dimension count, item partition table, and lavaan-style measurement syntax.}
#'     \item{cfa}{When `fit_cfa = TRUE`: the partition source used, the model syntax, the fitted [lavaan::cfa()] object (`fit`), a convergence flag, the inter-factor correlation matrix (`factor_cor`), and selected fit measures. `NULL` when `fit_cfa = FALSE`.}
#'     \item{hyper}{When `propose_hyper = TRUE`: per-strategy hyper-structure results (verdict, hyper-blocks, stand-alone dimensions, unidim values, permutation pair table and global p-value), an agreement flag, and a plain-language recommendation. `NULL` when `propose_hyper = FALSE`.}
#'     \item{plots}{A list containing `network` and `sem` ggplot/ggraph objects when `make_plots = TRUE`.}
#'     \item{graphs}{Underlying igraph objects for the item network and SEM schematic.}
#'     \item{data_info}{Basic information about the analyzed item matrix.}
#'   }
#'
#' @examples
#' \dontrun{
#' proposal <- DimSEM_propose(
#'   data = my_items,
#'   items = c("q1", "q2", "q3", "q4", "q5", "q6"),
#'   methods = c("EGA", "PA", "EKC"),
#'   pa_args = list(n.iter = 500, n_cores = 8, fa = "fa", plot = FALSE),
#'   ega_args = list(corr = "auto", model = "glasso", algorithm = "walktrap")
#' )
#'
#' proposal
#' proposal$selected$model_syntax
#' summary(proposal$cfa$fit, fit.measures = TRUE)
#' proposal$cfa$factor_cor
#' proposal$results$PA$partition_table
#' print(proposal$plots$network)
#' print(proposal$plots$sem)
#' }
#'
#' @export
DimSEM_propose <- function(data,
                           items = NULL,
                           methods = "EGA",
                           item_names = NULL,
                           ega_args = list(),
                           pa_args = list(),
                           ekc_args = list(),
                           partition_source = c("auto", "EGA", "PA", "EKC"),
                           assignment_args = list(),
                           fit_cfa = TRUE,
                           cfa_source = c("selected", "EGA", "PA", "EKC"),
                           cfa_args = list(),
                           propose_hyper = TRUE,
                           hyper_args = list(),
                           parallel_hyper = TRUE,
                           parallel_hyper_ncores = NULL,
                           make_plots = TRUE,
                           plot_args = list(),
                           seed = NULL,
                           verbose = TRUE) {
  partition_source <- match.arg(partition_source)
  x <- .dimsem_prepare_data(data, items = items, item_names = item_names)
  methods <- .dimsem_normalize_methods(methods)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  results <- list()

  if ("EGA" %in% methods) {
    if (isTRUE(verbose)) message("Running EGA...")
    results$EGA <- .dimsem_run_ega(x$data, ega_args = ega_args, seed = seed,
                                   verbose = verbose)
  }

  if ("PA" %in% methods) {
    if (isTRUE(verbose)) message("Running Parallel Analysis...")
    results$PA <- .dimsem_run_pa(x$data, pa_args = pa_args, seed = seed)
  }

  if ("EKC" %in% methods) {
    if (isTRUE(verbose)) message("Running Empirical Kaiser Criterion...")
    results$EKC <- .dimsem_run_ekc(x$data, ekc_args = ekc_args, seed = seed)
  }

  network_matrix <- .dimsem_get_network_matrix(
    data = x$data,
    results = results,
    item_names = x$item_names,
    assignment_args = assignment_args
  )

  results <- .dimsem_add_graph_partitions(
    results = results,
    network_matrix = network_matrix,
    item_names = x$item_names,
    assignment_args = assignment_args,
    seed = seed,
    verbose = verbose
  )

  selected_source <- .dimsem_select_partition_source(
    requested = partition_source,
    results = results
  )

  selected <- .dimsem_build_selected_partition(
    item_names = x$item_names,
    results = results,
    selected_source = selected_source
  )

  cfa_source <- match.arg(cfa_source)
  cfa <- NULL
  if (isTRUE(fit_cfa)) {
    if (isTRUE(verbose)) {
      message("Fitting CFA (lavaan) on the ",
              if (identical(cfa_source, "selected")) selected$source
              else cfa_source,
              " partition...")
    }
    cfa <- .dimsem_fit_cfa(
      data = x$data,
      results = results,
      selected = selected,
      cfa_source = cfa_source,
      cfa_args = cfa_args
    )
  }

  hyper <- NULL
  if (isTRUE(propose_hyper)) {
    hyper <- .dimsem_hyper_proposal(
      data = x$data,
      network_matrix = network_matrix,
      partition = selected$partition,
      hyper_args = hyper_args,
      seed = seed,
      verbose = verbose,
      parallel = parallel_hyper,
      ncores = parallel_hyper_ncores
    )
  }

  graphs <- list(
    network = .dimsem_network_graph(
      network_matrix = network_matrix,
      partition = selected$partition,
      hyper = hyper
    ),
    sem = .dimsem_sem_graph(
      partition_table = selected$partition_table,
      hyper = hyper
    )
  )

  plots <- list(network = NULL, sem = NULL)
  if (isTRUE(make_plots)) {
    plots$network <- .dimsem_plot_network(graphs$network, plot_args = plot_args)
    plots$sem <- .dimsem_plot_sem(graphs$sem, plot_args = plot_args)
  }

  out <- list(
    call = match.call(),
    methods = methods,
    results = results,
    selected = selected,
    cfa = cfa,
    hyper = hyper,
    plots = plots,
    graphs = graphs,
    data_info = list(
      n = nrow(x$data),
      p = ncol(x$data),
      item_names = x$item_names
    )
  )

  class(out) <- "DimSEM_proposal"
  out
}

#' @export
print.DimSEM_proposal <- function(x, ...) {
  cat("DimSEM proposal\n")
  cat("  Cases:", x$data_info$n, "\n")
  cat("  Items:", x$data_info$p, "\n")
  cat("  Selected source:", x$selected$source, "\n")
  cat("  Selected dimensions:", x$selected$n_dim, "\n\n")

  cat("Method summaries:\n")
  if (!is.null(x$results$EGA)) {
    cat("  EGA:", x$results$EGA$n_dim, "dimension(s)\n")
  }
  if (!is.null(x$results$PA)) {
    cat(
      "  PA:",
      x$results$PA$n_dim,
      "dimension(s); graph partition:",
      x$results$PA$n_dim_partition,
      "community/communities\n"
    )
  }
  if (!is.null(x$results$EKC)) {
    cat(
      "  EKC:",
      x$results$EKC$n_dim_raw,
      "dimension(s) retained; graph partition:",
      x$results$EKC$n_dim_partition,
      "community/communities\n"
    )
  }

  if (!is.null(x$cfa)) {
    cat("\nCFA (lavaan), based on", x$cfa$source, "partition:\n")
    if (!is.null(x$cfa$error)) {
      cat("  Not fitted:", x$cfa$error, "\n")
    } else if (!isTRUE(x$cfa$converged)) {
      cat("  Estimation did not converge; see `$cfa$fit`.\n")
    } else {
      fm <- x$cfa$fit_measures
      if (!is.null(fm)) {
        if (isTRUE(attr(fm, "scaled"))) {
          cat("  Robust (scaled) fit statistics:\n")
        }
        cat(sprintf(
          "  chisq(%d) = %.2f, p = %.3f | CFI = %.3f, TLI = %.3f, RMSEA = %.3f, SRMR = %.3f\n",
          as.integer(fm[["df"]]), fm[["chisq"]], fm[["pvalue"]],
          fm[["cfi"]], fm[["tli"]], fm[["rmsea"]], fm[["srmr"]]
        ))
      }
      if (!is.null(x$cfa$factor_cor) && ncol(x$cfa$factor_cor) > 1) {
        cat("  Inter-factor correlations:\n")
        print(round(x$cfa$factor_cor, 3))
      } else {
        cat("  Single latent factor: no inter-factor correlations.\n")
      }
    }
  }

  if (!is.null(x$hyper)) {
    cat("\nHyper-structure (network-based):\n")
    fmt_blocks <- function(res) {
      if (length(res$blocks) == 0) return("-")
      paste(vapply(seq_along(res$blocks), function(i) {
        paste0(names(res$blocks)[i], " = {",
               paste(res$blocks[[i]], collapse = ", "), "}")
      }, character(1)), collapse = "; ")
    }
    if (!is.null(x$hyper$tga)) {
      cat("  TGA-style:", x$hyper$tga$verdict,
          "| blocks:", fmt_blocks(x$hyper$tga),
          "| unidim u =", round(x$hyper$tga$unidim, 3), "\n")
    }
    if (!is.null(x$hyper$permutation)) {
      cat("  Permutation:", x$hyper$permutation$verdict,
          "| blocks:", fmt_blocks(x$hyper$permutation),
          "| global p =", round(x$hyper$permutation$global_p, 4), "\n")
    }
    cat("  Recommendation:", x$hyper$recommendation, "\n")
  }

  cat("\nItem allocation:\n")
  print(x$selected$partition_table, row.names = FALSE)

  invisible(x)
}
