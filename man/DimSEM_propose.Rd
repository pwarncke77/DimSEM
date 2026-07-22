\name{DimSem_propose}
\alias{DimSem_propose}
\alias{print.DimSem_proposal}
\title{Propose a dimensionality-sensitive SEM measurement schematic}
\usage{
DimSem_propose(
  data,
  methods = "EGA",
  item_names = NULL,
  ega_args = list(),
  pa_args = list(),
  tnn_args = list(),
  partition_source = c("auto", "EGA", "PA", "TNN"),
  assignment_args = list(),
  make_plots = TRUE,
  plot_args = list(),
  seed = NULL,
  verbose = TRUE
)

\method{print}{DimSem_proposal}(x, ...)
}
\arguments{
\item{data}{A numeric data frame or matrix with rows as cases/respondents and columns as manifest items.}
\item{methods}{Character vector. Methods to run. Any combination of \code{"EGA"}, \code{"PA"}, \code{"TNN"}, or \code{"all"}. Defaults to \code{"EGA"}.}
\item{item_names}{Optional character vector of item names. Defaults to \code{colnames(data)} or generated names.}
\item{ega_args}{Optional named list of arguments passed to \code{EGAnet::EGA()}. Defaults are merged with \code{list(corr = "auto", model = "glasso", algorithm = "walktrap", uni.method = "LE", plot.EGA = FALSE)}.}
\item{pa_args}{Optional named list of arguments passed to \code{psych::fa.parallel()}. Set \code{n_cores} inside this list to temporarily set \code{options(mc.cores = n)} during the PA call. Defaults are merged with \code{list(fm = "minres", fa = "fa", n.iter = 100, quant = .95, cor = "cor", use = "pairwise", plot = FALSE)}.}
\item{tnn_args}{Optional named list of arguments passed to \code{intRinsic::twonn()}. Convenience arguments intercepted before the call are \code{scale} and \code{missing}. The returned object includes both the raw TNN estimate and \code{ceiling(raw)}.}
\item{partition_source}{Character scalar. Which method should define the final item partition? \code{"auto"} chooses EGA when available, then PA, then TNN.}
\item{assignment_args}{Optional named list of arguments passed to \code{psych::fa()} for fallback item assignment when PA or TNN defines the dimensionality.}
\item{make_plots}{Logical. If \code{TRUE}, create a community/network plot and a directed SEM schematic.}
\item{plot_args}{Optional named list controlling plotting. Supports \code{network_layout}, \code{sem_layout}, \code{node_size}, and \code{label_size}.}
\item{seed}{Optional integer seed applied before stochastic procedures.}
\item{verbose}{Logical. If \code{TRUE}, print method-progress messages.}
\item{x}{A \code{DimSem_proposal} object.}
\item{...}{Ignored.}
}
\description{
Runs one or more dimensionality-detection procedures and returns a proposed item-to-dimension partition that can later be converted into a confirmatory SEM measurement model. By default, the proposal is based on Exploratory Graph Analysis (EGA), which directly returns item communities. If only Parallel Analysis (PA) or Two-Nearest Neighbors (TNN) is requested, the function estimates the number of dimensions and then uses exploratory factor loadings to assign each manifest item to its strongest latent dimension.
}
\value{
An object of class \code{DimSem_proposal}, a list with method-specific results, the selected item partition, lavaan-style measurement syntax, optional plots, and underlying igraph objects.
}
\examples{
\dontrun{
proposal <- DimSem_propose(
  data = my_items,
  methods = c("EGA", "PA", "TNN"),
  pa_args = list(n.iter = 500, n_cores = 8, fa = "fa", plot = FALSE),
  tnn_args = list(method = "mle", scale = TRUE),
  ega_args = list(corr = "auto", model = "glasso", algorithm = "walktrap")
)

proposal
proposal$selected$model_syntax
print(proposal$plots$network)
print(proposal$plots$sem)
}
}
