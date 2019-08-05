% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genscore.R
\name{get_elts_gauss}
\alias{get_elts_gauss}
\title{The R implementation to get the elements necessary for calculations for the untruncated gaussian setting.}
\usage{
get_elts_gauss(x, centered = TRUE, profiled_if_noncenter = TRUE,
  scale = "norm", diagonal_multiplier = 1)
}
\arguments{
\item{x}{A matrix, the data matrix.}

\item{centered}{A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.}

\item{profiled_if_noncenter}{A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if \code{centered==TRUE}. Default to \code{TRUE}.}

\item{scale}{A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.}

\item{diagonal_multiplier}{A number >= 1, the diagonal multiplier.}
}
\value{
A list that contains the elements necessary for estimation.
  \item{n}{The sample size.}
  \item{p}{The dimension.}
  \item{centered}{The centered setting or not. Same as input.}
  \item{scale}{The scaling method. Same as input.}
  \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
  \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
  \item{setting}{The setting \code{"gaussian"}.}
  \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}. Except for the \emph{profiled} setting, this is \eqn{\mathbf{xx}^{\top}/n}{xx'/n}.}
  \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}. The minus column means of \code{x}.}
  \item{t1,t2}{Returned in the profiled non-centered setting, where the\eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
}
\description{
The R implementation to get the elements necessary for calculations for the untruncated gaussian setting.
}
\details{
For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
}
\examples{
if (!requireNamespace("mvtnorm", quietly = TRUE)){
  stop("Please install package \\"mvtnorm\\" first.", call. = FALSE)
}
require(mvtnorm)
n <- 50
p <- 30
h_hp <- get_h_hp("min_pow", 1, 3)
mu <- rep(0, p)
K <- diag(p)
x <- mvtnorm::rmvnorm(n, mean=mu, sigma=solve(K))
get_elts_gauss(x, centered=TRUE, scale="norm", diag=1.5)
get_elts_gauss(x, centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
}