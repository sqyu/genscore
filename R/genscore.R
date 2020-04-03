## Package:
#library(knitr); library(rmarkdown); library(devtools); library(roxygen2)
#tools::package_native_routine_registration_skeleton(".") // Copy as src/genscore_init.c and manually add dist, rab_arms
#document(); build(); install(); devtools::check()
#source("vignettes/precompile.R")
#check_win_devel(); check_win_release(); check_rhub()
#rhub::check_on_macos()


#' Finds the distance of each element in a matrix x to the its boundary of the domain while fixing the others in the same row (dist(x, domain)), and calculates element-wise h(dist(x, domain)) and h\'(dist(x, domain)) (w.r.t. each element in x).
#'
#' Finds the distance of each element in a matrix  \code{x} to its boundary of the \code{domain} while fixing the others in the same row (\code{dist(x, domain)}), and calculates element-wise \code{h(dist(x, domain))} and \code{h\'(dist(x, domain))} (w.r.t. each element in \code{x}).
#'
#' @param h_hp A function, the \eqn{h} and \eqn{hp} (the derivative of \code{h}) functions. \code{h_hp(x)} should return a list of elements \code{hx} (\code{h(x)}) and \code{hpx} (\code{hp(x)}), both of which have the same size as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @return A list that contains \code{h(dist(x, domain))} and \code{h\'(dist(x, domain))}.
#'   \item{hdx}{\code{h(dist(x, domain))}.}
#'   \item{hpdx}{\code{hp(dist(x, domain))}.}
#' @details
#' Define \code{dist(x, domain)} as the matrix whose \code{i,j}-th component is the distance of \eqn{x_{i,j}} to the boundary of \code{domain}, assuming \eqn{x_{i,-j}} are fixed. The matrix has the same size of \code{x} (\code{n} by \code{p}), or if \code{domain$type == "simplex"} and \code{x} has full dimension \code{p}, it has \code{p-1} columns.\cr
#' Define \code{dist\'(x, domain)} as the component-wise derivative of \code{dist(x, domain)} in its components. That is, its \code{i,j}-th component is 0 if \eqn{x_{i,j}} is unbounded or is bounded from both below and above or is at the boundary, or -1 if \eqn{x_{i,j}} is closer to its lower boundary (or if its bounded from below but unbounded from above), or 1 otherwise.\cr
#' \code{h_of_dist(h_hp, x, domain)} simply returns \code{h_hp(dist(x, domain))$hx} and \code{h_hp(dist(x, domain))$hpx * dist\'(x, domain)} (element-wise derivative of \code{h_hp(dist(x, domain))$hx} w.r.t. \code{x}).
#' @examples
#' n <- 20
#' p <- 10
#' eta <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#'
#' # Gaussian on R^p:
#' domain <- make_domain("R", p=p)
#' x <- mvtnorm::rmvnorm(n, mean=solve(K, eta), sigma=solve(K))
#' # Equivalently:
#' \donttest{
#' x2 <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, burn_in=1000, thinning=100)
#' }
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx is all Inf and hpdx is all 0 since each coordinate is unbounded with domain R
#' c(all(is.infinite(hd$hdx)), all(hd$hpdx==0))
#'
#'
#' # exp on R_+^p:
#' domain <- make_domain("R+", p=p)
#' x <- gen(n, setting="exp", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx is x^2 and hpdx is 2*x; with domain R+, the distance of x to the boundary is just x itself
#' c(max(abs(hd$hdx - x^2)), max(abs(hd$hpdx - 2*x)))
#'
#'
#' # Gaussian on sum(x^2) > p with x allowed to be negative
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste("sum(x^2)>", p), abs=FALSE, nonnegative=FALSE)))
#' x <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' quota <- p - (rowSums(x^2) - x^2) # How much should xij^2 at least be so that sum(xi^2) > p?
#' # How far is xij from +/-sqrt(quota), if quota >= 0?
#' dist_to_bound <- abs(x[quota >= 0]) - abs(sqrt(quota[quota >= 0]))
#' # Should be equal to our own calculations
#' max(abs(dist$dx[is.finite(dist$dx)] - dist_to_bound))
#' # dist'(x) should be the same as the sign of x
#' all(dist$dpx[is.finite(dist$dx)] == sign(x[quota >= 0]))
#' # quota is negative <-> sum of x_{i,-j}^2 already > p <-> xij unbounded given others
#' #      <-> distance to boundary is Inf
#' all(quota[is.infinite(dist$dx)] < 0)
#'
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' # Now confirm that h_of_dist indeed applies h and hp to dists
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx = dist ^ 2
#' print(max(abs(hd$hdx[is.finite(dist$dx)] - dist$dx[is.finite(dist$dx)]^2)))
#' # hdx = Inf if dist = Inf
#' print(all(is.infinite(hd$hdx[is.infinite(dist$dx)])))
#'  # hpdx = 2 * dist' * dist
#' print(max(abs(hd$hpdx[is.finite(dist$dx)] - 2*(dist$dpx*dist$dx)[is.finite(dist$dx)])))
#' print(all(hd$hpdx[is.infinite(dist$dx)] == 0)) # hpdx = 0 if dist = Inf
#'
#'
#' # gamma on ([0, 1] v [2,3])^p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#' x <- gen(n, setting="gamma", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' # If 0 <= xij <= 1, distance to boundary is min(x-0, 1-x)
#' max(abs(dist$dx - pmin(x, 1-x))[x >= 0 & x <= 1])
#' # If 0 <= xij <= 1, dist'(xij) is 1 if it is closer to 0, or -1 if it is closer 1,
#' #   assuming xij %in% c(0, 0.5, 1) with probability 0
#' all((dist$dpx == 2 * (1-x > x) - 1)[x >= 0 & x <= 1])
#' # If 2 <= xij <= 3, distance to boundary is min(x-2, 3-x)
#' max(abs(dist$dx - pmin(x-2, 3-x))[x >= 2 & x <= 3])
#' # If 2 <= xij <= 3, dist'(xij) is 1 if it is closer to 2, or -1 if it is closer 3,
#' #   assuming xij %in% c(2, 2.5, 3) with probability 0
#' all((dist$dpx == 2 * (3-x > x-2) - 1)[x >= 2 & x <= 3])
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' # Now confirm that h_of_dist indeed applies h and hp to dists
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx = dist ^ 2
#' print(max(abs(hd$hdx - dist$dx^2)))
#' # hpdx = 2 * dist' * dist
#' print(max(abs(hd$hpdx - 2*dist$dpx*dist$dx)))
#'
#'
#' # a0.6_b0.7 on {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
#' set.seed(1)
#' xinit <- c(1.5, 0.5, abs(stats::rnorm(p-2)) + log(1.3))
#' x <- gen(n, setting="ab_3/5_7/10", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' # x_{i1} has uniform bound [1, +Inf), so its distance to its boundary is x_{i1} - 1
#' max(abs(dist$dx[,1] - (x[,1] - 1)))
#' # x_{i2} has uniform bound [log(1.3), 1], so its distance to its boundary
#' #    is min(x_{i2} - log(1.3), 1 - x_{i2})
#' max(abs(dist$dx[,2] - pmin(x[,2] - log(1.3), 1 - x[,2])))
#' # x_{ij} for j >= 3 has uniform bound [log(1.3), +Inf), so its distance to its boundary is
#' #    simply x_{ij} - log(1.3)
#' max(abs(dist$dx[,3:p] - (x[,3:p] - log(1.3))))
#' # dist'(xi2) is 1 if it is closer to log(1.3), or -1 if it is closer 1,
#' #    assuming x_{i2} %in% c(log(1.3), (1+log(1.3))/2, 1) with probability 0
#' all((dist$dpx[,2] == 2 * (1 - x[,2] > x[,2] - log(1.3)) - 1))
#' all(dist$dpx[,-2] == 1) # x_{ij} for j != 2 is bounded from below but unbounded from above,
#' #    so dist'(xij) is always 1
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' # Now confirm that h_of_dist indeed applies h and hp to dists
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx = dist ^ 2
#' print(max(abs(hd$hdx - dist$dx^2)))
#' # hpdx = 2 * dist' * dist
#' print(max(abs(hd$hpdx - 2*dist$dpx*dist$dx)))
#'
#'
#' # log_log model on {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' x <- gen(n, setting="log_log", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' # Upper bound for j * xij so that sum_j j * xij <= 1
#' quota <- 1 - (rowSums(t(t(x) * 1:p)) - t(t(x) * 1:p))
#' # Distance of xij to its boundary is min(xij - 0, quota_{i,j} / j - xij)
#' max(abs(dist$dx - pmin((t(t(quota) / 1:p) - x), x)))
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' # Now confirm that h_of_dist indeed applies h and hp to dists
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx = dist ^ 2
#' print(max(abs(hd$hdx - dist$dx^2)))
#' # hpdx = 2 * dist' * dist
#' print(max(abs(hd$hpdx - 2*dist$dpx*dist$dx)))
#'
#'
#' # log_log_sum0 model on the simplex with K having row and column sums 0 (Aitchison model)
#' domain <- make_domain("simplex", p=p)
#' K <- -cov_cons("band", p=p, spars=3, eig=1)
#' diag(K) <- diag(K) - rowSums(K) # So that rowSums(K) == colSums(K) == 0
#' eigen(K)$val[(p-1):p] # Make sure K has one 0 and p-1 positive eigenvalues
#' x <- gen(n, setting="log_log_sum0", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' # Note that dist$dx and dist$dpx only has p-1 columns -- excluding the last coordinate in x
#' dist <- get_dist(x, domain)
#' # Upper bound for x_{i,j} so that x_{i,1} + ... + x_{i,p-1} <= 1
#' quota <- 1 - (rowSums(x[,-p]) - x[,-p])
#' # Distance of x_{i,j} to its boundary is min(xij - 0, quota_{i,j} - xij)
#' max(abs(dist$dx - pmin(quota - x[,-p], x[,-p])))
#' h_hp <- get_h_hp("pow", 2) # For demonstration only
#' # Now confirm that h_of_dist indeed applies h and hp to dists
#' hd <- h_of_dist(h_hp, x, domain)
#' # hdx = dist ^ 2
#' print(max(abs(hd$hdx - dist$dx^2)))
#' # hpdx = 2 * dist' * dist
#' print(max(abs(hd$hpdx - 2*dist$dpx*dist$dx)))
#' @export
h_of_dist <- function(h_hp, x, domain) {
  dists <- get_dist(x, domain)
  dx <- dists$dx # Distance to boundary
  dpx <- dists$dpx # Derivative of distance to boundary, 1 if closer to left boundary, -1 if right, 0 if equal or at boundary of domain is entire R
  h_hp_dx <- h_hp(dx) # h(dist to boundary) and h'(dist to boundary)
  hdx <- h_hp_dx$hx # h(dist to boundary)

  hpdx <- h_hp_dx$hpx # h'(dist to boundary(x)) = h'(dist to boundary) * dist'(x)
  hpdx[dpx == 0] <- 0 # Instead of doing h_hp_dx$hpx * dpx, this would avoid Inf * 0 = NaN
  hpdx[dpx == -1] <- -hpdx[dpx == -1]
  return (list("hdx"=hdx, "hpdx"=hpdx))
}

#' Parses an ab setting into rational numbers a and b.
#'
#' Parses an ab setting into rational numbers a and b.
#'
#' @param s A string starting with "ab_", followed by rational numbers a and b separated by "_". a and b must be integers or rational numbers of the form "int/int". See examples.
#' @return A list of the following elements:
#'    \item{a_numer}{The numerator of \code{a}.}
#'    \item{a_denom}{The denominator of \code{a}.}
#'    \item{b_numer}{The numerator of \code{b}.}
#'    \item{b_denom}{The denominator of \code{b}.}
#' @examples
#' parse_ab("ab_1_1") # gaussian: a = 1, b = 1
#' parse_ab("ab_2_5/4") # a = 2, b = 5/4
#' parse_ab("ab_5/4_3/2") # a = 5/4, b = 3/2
#' parse_ab("ab_3/2_0/0") # a = 3/2, b = 0/0 (log)
#' parse_ab("ab_1/2_0/0") # exp: a = 1/2, b = 0/0 (log)
#' @export
parse_ab <- function(s) {
  if (!startsWith(s, "ab_"))
    stop("To be parsed as an ab setting, s must start with ab_. Got '", s, "'.")
  s <- substr(s, 4, nchar(s))
  if (!grepl("_", s))
    stop("a and b must be separated by an underscore _. For example, ab_2_2, ab_2_5/4 or ab_2/3_1/2.")
  ab <- strsplit(s, "_")[[1]]
  get_numer_denom <- function(ss) {
    if (grepl("/", ss)) {
      numer_denom <- as.integer(strsplit(ss, "/")[[1]])
      if (any(is.na(numer_denom)))
        stop("a and b must both be integers or two integers separated by /. For example, ab_2_2, ab_2_5/4 or ab_2/3_1/2.")
      return (numer_denom)
    } else {
      numer <- as.integer(ss)
      if (is.na(numer))
        stop("a and b must both be integers or two integers separated by /. For example, ab_2_2, ab_2_5/4 or ab_2/3_1/2.")
      return (c(numer, 1))
    }
  }
  a <- get_numer_denom(ab[1]); b <- get_numer_denom(ab[2])
  return (list("a_numer"=a[1], "a_denom"=a[2], "b_numer"=b[1], "b_denom"=b[2]))
}

#' The R implementation to get the elements necessary for calculations for general a and b.
#'
#' The R implementation to get the elements necessary for calculations for general \eqn{a} and \eqn{b}.
#'
#' @param hdx A matrix, \eqn{h(\mathbf{x})}{h(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param hpdx A matrix, \eqn{h'(\mathbf{x})}{h\'(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param a A number, must be strictly larger than \eqn{b/2}.
#' @param b A number, must be >= 0.
#' @param setting A string that indicates the distribution type. Returned without being checked or used in the function body.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if \code{centered=TRUE}. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The setting. Same as input.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details
#' Computes the \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix and the \eqn{\boldsymbol{g}}{g} vector for generalized score matching.
#'
#' Here, \eqn{\boldsymbol{\Gamma}}{\Gamma} is block-diagonal, and in the non-profiled non-centered setting, the \eqn{j}-th block is composed of \eqn{\boldsymbol{\Gamma}_{\mathbf{KK},j}}{\Gamma_{KK,j}}, \eqn{\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}}{\Gamma_{K\eta,j}} and its transpose, and finally \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta},j}}{\Gamma_{\eta\eta,j}}. In the centered case, only \eqn{\boldsymbol{\Gamma}_{\mathbf{KK},j}}{\Gamma_{KK,j}} is computed. In the profiled non-centered case, \deqn{\boldsymbol{\Gamma}_{j}\equiv\boldsymbol{\Gamma}_{\mathbf{KK},j}-\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}\boldsymbol{\Gamma}_{\boldsymbol{\eta}\boldsymbol{\eta},j}^{-1}\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta}}^{\top}.}{\Gamma_j=\Gamma_{KK,j}-\Gamma_{K\eta,j}\Gamma_{\eta\eta,j}^(-1)\Gamma_{K\eta}'.}
#'
#' Similarly, in the non-profiled non-centered setting, \eqn{\boldsymbol{g}}{g} can be partitioned \eqn{p} parts, each with a \eqn{p}-vector \eqn{\boldsymbol{g}_{\mathbf{K},j}}{g_{K,j}} and a scalar \eqn{g_{\boldsymbol{\eta},j}}{g_{\eta,j}}. In the centered setting, only \eqn{\boldsymbol{g}_{\mathbf{K},j}}{g_{K,j}} is needed. In the profiled non-centered case, \deqn{\boldsymbol{g}_j\equiv\boldsymbol{g}_{\mathbf{K},j}-\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta},j}^{-1}g_{\boldsymbol{\eta},j}.}{g_j=g_{K,j}-\Gamma_{K\eta,j}\Gamma_{\eta\eta,j}^(-1)g_{\eta,j}.}
#'
#' The formulae for the pieces above are
#' \deqn{\boldsymbol{\Gamma}_{\mathbf{KK},j}\equiv\frac{1}{n}\sum_{i=1}^nh\left(X_j^{(i)}\right){X_j^{(i)}}^{2a-2}{\boldsymbol{X}^{(i)}}^a{{\boldsymbol{X}^{(i)}}^a}^{\top},}{\Gamma_{KK,j}=1/n*\sum_{i=1}^n h(Xij)*Xij^(2a-2)*Xi^a*(Xi^a)',}
#' \deqn{\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}\equiv-\frac{1}{n}\sum_{i=1}^nh\left(X_j^{(i)}\right){X_j^{(i)}}^{a+b-2}{\boldsymbol{X}^{(i)}}^a,}{\Gamma_{K\eta,j}=-1/n*\sum_{i=1}^n h(Xij)*Xij^(a+b-2)*Xi^a,}
#' \deqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta},j}\equiv\frac{1}{n}\sum_{i=1}^nh\left(X_j^{(i)}\right){X_j^{(i)}}^{2b-2},}{\Gamma_{\eta\eta,j}=1/n*\sum_{i=1}^n h(Xij)*Xij^(2b-2),}
#' \deqn{\boldsymbol{g}_{\mathbf{K},j}\equiv\frac{1}{n}\sum_{i=1}^n\left(h'\left(X_j^{(i)}\right){X_j^{(i)}}^{a-1}+(a-1)h\left(X_j^{(i)}\right){X_j^{(i)}}^{a-2}\right){\boldsymbol{X}^{(i)}}^a+ah\left(X_j^{(i)}\right){X_j^{(i)}}^{2a-2}\boldsymbol{e}_{j,p},}{g_{K,j}=1/n*\sum_{i=1}^n (h'(Xij)*Xij^(a-1)+(a-1)*h(Xij)*Xij^(a-2))*Xi^a+a*h(Xij)*Xij^(2a-2)*e_{j,p},}
#' \deqn{\boldsymbol{g}_{\boldsymbol{\eta},j}\equiv\frac{1}{n}\sum_{i=1}^n-h'\left(X_j^{(i)}\right){X_j^{(i)}}^{b-1}-(b-1)h\left(X_j^{(i)}\right){X_j^{(i)}}^{b-2},}{g_{\eta,j}=1/n*\sum_{i=1}^n -h'(Xij)*Xij^(b-1)-(b-1)*h(Xij)*Xij^(b-2)),}
#' where \eqn{\boldsymbol{e}_{j,p}}{e_{j,p}} is the \eqn{p}-vector with 1 at the \eqn{j}-th position and 0 elsewhere.
#'
#' In the profiled non-centered setting, the function also returns \eqn{t_1}{t1} and \eqn{t_2}{t2} defined as
#' \deqn{\boldsymbol{t}_1\equiv\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta}}^{-1}\boldsymbol{g}_{\boldsymbol{\eta}},\quad\boldsymbol{t}_2\equiv\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta}}^{-1}\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta}}^{\top},}{t1=\Gamma_{\eta\eta}^(-1)g_{\eta}, t2=\Gamma_{\eta\eta}^(-1)\Gamma_{K\eta}',}
#' so that \eqn{\hat{\boldsymbol{\eta}}=\boldsymbol{t}_1-\boldsymbol{t}_2\mathrm{vec}(\hat{\mathbf{K}}).}{\hat{\eta}=t1-t2*vec(\hat{K}).}
#' @examples
#' n <- 50
#' p <- 30
#' eta <- rep(0, p)
#' K <- diag(p)
#' domain <- make_domain("R+", p=p)
#' x <- gen(n, setting="ab_1/2_7/10", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("min_pow", 1.5, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#' get_elts_ab(h_hp_dx$hdx, h_hp_dx$hpdx, x, a=0.5, b=0.7, setting="ab_1/2_7/10",
#'             centered=TRUE, scale="norm", diag=1.5)
#' get_elts_ab(h_hp_dx$hdx, h_hp_dx$hpdx, x, a=0.5, b=0.7, setting="ab_1/2_7/10",
#'             centered=FALSE, profiled_if_noncenter=TRUE, scale="norm", diag=1.7)
#' get_elts_ab(h_hp_dx$hdx, h_hp_dx$hpdx, x, a=0.5, b=0.7, setting="ab_1/2_7/10",
#'             centered=FALSE, profiled_if_noncenter=FALSE, scale="norm", diag=1.9)
#' @export
get_elts_ab <- function(hdx, hpdx, x, a, b, setting,
                        centered=TRUE, profiled_if_noncenter=TRUE,
                        scale="", diagonal_multiplier=1){
  ### centered and profiled_if_noncenter IGNORED, just for compatibility
  if (b < 0 || 2*a <= b) {stop("a and b must satisfy 2a > b >= 0.")}
  n <- dim(x)[1]; p <- dim(x)[2]
  xa <- x^a
  if (a == 1/2){xa_1 <- 1/xa
  } else {xa_1 <- xa/x}
  if (b-a == 1) {xb_1 <- xa
  } else {xb_1 <- x^(b-1)}
  g_K <- c(crossprod(xa, hpdx*xa_1+(a-1)*hdx*xa_1/x)/n + diag(a*colMeans(hdx*xa_1*xa_1)))
  if (centered){
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hdx[,j]/n)*xa_1[,j]*xa)}))
    diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
  }
  #x2 <- x^2
  #Gamma_K <- Reduce(cbind, lapply(1:p, function(j){sweep(t(xa), 2, hdx[,j]*xa[,j]*xa[,j]/x2[,j]/n, "*") %*% xa}))
  #Gamma_K_eta <- crossprod(xa, hdx*xa*xb_1/x)/n
  #Gamma_eta <- colMeans(hdx*xb_1*xb_1)
  Gamma0 <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hdx[,j]/n)*cbind(-xa_1[,j]*xa, xb_1[,j]))}))
  Gamma_K <- Gamma0[-p-1, -c(1:p)*(p+1)]
  Gamma_K_eta <- Gamma0[-p-1, c(1:p)*(p+1)]
  Gamma_eta <- Gamma0[p+1, c(1:p)*(p+1)]
  remove(Gamma0)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  g_eta <- -colMeans(hpdx*xb_1+hdx*(b-1)*xb_1/x)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
}

#' The R implementation to get the elements necessary for calculations for the exponential square-root setting (a=0.5, b=0.5).
#'
#' @param hdx A matrix, \eqn{h(\mathbf{x})}{h(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param hpdx A matrix, \eqn{h'(\mathbf{x})}{h\'(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if centered=TRUE. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The setting \code{"exp"}.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' n <- 50
#' p <- 30
#' eta <- rep(0, p)
#' K <- diag(p)
#' domain <- make_domain("R+", p=p)
#' x <- gen(n, setting="exp", abs=FALSE, eta=eta, K=K, domain=domain, finite_infinity=100,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#' get_elts_exp(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_exp(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=FALSE, profiled_if_noncenter=TRUE,
#'       scale="norm", diag=1.7)
#' get_elts_exp(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=FALSE, profiled_if_noncenter=FALSE,
#'       scale="norm", diag=1.7)
#' @export
get_elts_exp <- function(hdx, hpdx, x, centered=TRUE, profiled_if_noncenter=TRUE, scale="", diagonal_multiplier=1){
  ## scale only used in the returned list
  ## Note that elts_exp with centered=TRUE is the same as elts_gamma with centered=TRUE
  n <- dim(x)[1]; p <- dim(x)[2]
  xsqrt <- sqrt(x)
  tmp <- (hpdx-0.5*hdx/x)/xsqrt
  g_K <- c(crossprod(xsqrt, tmp) / n + diag(colMeans(hdx/x)/2))
  if (centered){
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hdx[,j]/n)/xsqrt[,j]*xsqrt)}))
    diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="exp"))
  }
  xsqrt <- sqrt(x)
  Gamma0 <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hdx[,j]/n)/xsqrt[,j]*cbind(-xsqrt,1))}))
  Gamma_K <- Gamma0[-p-1, -c(1:p)*(p+1)]
  Gamma_K_eta <- Gamma0[-p-1, c(1:p)*(p+1)]
  Gamma_eta <- Gamma0[p+1, c(1:p)*(p+1)]
  remove(Gamma0)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  g_eta <- -colMeans(tmp)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="exp"))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="exp"))
}

#' The R implementation to get the elements necessary for calculations for the gamma setting (a=0.5, b=0).
#'
#' @param hdx A matrix, \eqn{h(\mathbf{x})}{h(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param hpdx A matrix, \eqn{h'(\mathbf{x})}{h\'(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if centered=TRUE. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The setting \code{"gamma"}.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' n <- 50
#' p <- 30
#' eta <- rep(0, p)
#' K <- diag(p)
#' domain <- make_domain("R+", p=p)
#' x <- gen(n, setting="gamma", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("min_pow", 1.5, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#' get_elts_gamma(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_gamma(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=FALSE, profiled_if_noncenter=TRUE,
#'        scale="norm", diag=1.7)
#' get_elts_gamma(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=FALSE, profiled_if_noncenter=FALSE,
#'        scale="norm", diag=1.9)
#' @export
get_elts_gamma <- function(hdx, hpdx, x, centered=TRUE, profiled_if_noncenter=TRUE, scale="", diagonal_multiplier=1){
  ## scale only used in the returned list
  ## Note that elts_exp with centered=TRUE is the same as elts_gamma with centered=TRUE
  n <- dim(x)[1]; p <- dim(x)[2]
  xsqrt <- sqrt(x)
  g_K <- c(crossprod(xsqrt, (hpdx-0.5*hdx/x)/xsqrt)/n + diag(0.5*colMeans(hdx/x)))
  if (centered){
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hdx[,j]/n)/xsqrt[,j]*xsqrt)}))
    diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gamma"))
  }
  Gamma0 <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hdx[,j]/n)*cbind(-xsqrt/xsqrt[,j], 1/x[,j]))}))
  Gamma_K <- Gamma0[-p-1, -c(1:p)*(p+1)]
  Gamma_K_eta <- Gamma0[-p-1, c(1:p)*(p+1)]
  Gamma_eta <- Gamma0[p+1, c(1:p)*(p+1)]
  remove(Gamma0)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  g_eta <- -colMeans((hpdx-hdx/x)/x)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gamma"))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gamma"))
}

#' The R implementation to get the elements necessary for calculations for the gaussian setting (a=1, b=1) on domains other than R^p.
#'
#' @param hdx A matrix, \eqn{h(\mathbf{x})}{h(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param hpdx A matrix, \eqn{h'(\mathbf{x})}{h\'(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if centered=TRUE. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The setting \code{"gaussian"}.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' n <- 50
#' p <- 30
#' mu <- rep(0, p)
#' K <- diag(p)
#' eta <- K %*% mu
#' domain <- make_domain("R+", p=p)
#' x <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#' get_elts_trun_gauss(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_trun_gauss(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=FALSE,
#'        profiled_if_noncenter=TRUE, scale="norm", diag=1.7)
#' get_elts_trun_gauss(h_hp_dx$hdx, h_hp_dx$hpdx, x, centered=FALSE,
#'        profiled_if_noncenter=FALSE, scale="norm", diag=1.9)
#' @export
get_elts_trun_gauss <- function(hdx, hpdx, x, centered=TRUE, profiled_if_noncenter=TRUE, scale="", diagonal_multiplier=1){
  n <- dim(x)[1]; p <- dim(x)[2]
  g_K <- crossprod(hpdx, x)/n
  Gamma_K <- do.call("cbind", lapply(1:p, function(k){1/n * crossprod(x, hdx[,k]*x)}))
  Gamma_eta <- colMeans(hdx)
  diag(g_K) <- diag(g_K) + Gamma_eta
  g_K <- c(g_K)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  if (centered){
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gaussian"))
  }
  Gamma_K_eta <- -crossprod(x, hdx)/n
  g_eta <- -colMeans(hpdx)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gaussian"))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gaussian"))
}

#' The R implementation to get the elements necessary for calculations for the gaussian setting on R^p.
#'
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if \code{centered==TRUE}. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The setting \code{"gaussian"}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}. Except for the \emph{profiled} setting, this is \eqn{\mathbf{xx}^{\top}/n}{xx'/n}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}. The minus column means of \code{x}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the\eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' n <- 50
#' p <- 30
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- mvtnorm::rmvnorm(n, mean=mu, sigma=solve(K))
#' # Equivalently:
#' \donttest{
#' x2 <- gen(n, setting="gaussian", abs=FALSE, eta=c(K%*%mu), K=K, domain=make_domain("R",p),
#'        finite_infinity=100, xinit=NULL, burn_in=1000, thinning=100)
#' }
#' get_elts_gauss(x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_gauss(x, centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#' @export
get_elts_gauss <- function(x, centered=TRUE, profiled_if_noncenter=TRUE, scale="", diagonal_multiplier=1){
  n <- dim(x)[1]; p <- dim(x)[2]
  Gamma_K <- crossprod(x)/n # p copies of this
  #g_K <- c(diag(p))
  diagonals_with_multiplier <- diag(Gamma_K)*diagonal_multiplier
  if (centered){
    return (list("n"=n, "p"=p, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gaussian"))
  }
  #Gamma_eta <- 1 # p copies of this
  Gamma_K_eta <- -colMeans(x) # p copies of this
  #g_eta <- rep(0,p)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gaussian"))
  tmp <- tcrossprod(Gamma_K_eta)
  Gamma_K <- Gamma_K - tmp
  diagonals_with_multiplier <- diagonals_with_multiplier - diag(tmp)
  #g_K does not change since g_eta = 0
  return (list("n"=n, "p"=p, "Gamma_K"=Gamma_K, "t1"=0, "t2"=Gamma_K_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gaussian"))
}

#' The R implementation to get the elements necessary for calculations for the log-log setting (a=0, b=0).
#'
#' @param hdx A matrix, \eqn{h(\mathbf{x})}{h(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param hpdx A matrix, \eqn{h'(\mathbf{x})}{h\'(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param setting A string, log_log.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if centered=TRUE. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The same setting as in the function argument.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' n <- 50
#' p <- 30
#' eta <- rep(0, p)
#' K <- diag(p)
#' domain <- make_domain("uniform", p=p, lefts=c(0), rights=c(1))
#' x <- gen(n, setting="log_log", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' h_hp <- get_h_hp("min_pow", 1.5, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#' get_elts_loglog(h_hp_dx$hdx, h_hp_dx$hpdx, x, setting="log_log", centered=TRUE,
#'        scale="", diag=1.5)
#' get_elts_loglog(h_hp_dx$hdx, h_hp_dx$hpdx, x, setting="log_log", centered=FALSE,
#'        profiled_if_noncenter=TRUE, scale="", diag=1.7)
#' get_elts_loglog(h_hp_dx$hdx, h_hp_dx$hpdx, x, setting="log_log", centered=FALSE,
#'        profiled_if_noncenter=FALSE, scale="", diag=1.9)
#' @export
get_elts_loglog <- function(hdx, hpdx, x, setting, centered=TRUE,
                            profiled_if_noncenter=TRUE,
                            scale="", diagonal_multiplier=1){
  ### centered and profiled_if_noncenter IGNORED, just for compatibility
  #offdiag_only <- FALSEE
  n <- dim(x)[1]; p <- dim(x)[2]
  if (any(x <= 0))
    stop("All entries in x must be positive.")
  logx <- log(x)
  h_over_xsq <- hdx / x^2
  hp_over_x <- hpdx / x
  g_K <- crossprod(logx, hp_over_x - h_over_xsq)/n + diag(colMeans(h_over_xsq))
  #if (offdiag_only)
  #  g_K <- t(t(g_K) - diag(g_K)) # g_j <- g_j - g_jj
  g_K <- c(g_K)
  if (centered) {
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){t(logx) %*% diag(h_over_xsq[,j]) %*% logx / n}))
  } else {
    logx_m1 <- cbind(logx, -1)
    Gamma0 <- Reduce(cbind, lapply(1:p, function(j){t(logx_m1) %*% diag(h_over_xsq[,j]) %*% logx_m1 / n}))
    Gamma_K <- Gamma0[-p-1, -c(1:p)*(p+1)]
  }
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p + rep(1:p,p)] * diagonal_multiplier
  if (centered)
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
  Gamma_K_eta <- Gamma0[-p-1, c(1:p)*(p+1)]
  #if (offdiag_only)
  #  Gamma_K_eta <- t(t(Gamma_K_eta) - diag(Gamma_K_eta)) # Gamma_K_eta_j <- Gamma_K_eta_j - Gamma_K_eta_jj
  Gamma_eta <- Gamma0[p+1, c(1:p)*(p+1)]
  remove(Gamma0)
  g_eta <- colMeans(-hp_over_x + h_over_xsq)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
}

#' The R implementation to get the elements necessary for calculations for the log-log setting (a=0, b=0) on the p-simplex.
#'
#' @param hdx A matrix, \eqn{h(\mathbf{x})}{h(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param hpdx A matrix, \eqn{h'(\mathbf{x})}{h\'(x)} applied to the distance of x from the boundary of the domain, should be of the same dimension as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param setting A string, log_log or log_log_sum0. If log_log_sum0, assumes that the true K has row and column sums 0 (see the A^d model), so only the off-diagonal entries will be estimated; the diagonal entries will be profiled out in the loss), so elements corresponding to the diagonals of K will be set to 0, and the loss will be rewritten in the off-diagonal entries only.
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if centered=TRUE. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. Returned without being checked or used in the function body. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The same setting as in the function argument.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' n <- 50
#' p <- 30
#' eta <- rep(0, p)
#' K <- -cov_cons("band", p=p, spars=3, eig=1)
#' diag(K) <- diag(K) - rowSums(K) # So that rowSums(K) == colSums(K) == 0
#' eigen(K)$val[(p-1):p] # Make sure K has one 0 and p-1 positive eigenvalues
#' domain <- make_domain("simplex", p=p)
#' x <- gen(n, setting="log_log_sum0", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' h_hp <- get_h_hp("min_pow", 2, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#'
#' elts_simplex_0 <- get_elts_loglog_simplex(h_hp_dx$hdx, h_hp_dx$hpdx, x,
#'        setting="log_log", centered=FALSE, profiled=FALSE, scale="", diag=1.5)
#'
#' # If want K to have row sums and column sums equal to 0; estimate off-diagonals only
#' elts_simplex_1 <- get_elts_loglog_simplex(h_hp_dx$hdx, h_hp_dx$hpdx, x,
#'        setting="log_log_sum0", centered=FALSE, profiled=FALSE, scale="", diag=1.5)
#' # All entries corresponding to the diagonals of K should be 0:
#' max(abs(sapply(1:p, function(j){c(elts_simplex_1$Gamma_K[j, (j-1)*p+1:p],
#'        elts_simplex_1$Gamma_K[, (j-1)*p+j])})))
#' max(abs(diag(elts_simplex_1$Gamma_K_eta)))
#' max(abs(diag(matrix(elts_simplex_1$g_K, nrow=p))))
#' @export
get_elts_loglog_simplex <- function(hdx, hpdx, x, setting,
                                    centered=TRUE,
                                    profiled_if_noncenter=TRUE,
                                    scale="", diagonal_multiplier=1){
  ### centered and profiled_if_noncenter IGNORED, just for compatibility
  n <- dim(x)[1]; p <- dim(x)[2]
  if (any(dim(hdx) != c(n, p-1)) || any(dim(hpdx) != c(n, p-1)))
    stop("hdx and hpdx must have dimension n x (p-1), i.e. xp must be excluded.")
  if (any(x <= 0))
    stop("All entries in x must be positive.")
  if (setting != "log_log" && setting != "log_log_sum0")
    stop("Setting must be log_log or log_log_sum0.")
  sum_to_zero <- (setting == "log_log_sum0")
  logx <- log(x)
  h_over_xsq_nop <- hdx / x[,-p]^2
  minus_h_over_x_xp_nop <- -hdx / x[,-p] / x[,p]
  sum_h_over_xmsq <- rowSums(hdx) / x[,p]^2
  hp_over_x_nop <- hpdx / x[,-p]
  sum_hp_over_xm <- rowSums(hpdx) / x[,p]

  g_K <- crossprod(logx, cbind(hp_over_x_nop - h_over_xsq_nop, -sum_h_over_xmsq - sum_hp_over_xm))/n
  g_K[-p,-p] <- g_K[-p,-p] + diag(colMeans(h_over_xsq_nop))
  tmp <- colMeans(minus_h_over_x_xp_nop)
  g_K[p,-p] <- g_K[p,-p] + tmp
  g_K[-p,p] <- g_K[-p,p] + tmp
  g_K[p,p] <- g_K[p,p] + mean(sum_h_over_xmsq)

  if (centered) {
    Gamma_K <- Reduce(cbind, lapply(1:(p-1), function(j){t(logx) %*% diag(h_over_xsq_nop[,j]) %*% logx / n}))
    Gamma_K <- cbind(Gamma_K, t(logx) %*% diag(sum_h_over_xmsq) %*% logx / n)
    Gamma_K_jp <- Reduce(cbind, lapply(1:(p-1), function(j){t(logx) %*% diag(minus_h_over_x_xp_nop[,j]) %*% logx / n}))
  } else {
    logx_m1 <- cbind(logx, -1)
    Gamma_K0 <- Reduce(cbind, lapply(1:(p-1), function(j){t(logx_m1) %*% diag(h_over_xsq_nop[,j]) %*% logx_m1 / n}))
    Gamma_K <- cbind(Gamma_K0[-p-1, -c(1:(p-1))*(p+1)]) # p * p*(p-1), interaction matrix for K_j (j!=p)
    Gamma_K_eta <- cbind(Gamma_K0[-p-1, c(1:(p-1))*(p+1)]) # p * (p-1), interactions between K_j and eta_j (j!=p)
    Gamma_eta <- c(Gamma_K0[p+1, c(1:(p-1))*(p+1)]) # p-1, coefficient on eta_j^2 (j!=p)
    remove(Gamma_K0)

    Gamma_K0_jp <- Reduce(cbind, lapply(1:(p-1), function(j){t(logx_m1) %*% diag(minus_h_over_x_xp_nop[,j]) %*% logx_m1 / n}))
    Gamma_K_jp <- Gamma_K0_jp[-p-1, -c(1:(p-1))*(p+1)]  # p * p(p-1), interactions between K_j (j!=p) and K_p, (p-1) blocks
    Gamma_K_eta_jp <- Gamma_K0_jp[-p-1, c(1:(p-1))*(p+1)] # p * (p-1), interactions between Kj (j!=p) and etap, or Kp and etaj (j!=p)
    Gamma_eta_jp <- Gamma_K0_jp[p+1, c(1:(p-1))*(p+1)] # p-1, interaction between eta_j (j!=p) and eta_p
    remove(Gamma_K0_jp)

    Gamma_Kp0 <- t(logx_m1) %*% diag(sum_h_over_xmsq) %*% logx_m1 / n
    Gamma_K <- cbind(Gamma_K, Gamma_Kp0[1:p, 1:p]) # p * p^2, interaction matrix for K_j (including j=p)
    Gamma_K_eta <- cbind(Gamma_K_eta, Gamma_Kp0[1:p, p+1]) # p * p, interactions between K_j and eta_j (including j=p)
    Gamma_eta <- c(Gamma_eta, Gamma_Kp0[p+1, p+1]) # p, coefficient on eta_j^2 (including j=p)
    remove(Gamma_Kp0)

    g_eta <- c(colMeans(h_over_xsq_nop - hp_over_x_nop),
               mean(sum_h_over_xmsq + sum_hp_over_xm))
  }

  if (sum_to_zero) {
    g_K <- t(t(g_K) - diag(g_K)) # g_j <- g_j - g_jj
    # Subtract j-th column from all columns of mat
    eliminate_col <- function(mat, j) {return (mat - mat[,j])}
    # Subtract j-th row from all rows of mat
    eliminate_row <- function(mat, j) {return (t(eliminate_col(t(mat), j)))}
    # Subtract j-th row from all rows of mat and then k-th column from all columns
    eliminate_row_col <- function(mat, j, k) {return (eliminate_row(eliminate_col(mat, k), j))}
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){
      eliminate_row_col(Gamma_K[,(j-1)*p+1:p], j, j)})) # Each block is interaction between Kj and Kj, so reduce row j and col j
    Gamma_K_jp <- Reduce(cbind, lapply(1:(p-1), function(j){
      eliminate_row_col(Gamma_K_jp[,(j-1)*p+1:p], j, p) # Each block is interaction betweeen Kj and Kp, so reduce row j and col p (symmetric)
    }))
  }
  g_K <- c(g_K)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p + rep(1:p,p)] * diagonal_multiplier
  if (centered)
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "Gamma_K_jp"=Gamma_K_jp,
                 "centered"=TRUE, "scale"=scale, "sum_to_zero"=sum_to_zero,
                 "diagonal_multiplier"=diagonal_multiplier,
                 "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
  if (sum_to_zero) {
    Gamma_K_eta <- t(t(Gamma_K_eta) - diag(Gamma_K_eta)) # sapply(1:p, function(j){Gamma_K_eta[,j]-Gamma_K_eta[j,j]}) # Each column is interaction between Kj and etaj, so reduce row j
    Gamma_Kj_etap <- sapply(1:(p-1), function(j){Gamma_K_eta_jp[,j]-Gamma_K_eta_jp[j,j]}) # Each column is interaction between Kj and etap, so reduce row j
    Gamma_Kp_etaj <- eliminate_row(Gamma_K_eta_jp, p) # Each column is interaction between Kp and etaj, so reduce row p
  } else {
    Gamma_Kj_etap <- Gamma_Kp_etaj <- Gamma_K_eta_jp
  }
  remove(Gamma_K_eta_jp)
  if (profiled_if_noncenter)
    stop("Profiled estimator not available for log_log on simplex domains.")
  return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta,
               "Gamma_K_jp"=Gamma_K_jp, "Gamma_Kj_etap"=Gamma_Kj_etap, "Gamma_Kp_etaj"=Gamma_Kp_etaj,
               "Gamma_eta"=Gamma_eta, "Gamma_eta_jp"=Gamma_eta_jp, "centered"=FALSE,
               "scale"=scale, "profiled_if_noncenter"=FALSE, "sum_to_zero"=sum_to_zero,
               "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier,
               setting=setting))
}

#' The function wrapper to get the elements necessary for calculations for all settings.
#'
#' @param h_hp A function that returns a list containing \code{hx=h(x)} (element-wise) and \code{hpx=hp(x)} (element-wise derivative of \eqn{h}) when applied to a vector or a matrix \code{x}, both of which has the same shape as \code{x}.
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param setting A string that indicates the distribution type, must be one of \code{"exp"}, \code{"gamma"}, \code{"gaussian"}, \code{"log_log"}, \code{"log_log_sum0"}, or of the form \code{"ab_NUM1_NUM2"}, where \code{NUM1} is the \code{a} value and \code{NUM2} is the \code{b} value, and \code{NUM1} and \code{NUM2} must be integers or two integers separated by "/", e.g. "ab_2_2", "ab_2_5/4" or "ab_2/3_1/2". If \code{domain$type == "simplex"}, only \code{"log_log"} and \code{"log_log_sum0"} are supported, and on the other hand \code{"log_log_sum0"} is supported for \code{domain$type == "simplex"} only.
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @param centered A boolean, whether in the centered setting(assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if \code{centered=TRUE}. Default to \code{TRUE}. Can only be \code{FALSE} if \code{setting == "log_log_sum0" && centered == FALSE}.
#' @param scale A string indicating the scaling method. If contains \code{"sd"}, columns are scaled by standard deviation; if contains \code{"norm"}, columns are scaled by l2 norm; if contains \code{"center"} and \code{setting == "gaussian" && domain$type == "R"}, columns are centered to have mean zero. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @param use_C Optional. A boolean, use C (\code{TRUE}) or R (\code{FALSE}) functions for computation. Default to \code{TRUE}. Ignored if \code{setting == "gaussian" && domain$type == "R"}.
#' @param tol Optional. A positive number. If \code{setting != "gaussian" || domain$type != "R"}, function stops if any entry if smaller than -tol, and all entries between -tol and 0 are set to tol, for numerical stability and to avoid violating the assumption that \eqn{h(\mathbf{x})>0}{h(x)>0} almost surely.
#' @param unif_dist Optional, defaults to \code{NULL}. If not \code{NULL}, \code{h_hp} must be \code{NULL} and \code{unif_dist(x)} must return a list containing \code{"g0"} of length \code{nrow(x)} and \code{"g0d"} of dimension \code{dim(x)}, representing the l2 distance and the gradient of the l2 distance to the boundary: the true l2 distance function to the boundary is used for all coordinates in place of h_of_dist; see "Estimating Density Models with Complex Truncation Boundaries" by Liu et al, 2019. That is, \eqn{(h_j\circ \phi_j)(x_i)}{(h_j\circ phi_j)(xi)} in the score-matching loss is replaced by \eqn{g_0(x_i)}{g0(xi)}, the l2 distance of xi to the boundary of the domain.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{domain_type}{The domain type. Same as domain$type in the input.}
#'   \item{setting}{The setting. Same as input.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}. A \eqn{p^2}-vector. Not returned if \code{setting == "gaussian" && domain$type == "R"} since it is just \eqn{diag(p)}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}. A vector of length \eqn{p^2} if \code{setting == "gaussian" && domain$type == "R"} or \eqn{p^3} otherwise.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}. A \eqn{p}-vector. Not returned if \code{setting == "gaussian" && domain$type == "R"} since it is just \eqn{numeric(p)}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}. If \code{setting == "gaussian" && domain$type == "R"}, returns a vector of length \eqn{p}, or \eqn{p^2} otherwise.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}. A \eqn{p}-vector. Not returned if \code{setting == "gaussian" && domain$type == "R"} since it is just \code{rep(1,p)}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' If \code{domain$type == "simplex", the following are also returned.}
#'   \item{Gamma_K_jp}{A matrix of size \code{p} by \code{p(p-1)}. The \code{(j-1)*p+1} through \code{j*p} columns represent the interaction matrix between the \code{j}-th column and the \code{m}-th column of \code{K}.}
#'   \item{Gamma_Kj_etap}{Non-centered only. A matrix of size \code{p} by \code{p(p-1)}. The \code{j}-th column represents the interaction between the \code{j}-th column of \code{K} and \code{eta[p]}.}
#'   \item{Gamma_Kp_etaj}{Non-centered only. A matrix of size \code{p} by \code{p(p-1)}. The \code{j}-th column represents the interaction between the \code{p}-th column of \code{K} and \code{eta[j]}. Note that it is equal to \code{Gamma_Kj_etap} if \code{setting} does not contain substring \code{"sum0"}.}
#'   \item{Gamma_eta_jp}{Non-centered only. A vector of size \code{p-1}. The \code{j}-th component represents the interaction between \code{eta[j]} and \code{eta[p]}.}
#' @details
#' Computes the \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix and the \eqn{\boldsymbol{g}}{g} vector for generalized score matching.
#'
#' Here, \eqn{\boldsymbol{\Gamma}}{\Gamma} is block-diagonal, and in the non-profiled non-centered setting, the \eqn{j}-th block is composed of \eqn{\boldsymbol{\Gamma}_{\mathbf{KK},j}}{\Gamma_{KK,j}}, \eqn{\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}}{\Gamma_{K\eta,j}} and its transpose, and finally \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta},j}}{\Gamma_{\eta\eta,j}}. In the centered case, only \eqn{\boldsymbol{\Gamma}_{\mathbf{KK},j}}{\Gamma_{KK,j}} is computed. In the profiled non-centered case, \deqn{\boldsymbol{\Gamma}_{j}\equiv\boldsymbol{\Gamma}_{\mathbf{KK},j}-\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}\boldsymbol{\Gamma}_{\boldsymbol{\eta}\boldsymbol{\eta},j}^{-1}\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta}}^{\top}.}{\Gamma_j=\Gamma_{KK,j}-\Gamma_{K\eta,j}\Gamma_{\eta\eta,j}^(-1)\Gamma_{K\eta}'.}
#'
#' Similarly, in the non-profiled non-centered setting, \eqn{\boldsymbol{g}}{g} can be partitioned \eqn{p} parts, each with a \eqn{p}-vector \eqn{\boldsymbol{g}_{\mathbf{K},j}}{g_{K,j}} and a scalar \eqn{g_{\boldsymbol{\eta},j}}{g_{\eta,j}}. In the centered setting, only \eqn{\boldsymbol{g}_{\mathbf{K},j}}{g_{K,j}} is needed. In the profiled non-centered case, \deqn{\boldsymbol{g}_j\equiv\boldsymbol{g}_{\mathbf{K},j}-\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta},j}^{-1}g_{\boldsymbol{\eta},j}.}{g_j=g_{K,j}-\Gamma_{K\eta,j}\Gamma_{\eta\eta,j}^(-1)g_{\eta,j}.}
#'
#' The formulae for the pieces above are
#' \deqn{\boldsymbol{\Gamma}_{\mathbf{KK},j}\equiv\frac{1}{n}\sum_{i=1}^nh\left(X_j^{(i)}\right){X_j^{(i)}}^{2a-2}{\boldsymbol{X}^{(i)}}^a{{\boldsymbol{X}^{(i)}}^a}^{\top},}{\Gamma_{KK,j}=1/n*\sum_{i=1}^n h(Xij)*Xij^(2a-2)*Xi^a*(Xi^a)',}
#' \deqn{\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta},j}\equiv-\frac{1}{n}\sum_{i=1}^nh\left(X_j^{(i)}\right){X_j^{(i)}}^{a+b-2}{\boldsymbol{X}^{(i)}}^a,}{\Gamma_{K\eta,j}=-1/n*\sum_{i=1}^n h(Xij)*Xij^(a+b-2)*Xi^a,}
#' \deqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta},j}\equiv\frac{1}{n}\sum_{i=1}^nh\left(X_j^{(i)}\right){X_j^{(i)}}^{2b-2},}{\Gamma_{\eta\eta,j}=1/n*\sum_{i=1}^n h(Xij)*Xij^(2b-2),}
#' \deqn{\boldsymbol{g}_{\mathbf{K},j}\equiv\frac{1}{n}\sum_{i=1}^n\left(h'\left(X_j^{(i)}\right){X_j^{(i)}}^{a-1}+(a-1)h\left(X_j^{(i)}\right){X_j^{(i)}}^{a-2}\right){\boldsymbol{X}^{(i)}}^a+ah\left(X_j^{(i)}\right){X_j^{(i)}}^{2a-2}\boldsymbol{e}_{j,p},}{g_{K,j}=1/n*\sum_{i=1}^n (h'(Xij)*Xij^(a-1)+(a-1)*h(Xij)*Xij^(a-2))*Xi^a+a*h(Xij)*Xij^(2a-2)*e_{j,p},}
#' \deqn{\boldsymbol{g}_{\boldsymbol{\eta},j}\equiv\frac{1}{n}\sum_{i=1}^n-h'\left(X_j^{(i)}\right){X_j^{(i)}}^{b-1}-(b-1)h\left(X_j^{(i)}\right){X_j^{(i)}}^{b-2},}{g_{\eta,j}=1/n*\sum_{i=1}^n -h'(Xij)*Xij^(b-1)-(b-1)*h(Xij)*Xij^(b-2)),}
#' where \eqn{\boldsymbol{e}_{j,p}}{e_{j,p}} is the \eqn{p}-vector with 1 at the \eqn{j}-th position and 0 elsewhere.
#'
#' In the profiled non-centered setting, the function also returns \eqn{t_1}{t1} and \eqn{t_2}{t2} defined as
#' \deqn{\boldsymbol{t}_1\equiv\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta}}^{-1}\boldsymbol{g}_{\boldsymbol{\eta}},\quad\boldsymbol{t}_2\equiv\boldsymbol{\Gamma}_{\boldsymbol{\eta\eta}}^{-1}\boldsymbol{\Gamma}_{\mathbf{K}\boldsymbol{\eta}}^{\top},}{t1=\Gamma_{\eta\eta}^(-1)g_{\eta}, t2=\Gamma_{\eta\eta}^(-1)\Gamma_{K\eta}',}
#' so that \eqn{\hat{\boldsymbol{\eta}}=\boldsymbol{t}_1-\boldsymbol{t}_2\mathrm{vec}(\hat{\mathbf{K}}).}{\hat{\eta}=t1-t2*vec(\hat{K}).}
#' @examples
#' n <- 30
#' p <- 10
#' eta <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#'
#' # Gaussian on R^p:
#' domain <- make_domain("R", p=p)
#' x <- mvtnorm::rmvnorm(n, mean=solve(K, eta), sigma=solve(K))
#' # Equivalently:
#' \donttest{
#' x2 <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, burn_in=1000, thinning=100)
#' }
#' get_elts(NULL, x, "gaussian", domain, centered=TRUE, scale="norm", diag=dm)
#' get_elts(NULL, x, "gaussian", domain, FALSE, profiled=FALSE, scale="sd", diag=dm)
#'
#' # Gaussian on R_+^p:
#' domain <- make_domain("R+", p=p)
#' x <- tmvtnorm::rtmvnorm(n, mean = solve(K, eta), sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' # Equivalently:
#' \donttest{
#' x2 <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, burn_in=1000, thinning=100)
#' }
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' get_elts(h_hp, x, "gaussian", domain, centered=TRUE, scale="norm", diag=dm)
#'
#' # Gaussian on sum(x^2) > 1 && sum(x^(1/3)) > 1 with x allowed to be negative
#' domain <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list("expression"="sum(x^2)>1", abs=FALSE, nonnegative=FALSE),
#'                       list("expression"="sum(x^(1/3))>1", abs=FALSE, nonnegative=FALSE)))
#' xinit <- rep(sqrt(2/p), p)
#' x <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' get_elts(h_hp, x, "gaussian", domain, centered=FALSE,
#'        profiled_if_noncenter=TRUE, scale="", diag=dm)
#'
#' # exp on ([0, 1] v [2,3])^p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#' x <- gen(n, setting="exp", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' h_hp <- get_h_hp("min_pow", 1.5, 3)
#' get_elts(h_hp, x, "exp", domain, centered=TRUE, scale="", diag=dm)
#' get_elts(h_hp, x, "exp", domain, centered=FALSE,
#'        profiled_if_noncenter=FALSE, scale="", diag=dm)
#'
#' # gamma on {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=TRUE)))
#' set.seed(1)
#' xinit <- c(1.5, 0.5, abs(stats::rnorm(p-2))+log(1.3))
#' x <- gen(n, setting="gamma", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("min_pow", 1.5, 3)
#' get_elts(h_hp, x, "gamma", domain, centered=TRUE, scale="", diag=dm)
#' get_elts(h_hp, x, "gamma", domain, centered=FALSE,
#'        profiled_if_noncenter=FALSE, scale="", diag=dm)
#'
#' # a0.6_b0.7 on {x in R_+^p: sum(log(x))<2 || (x1^(2/3)-1.3x2^(-3)<1 && exp(x1)+2.3*x2>2)}
#' domain <- make_domain("polynomial", p=p, rule="1 || (2 && 3)",
#'        ineqs=list(list("expression"="sum(log(x))<2", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x1^(2/3)-1.3x2^(-3)<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x1)+2.3*x2^2>2", abs=FALSE, nonnegative=TRUE)))
#' xinit <- rep(1, p)
#' x <- gen(n, setting="ab_3/5_7/10", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#' h_hp <- get_h_hp("min_pow", 1.4, 3)
#' get_elts(h_hp, x, "ab_3/5_7/10", domain, centered=TRUE, scale="", diag=dm)
#' get_elts(h_hp, x, "ab_3/5_7/10", domain, centered=FALSE,
#'        profiled_if_noncenter=TRUE, scale="", diag=dm)
#'
#' # log_log model on {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' x <- gen(n, setting="log_log", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100,
#'        verbose=TRUE)
#' h_hp <- get_h_hp("min_pow", 2, 3)
#' get_elts(h_hp, x, "log_log", domain, centered=TRUE, scale="", diag=dm)
#' get_elts(h_hp, x, "log_log", domain, centered=FALSE,
#'        profiled_if_noncenter=FALSE, scale="", diag=dm)
#' # Example of using the uniform distance function to boundary as in Liu (2019)
#' g0 <- function(x) {
#'        row_min <- apply(x, 1, min)
#'        row_which_min <- apply(x, 1, which.min)
#'        dist_to_sum_boundary <- apply(x, 1, function(xx){
#'                    (1 - sum(1:p * xx)) / sqrt(p*(p+1)*(2*p+1)/6)})
#'        grad_sum_boundary <- -(1:p) / sqrt(p*(p+1)*(2*p+1)/6)
#'        g0 <- pmin(row_min, dist_to_sum_boundary)
#'        g0d <- t(sapply(1:nrow(x), function(i){
#'           if (row_min[i] < dist_to_sum_boundary[i]){
#'              tmp <- numeric(ncol(x)); tmp[row_which_min[i]] <- 1
#'           } else {tmp <- grad_sum_boundary}
#'           tmp
#'        }))
#'        list("g0"=g0, "g0d"=g0d)
#' }
#' get_elts(NULL, x, "exp", domain, centered=TRUE, profiled_if_noncenter=FALSE,
#'        scale="", diag=dm, unif_dist=g0)
#'
#' # log_log_sum0 model on the simplex with K having row and column sums 0 (Aitchison model)
#' domain <- make_domain("simplex", p=p)
#' K <- -cov_cons("band", p=p, spars=3, eig=1)
#' diag(K) <- diag(K) - rowSums(K) # So that rowSums(K) == colSums(K) == 0
#' eigen(K)$val[(p-1):p] # Make sure K has one 0 and p-1 positive eigenvalues
#' x <- gen(n, setting="log_log_sum0", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' h_hp <- get_h_hp("min_pow", 2, 3)
#' h_hp_dx <- h_of_dist(h_hp, x, domain) # h and h' applied to distance from x to boundary
#'
#' # Does not assume K has 0 row and column sums
#' elts_simplex_0 <- get_elts(h_hp, x, "log_log", domain, centered=TRUE, profiled=FALSE,
#'        scale="", diag=1.5)
#'
#' # If want K to have row sums and column sums equal to 0 (Aitchison); estimate off-diagonals only
#' elts_simplex_1 <- get_elts(h_hp, x, "log_log_sum0", domain, centered=FALSE,
#'        profiled=FALSE, scale="", diag=1.5)
#' # All entries corresponding to the diagonals of K should be 0:
#' max(abs(sapply(1:p, function(j){c(elts_simplex_1$Gamma_K[j, (j-1)*p+1:p],
#'        elts_simplex_1$Gamma_K[, (j-1)*p+j])})))
#' max(abs(diag(elts_simplex_1$Gamma_K_eta)))
#' max(abs(diag(matrix(elts_simplex_1$g_K, nrow=p))))
#' @export
#' @useDynLib genscore, .registration = TRUE
get_elts <- function(h_hp, x, setting, domain, centered=TRUE, profiled_if_noncenter=TRUE,
                     scale="", diagonal_multiplier=1, use_C=TRUE, tol=.Machine$double.eps^0.5,
                     unif_dist=NULL){
  ## Note that in the result the diagonals of elts$Gamma_K are without multipliers.
  ## The diagonal entries with multipliers are stored in elts$diagonals_with_multiplier
  if (!(setting %in% c("exp", "gamma", "gaussian", "log_log", "log_log_sum0") || startsWith(setting, "ab_"))){
    stop("\"setting\" parameter must be one of exp, gamma, gaussian, log_log, log_log_sum0, or ab_A_B (e.g. ab_2_2, ab_2_5/4 or ab_2/3_1/2).")
  }
  if (!"checked" %in% names(domain))
    stop("domain must be an object returned by make_domain().")
  if (domain$type == "simplex") {
    if (!grepl("log_log", setting))
      stop("Currently only log_log and log_log_sum0 are supported for simplex domains.")
  } else if (setting == "log_log_sum0")
    stop("log_log_sum0 only supported for simplex domains.")

  n <- dim(x)[1]; p <- dim(x)[2]
  if (tol <= 0) {stop("tol must be >= 0.")}

  ### Violates the assumption that h(x)>0 almost surely since each column will have at least one 0
  #if (substr(scale, 1,3) == "min")
  #  x <- t(t(x)-apply(x, 2, min))
  if (grepl("center", scale)){
    if (setting == "gaussian" && domain$type == "R") {x <- scale(x, center=TRUE, scale=FALSE)
    } else {warning("\"center\" in scale ignored for non-gaussian settings.")}
  }
  if (grepl("sd", scale) && grepl("norm", scale)){
    stop("scale cannot contain both sd and norm.")
  }
  if (grepl("sd", scale)){
    sdx <- apply(x, 2, stats::sd)
    x <- sweep(x, MARGIN=2, sdx, `/`)
  } else if (grepl("norm", scale)){
    x <- scale(x, center=FALSE)
  }
  if (setting != "gaussian" || domain$type != "R"){
    if (is.null(unif_dist)) {
      h_hp_dx <- h_of_dist(h_hp, x, domain)
      hdx <- h_hp_dx$hdx # h(dist to boundary)
      hpdx <- h_hp_dx$hpdx # h'(dist to boundary(x)) = h'(dist to boundary) * dist'(x)
    } else {
      if (!is.null(h_hp))
        stop("Exactly one of h_hp and unif_dist should be NULL.")
      tmp <- unif_dist(x)
      if (length(tmp$g0) != nrow(x))
        stop("g0 returned by unif_dist(x) must have length equal to nrow(x).")
      if ((nrow(tmp$g0d) != nrow(x) || ncol(tmp$g0d) != domain$p_deemed))
        stop("g0d returned by unif_dist(x) must have nrow(x) rows (", nrow(x), ") and domain$p_deemed (", domain$p_deemed, ") columns.")
      hdx <- matrix(tmp$g0, nrow=nrow(x), ncol=ncol(x))
      hpdx <- tmp$g0d
    }
    if (any(is.infinite(hdx)) || any(is.infinite(hpdx)))
      stop("Infinite values encountered in h(get_dist(x, domain)) and/or hp(get_dist(x, domain)). Please use a bounded h function with one of the following modes (min_pow is recommended): 1, mcp, min_asinh, min_cosh, min_exp, min_log_pow, min_pow, min_sinh, min_softplus, tanh, truncated_sin, truncated_tan.")
    if (any(is.nan(hdx)) || any(is.nan(hpdx)))
      stop("NaN values encountered in h(get_dist(x, domain)) and/or hp(get_dist(x, domain)). If you are using a custom h/hp, please make sure it allows 0 and Inf as inputs. You may set h(Inf) to Inf but this is only allowed if all coordinates are bounded.")
  }
  if (setting == "exp"){
    if (use_C){
      if (centered) {
        res <- .C("elts_exp_c", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_exp_np", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      } else {
        res <- .C("elts_exp_p", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      }
    } else {return (c(get_elts_exp(hdx, hpdx, x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))}
  } else if (setting == "gamma") {
    if (use_C){
      if (centered) {
        res <- .C("elts_exp_c", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_gamma_np", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      } else {
        res <- .C("elts_gamma_p", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      }
    } else {return (c(get_elts_gamma(hdx, hpdx, x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))}
  } else if (setting == "gaussian" && domain$type != "R") {
    if (use_C){
      if (centered) {
        res <- .C("elts_gauss_c", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_gauss_np", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      } else {
        res <- .C("elts_gauss_p", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, "setting"=setting, "domain_type"=domain$type))
      }
    } else {return (c(get_elts_trun_gauss(hdx, hpdx, x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))}
  } else if (setting == "gaussian" && domain$type == "R") {
    #use_C ignored
    return (c(get_elts_gauss(x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))
  } else if (setting == "log_log" && domain$type != "simplex") {
    if (use_C){
      if (centered){
        res <- .C("elts_loglog_c", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x),
                  g1=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)),
                  diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), logx=as.double(numeric(n*p)),
                  h_over_xsq=as.double(numeric(n*p)), hp_over_x=as.double(numeric(n*p)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting, "domain_type"=domain$type))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_loglog_np", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x),
                  g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting, "domain_type"=domain$type))
      } else {
        res <- .C("elts_loglog_p", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x),
                  g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting, "domain_type"=domain$type))
      }
    } else {return (c(get_elts_loglog(hdx, hpdx, x, setting, centered=centered,
                                      profiled_if_noncenter=profiled_if_noncenter,
                                      scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))}
  } else if (domain$type == "simplex" && (setting == "log_log" || setting == "log_log_sum0")) {
    sum_to_zero <- (setting == "log_log_sum0")
    if (use_C){
      if (centered){
        res <- .C("elts_loglog_simplex_c", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx),
                  hpdx=as.double(hpdx), x=as.double(x), sum_to_zero=as.integer(sum_to_zero), g_K=as.double(numeric(p^2)),
                  Gamma_K=as.double(numeric(p^3)), Gamma_K_jp=as.double(numeric(p^2*(p-1))),
                  Gamma_eta=as.double(numeric(p-1)), Gamma_eta_jp=as.double(numeric(p-1)),
                  diagonal_multiplier=as.double(diagonal_multiplier),
                  diagonals_with_multiplier=as.double(numeric(p^2)),
                  logx=as.double(numeric(n*p)), h_over_xsq_nop=as.double(numeric(n*(p-1))),
                  minus_h_over_x_xp_nop=as.double(numeric(n*(p-1))), sum_h_over_xmsq=as.double(numeric(n)),
                  hp_over_x_nop=as.double(numeric(n*(p-1))), sum_hp_over_xm=as.double(numeric(n)),
                  mean_sum_h_over_xmsq=as.double(0), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g_K, "Gamma_K"=matrix(res$Gamma_K,nrow=p,ncol=p^2),
                     "Gamma_K_jp"=matrix(res$Gamma_K_jp,nrow=p,ncol=p*(p-1)),
                     "centered"=TRUE, "scale"=scale, "sum_to_zero"=sum_to_zero,
                     "diagonal_multiplier"=diagonal_multiplier,
                     "diagonals_with_multiplier"=res$diagonals_with_multiplier,
                     "setting"=setting, "domain_type"=domain$type))
      } else {
        if (profiled_if_noncenter)
          stop("Profiled estimator not available for log_log on simplex domains.")
        res <- .C("elts_loglog_simplex_np", nIn=as.integer(n), pIn=as.integer(p), hdx=as.double(hdx),
                  hpdx=as.double(hpdx), x=as.double(x), sum_to_zero=as.integer(sum_to_zero), g_K=as.double(numeric(p^2)),
                  g_eta=as.double(numeric(p)), Gamma_K=as.double(numeric(p^3)),
                  Gamma_K_eta=as.double(numeric(p^2)),
                  Gamma_K_jp=as.double(numeric(p^2*(p-1))),
                  Gamma_Kj_etap=as.double(numeric(p*(p-1))),
                  Gamma_Kp_etaj=as.double(numeric(p*(p-1))),
                  Gamma_eta=as.double(numeric(p)), Gamma_eta_jp=as.double(numeric(p-1)),
                  diagonal_multiplier=as.double(diagonal_multiplier),
                  diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g_K, "g_eta"=res$g_eta,
                     "Gamma_K"=matrix(res$Gamma_K,nrow=p,ncol=p^2),
                     "Gamma_K_eta"=matrix(res$Gamma_K_eta,nrow=p,ncol=p),
                     "Gamma_K_jp"=matrix(res$Gamma_K_jp,nrow=p,ncol=p*(p-1)),
                     "Gamma_Kj_etap"=matrix(res$Gamma_Kj_etap,nrow=p,ncol=p-1),
                     "Gamma_Kp_etaj"=matrix(res$Gamma_Kp_etaj,nrow=p,ncol=p-1),
                     "Gamma_eta"=res$Gamma_eta, "Gamma_eta_jp"=res$Gamma_eta_jp,
                     "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "sum_to_zero"=sum_to_zero,
                     "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier,
                     "setting"=setting, "domain_type"=domain$type))
      }
    } else {
      if (centered == FALSE && profiled_if_noncenter)
        stop("Profiled estimator not available for log_log on simplex domains.")
      return (c(get_elts_loglog_simplex(hdx, hpdx, x, setting,
                                        centered=centered, profiled_if_noncenter=FALSE,
                                        scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))
    }
  } else { # Only left option: ab models
    if (!startsWith(setting, "ab_"))
      stop("\"setting\" parameter must be one of exp, gamma, gaussian, log_log, log_log_sum0, or ab_A_B (e.g. ab_2_2, ab_2_5/4 or ab_2/3_1/2).")
    tmp <- parse_ab(setting)
    a <- tmp$a_numer / tmp$a_denom; b <- tmp$b_numer / tmp$b_denom #### !
    if (use_C){
      if (centered) {
        res <- .C("elts_ab_c", nIn=as.integer(n), pIn=as.integer(p), a=as.double(a), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting, "domain_type"=domain$type))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_ab_np", nIn=as.integer(n), pIn=as.integer(p), a=as.double(a), b=as.double(b), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting, "domain_type"=domain$type))
      } else {
        res <- .C("elts_ab_p", nIn=as.integer(n), pIn=as.integer(p), a=as.double(a), b=as.double(b), hdx=as.double(hdx), hpdx=as.double(hpdx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting, "domain_type"=domain$type))
      }
    } else {return (c(get_elts_ab(hdx, hpdx, x, a, b, setting, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier), domain_type=domain$type))}
  }
}


#' Random data generator from general \code{a}-\code{b} distributions with general domain types, assuming a and b are rational numbers.
#'
#' Random data generator from general \code{a}-\code{b} graphs with general domain types using adaptive rejection metropolis sampling (ARMS). x^(0/0) treated as log(x) and x^(n/0) as exp(x) for n non-zero. Density only guaranteed to be a proper density when 2*a > b >= 0 or when a = b = 0.
#'
#' @param n An integer, number of observations.
#' @param setting A string that indicates the distribution type, must be one of \code{"exp"}, \code{"gamma"}, \code{"gaussian"}, \code{"log_log"}, \code{"log_log_sum0"}, or of the form \code{"ab_NUM1_NUM2"}, where \code{NUM1} is the \code{a} value and \code{NUM2} is the \code{b} value, and \code{NUM1} and \code{NUM2} must be integers or two integers separated by "/", e.g. "ab_2_2", "ab_2_5/4" or "ab_2/3_1/2".
#' @param abs A boolean. If TRUE, density is rewritten as f(|x|), i.e. with |x|^(a_numer/a_denom) and |x|^(b_numer/b_denom)
#' @param eta A vector, the linear part in the distribution.
#' @param K A square matrix, the interaction matrix. There should exist some C > 0 such that \deqn{{\boldsymbol{x}^a}^{\top}\mathbf{K}{\boldsymbol{x}}^a/({\boldsymbol{x}^a}^{\top}{\boldsymbol{x}}^a) >= C} for all x in the domain (i.e. \code{K} is positive definite if \code{domain$type == "R"} and \code{K} is co-positive if \code{domain$type == "R+"}.). If \code{a_numer == a_denom == b_numer == b_denom == 0 && domain$type == "simplex"}, K can also have all row and column sums equal to 0 but have all but one eigenvalues (0) positive.
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @param finite_infinity A finite positive number. \code{Inf} in actual generation will be truncated to \code{finite_infinity} if applicable. Although the code will adaptively increase \code{finite_infinity}, the user should set it to a large number initially so that \code{abs(x) > finite_infinity} with very small probability.
#' @param xinit Optional. A \code{p}-vector, an initial point in the domain. If the domain is defined by more than one ineq or by one ineq containing negative coefficients, this must be provided. In the unlikely case where the function fails to automatically generate an initial point this should also be provided.
#' @param seed Optional. A number, the seed for the random generator.
#' @param burn_in Optional. A positive integer, the number of burn-in samples in ARMS to be discarded, meaning that samples from the first \code{burn_in} x \code{thinning} iterations will be discarded.
#' @param thinning Optional. A positive integer, thinning factor in ARMS. Samples are taken at iteration steps \eqn{(\mathrm{burn\_in}+1)\times\mathrm{thinning}}{(burn_in+1) x thinning}, ..., \eqn{(\mathrm{burn\_in}+n)\times\mathrm{thinning}}{(burn_in+n) x thinning}. Default to \code{100}.
#' @param verbose Optional. A boolean. If \code{TRUE}, prints a progress bar showing the progress. Defaults to \code{TRUE}.
#' @param remove_outofbound Optional. A logical, defaults to \code{TRUE}. If \code{TRUE}, a test whether each sample lies inside the domain will be done, which may take a while for larger sample sizes, and rows that do not lie in the domain will be removed (may happen for \code{domain$type == "polynomial"} with more than 1 ineq and an OR ("|") in \code{domain$rule}.).
#' @return An \eqn{n\times p}{n*p} matrix of samples, where \eqn{p} is the length of \code{eta}.
#' @details
#' NOTE: For polynomial domains with many ineqs and a rule containing "OR" ("|"), not all samples generated are guaranteed to be inside the domain. It is thus recommended to set \code{remove_outofbound} to \code{TRUE} and rerun the function with new initial points until the desired number of in-bound samples have been generated.
#'
#' Randomly generates \code{n} samples from the \code{p}-variate \code{a}-\code{b} distributions with parameters \eqn{\boldsymbol{\eta}}{\eta} and \eqn{\mathbf{K}}{K}, where \code{p} is the length of \eqn{\boldsymbol{\eta}}{\eta} or the dimension of the square matrix \eqn{\mathbf{K}}{K}.
#'
#' Letting \code{a=a_numer/a_denom} and \code{b=b_numer/b_denom}, the \code{a}-\code{b} distribution is proportional to
#' \deqn{\exp\left(-\frac{1}{2a}{\boldsymbol{x}^a}^{\top}\mathbf{K}{\boldsymbol{x}}^a+\boldsymbol{\eta}^{\top}\frac{\boldsymbol{x}^b-\mathbf{1}_p}{b}\right)}{exp(-x^a'Kx^a/(2a)+eta'(x^b-rep(1,p))/b)}.
#' Note that \code{x^(0/0)} is understood as \code{log(x)}, and \code{x^(n/0)} with nonzero \code{n} is \code{exp(n*x)}, and in both cases the \code{a} and \code{b} in the denominators in the density are treated as 1.
#'
#' @examples
#' n <- 20
#' p <- 10
#' eta <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#'
#' # Gaussian on sum(x^2) > 10 && sum(x^(1/3)) > 10 with x allowed to be negative
#' domain <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list("expression"="sum(x^2)>10", abs=FALSE, nonnegative=FALSE),
#'                       list("expression"="sum(x^(1/3))>10", abs=FALSE, nonnegative=FALSE)))
#' xinit <- rep(sqrt(20/p), p)
#' x <- gen(n, setting="gaussian", abs=FALSE, eta=eta,  K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#'
#' # exp on ([0, 1] v [2,3])^p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#' x <- gen(n, setting="exp", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#'
#' # gamma on {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
#' set.seed(1)
#' xinit <- c(1.5, 0.5, abs(stats::rnorm(p-2))+log(1.3))
#' x <- gen(n, setting="gamma", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#'
#' # a0.6_b0.7 on {x in R_+^p: sum(log(x))<2 || (x1^(2/3)-1.3x2^(-3)<1 && exp(x1)+2.3*x2>2)}
#' domain <- make_domain("polynomial", p=p, rule="1 || (2 && 3)",
#'        ineqs=list(list("expression"="sum(log(x))<2", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x1^(2/3)-1.3x2^(-3)<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x1)+2.3*x2^2>2", abs=FALSE, nonnegative=TRUE)))
#' xinit <- rep(1, p)
#' x <- gen(n, setting="ab_3/5_7/10", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=1e4, xinit=xinit, seed=2, burn_in=1000, thinning=100)
#'
#' # log_log model exp(-log(x) %*% K %*% log(x)/2 + eta %*% log(x)) on {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' x <- gen(n, setting="log_log", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#'
#' # log_log model on the simplex with K having row and column sums 0 (Aitchison model)
#' domain <- make_domain("simplex", p=p)
#' K <- -cov_cons("band", p=p, spars=3, eig=1)
#' diag(K) <- diag(K) - rowSums(K) # So that rowSums(K) == colSums(K) == 0
#' eigen(K)$val[(p-1):p] # Make sure K has one 0 and p-1 positive eigenvalues
#' x <- gen(n, setting="log_log_sum0", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#'
#' # Gumbel_Gumbel model exp(-exp(2x) %*% K %*% exp(2x)/2 + eta %*% exp(-3x)) on {sum(|x|) < 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"="sum(x)<1", abs=TRUE, nonnegative=FALSE)))
#' K <- diag(p)
#' x <- gen(n, setting="ab_2/0_-3/0", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' @export
#' @useDynLib genscore, .registration = TRUE
gen <- function(n, setting, abs, eta, K, domain, finite_infinity=NULL,
                xinit=NULL, seed=NULL, burn_in=1000, thinning=100,
                verbose=TRUE, remove_outofbound=TRUE) {
  if (!"checked" %in% names(domain))
    stop("domain must be an object returned by make_domain().")
  if (!domain$type %in% c("R", "R+", "uniform", "polynomial", "simplex"))
    stop("Domain type '", domain$type, "' not supported for data generation.")
  if (length(eta) != domain$p)
    stop("eta must have length ", domain$p, ".")
  if (nrow(K) != domain$p || ncol(K) != domain$p)
    stop("K must have dimension ", domain$p, "x", domain$p, ".")
  if (domain$type != "simplex" && any(diag(K) <= 0))
    stop("All diagonal entries of K must be positive except for simplex domains.")
  if (!is.null(seed)){set.seed(seed)
  } else {seed <- -1}

  if (setting == "exp") {
    a_numer <- b_numer <- 1; a_denom <- b_denom <- 2
  } else if (setting == "gamma") {
    a_numer <- 1; a_denom <- 2; b_numer <- b_denom <- 0
  } else if (setting == "gaussian") {
    a_numer <- a_denom <- b_numer <- b_denom <- 1
  } else if (startsWith(setting, "log_log")) {
    a_numer <- a_denom <- b_numer <- b_denom <- 0
  } else if (startsWith(setting, "ab_")) {
    tmp <- parse_ab(setting)
    a_numer <- tmp$a_numer; a_denom <- tmp$a_denom
    b_numer <- tmp$b_numer; b_denom <- tmp$b_denom
  } else
    stop("\"setting\" parameter must be one of exp, gamma, gaussian, log_log, log_log_sum0, or ab_A_B e.g. ab_2_2, ab_2_5/4 or ab_2/3_1/2.")
  if (a_numer == 1 && a_denom == 1 && b_numer == 1 && b_denom == 1)
    setting <- "gaussian"
  if ((a_numer == 0 && a_denom != 0) || (b_numer == 0 && b_denom != 0))
    stop("For a and b in the density, if the numerator is 0 then the denominator must also be 0 (-> log), otherwise the density is not defined (x^(0/n) for n!=0 always evaluates to 1 and we do not allow this; set eta to 0 instead).")

  if (a_denom != 1 || b_denom != 1) {
    if (domain$type == "R" || (domain$type == "uniform" && domain$lefts[1] < 0))
      stop("For settings with non-integer a and/or b, the domain must be non-negative.")
  }

  if (setting == "gaussian" && domain$type == "R") {
    Sigma <- solve(K)
    return (mvtnorm::rmvnorm(n, mean=c(Sigma%*%eta), sigma=Sigma))
  } else if (setting == "gaussian" && domain$type == "R+") {
    Sigma <- solve(K)
    return (tmvtnorm::rtmvnorm(n, mean=c(Sigma%*%eta), sigma=Sigma,
                               lower=rep(0, domain$p),
                               upper=rep(Inf, domain$p), algorithm = "gibbs",
                               burn.in.samples = burn_in, thinning = thinning))
  }

  if (domain$type == "uniform") {
    if (domain$left_inf || domain$right_inf) {
      if (is.null(finite_infinity))
        stop("You must specify finite_infinity for a uniform domain if any of the intervals are unbounded")
      finite_infinity <- update_finite_infinity_for_uniform(domain$lefts, domain$rights, finite_infinity)
    } else
      finite_infinity <- 10 * max(abs(domain$lefts[1]), abs(domain$rights[length(domain$rights)]))
  } else if (domain$type == "simplex") {
    finite_infinity <- 10
  } else {
    if (is.null(finite_infinity))
      stop("You must specify finite_infinity for R, R+, and polynomial-type domains.")
  }

  if (is.null(xinit)) {
    if (domain$type == "R")
      xinit <- stats::rnorm(domain$p)
    else if (domain$type == "R+")
      xinit <- abs(stats::rnorm(domain$p))
    else if (domain$type == "uniform")
      xinit <- random_init_uniform(domain$p, domain$lefts, domain$rights)
    else if (domain$type == "polynomial")
      xinit <- random_init_polynomial(domain)
    else if (domain$type == "simplex")
      xinit <- random_init_simplex(domain$p, tol=domain$simplex_tol)
  } else {
    if (domain$type == "simplex") {
      if (length(xinit) == domain$p_deemed)
        xinit <- c(xinit, 1 - sum(xinit))
      if (length(xinit) != domain$p)
        stop("For simplex, xinit must have length either domain$p_deemed = ", domain$p_deemed, " (1-sum(xinit) will be added as the last coordinate) OR domain$p = ", domain$p, ", in which case sum(xinit) must be 1.")
      if (abs(sum(xinit) - 1) > domain$simplex_tol)
        stop("sum(xinit) must be 1 (within a small tolerance) if xinit provided has length domain$p = ", domain$p, ".")
    } else if (length(xinit) != domain$p)
      stop("If specified, xinit must have length domain$p = ", domain$p, ".")
    if (!in_bound(xinit, domain))
      stop("xinit supplied is not in the domain.")
  }
  if (max(abs(xinit)) > finite_infinity / 10) {
    warning("Randomly generalized initial point has maximum absolute value ",
            max(abs(xinit)), "; using 10x of this as the finite infinity.")
    finite_infinity <- max(abs(xinit)) * 10
  }
  finite_infinity <- max(finite_infinity, 100)
  res <- do.call(".C", c(list("rab_arms",
                              nsamp=as.integer(n),
                              burnin=as.integer(burn_in),
                              p=as.integer(domain$p),
                              every=as.integer(thinning),
                              a_numer=as.integer(a_numer),
                              a_denom=as.integer(a_denom),
                              b_numer=as.integer(b_numer),
                              b_denom=as.integer(b_denom),
                              abs=as.integer(abs),
                              xinit=as.double(xinit),
                              xres=as.double(numeric(n*domain$p)),
                              eta=as.double(eta),
                              K=as.double(K),
                              #seed=as.integer(seed),
                              finite_infinity=as.double(finite_infinity)),
                         domain_for_C(domain),
                         list("verbose"=as.integer(verbose), "errno"=as.integer(0))
  ))
  if (res$errno)
    stop("Error occurred in C -> rab_arms() called in gen(). Check if the density is well-defined on your domain. (If any fractional powers are involved, i.e. a and/or b are non-integer, the domain must be non-negative.)")
  x <- t(matrix(res$xres, nrow=domain$p))
  if (remove_outofbound) {
    if (verbose)
      message("Testing if all samples are inside the domain.")
    inbound <- in_bound(x, domain)
    if (!all(inbound))
      warning("Row number(s) of generated samples not in domain: ", paste(which(!inbound), collapse=", "), ". Removed.")
    x <- x[which(inbound==1), ]
  }
  return (x)
}

#' Finds the distance of each element in a matrix x to the its boundary of the domain while fixing the others in the same row.
#'
#' Finds the distance of each element in a matrix \code{x} to its boundary of the \code{domain} while fixing the others in the same row.
#'
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @return A list that contains \code{h(dist(x, domain))} and \code{h\'(dist(x, domain))}.
#'   \item{dx}{Coordinate-wise distance to the boundary.}
#'   \item{dpx}{Coordinate-wise derivative of \code{dx}.}
#' @details
#' Returned matrix \code{dx} has its \code{i,j}-th component the distance of \eqn{x_{i,j}} to the boundary of \code{domain}, assuming \eqn{x_{i,-j}} are fixed. The matrix has the same size of \code{x} (\code{n} by \code{p}), or if \code{domain$type == "simplex"} and \code{x} has full dimension \code{p}, it has \code{p-1} columns.
#' Returned matrix \code{dpx} contains the component-wise derivatives of \code{dx} in its components. That is, its \code{i,j}-th component is 0 if \eqn{x_{i,j}} is unbounded or is bounded from both below and above or is at the boundary, or -1 if \eqn{x_{i,j}} is closer to its lower boundary (or if its bounded from below but unbounded from above), or 1 otherwise.
#' @examples
#' n <- 20
#' p <- 10
#' eta <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#'
#' # Gaussian on R^p:
#' domain <- make_domain("R", p=p)
#' x <- mvtnorm::rmvnorm(n, mean=solve(K, eta), sigma=solve(K))
#' # Equivalently:
#' \donttest{
#' x2 <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, burn_in=1000, thinning=100)
#' }
#' dist <- get_dist(x, domain)
#' # dx is all Inf and dpx is all 0 since each coordinate is unbounded with domain R
#' c(all(is.infinite(dist$dx)), all(dist$dpx==0))
#'
#' # exp on R_+^p:
#' domain <- make_domain("R+", p=p)
#' x <- gen(n, setting="exp", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' # dx is x and dpx is 1; with domain R+, the distance of x to the boundary is just x itself
#' c(max(abs(dist$dx - x))<.Machine$double.eps^0.5, all(dist$dpx == 1))
#'
#' # Gaussian on sum(x^2) > p with x allowed to be negative
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste("sum(x^2)>", p), abs=FALSE, nonnegative=FALSE)))
#' x <- gen(n, setting="gaussian", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' quota <- p - (rowSums(x^2) - x^2) # How much should xij^2 at least be so that sum(xi^2) > p?
#' # How far is xij from +/-sqrt(quota), if quota >= 0?
#' dist_to_bound <- abs(x[quota >= 0]) - abs(sqrt(quota[quota >= 0]))
#' max(abs(dist$dx[is.finite(dist$dx)] - dist_to_bound)) # Should be equal to our own calculations
#' # dist'(x) should be the same as the sign of x
#' all(dist$dpx[is.finite(dist$dx)] == sign(x[quota >= 0]))
#' # quota is negative <-> sum of x_{i,-j}^2 already > p <-> xij unbounded
#' #   given others <-> distance to boundary is Inf
#' all(quota[is.infinite(dist$dx)] < 0)
#'
#' # gamma on ([0, 1] v [2,3])^p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#' x <- gen(n, setting="gamma", abs=FALSE, eta=eta, K=K, domain=domain,
#'        xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' dist <- get_dist(x, domain)
#' # If 0 <= xij <= 1, distance to boundary is min(x-0, 1-x)
#' max(abs(dist$dx - pmin(x, 1-x))[x >= 0 & x <= 1])
#' # If 0 <= xij <= 1, dist'(xij) is 1 if it is closer to 0, or -1 if it is closer 1,
#' #   assuming xij %in% c(0, 0.5, 1) with probability 0
#' all((dist$dpx == 2 * (1-x > x) - 1)[x >= 0 & x <= 1])
#' # If 2 <= xij <= 3, distance to boundary is min(x-2, 3-x)
#' max(abs(dist$dx - pmin(x-2, 3-x))[x >= 2 & x <= 3])
#' # If 2 <= xij <= 3, dist'(xij) is 1 if it is closer to 2, or -1 if it is closer 3,
#' #   assuming xij %in% c(2, 2.5, 3) with probability 0
#' all((dist$dpx == 2 * (3-x > x-2) - 1)[x >= 2 & x <= 3])
#'
#' # a0.6_b0.7 on {x1 > 1 && 0 < x2 < 1 && x3 > 0 && ... && xp > 0}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
#' set.seed(1)
#' xinit <- c(1.5, 0.5, abs(stats::rnorm(p-2)) + log(1.3))
#' x <- gen(n, setting="ab_3/5_7/10", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=xinit, seed=2, burn_in=1000, thinning=100,
#'        verbose=TRUE)
#' dist <- get_dist(x, domain)
#' # x_{i1} has uniform bound [1, +Inf), so its distance to its boundary is x_{i1} - 1
#' max(abs(dist$dx[,1] - (x[,1] - 1)))
#' # x_{i2} has uniform bound [log(1.3), 1], so its distance to its boundary
#' #   is min(x_{i2} - log(1.3), 1 - x_{i2})
#' max(abs(dist$dx[,2] - pmin(x[,2] - log(1.3), 1 - x[,2])))
#' # x_{ij} for j >= 3 has uniform bound [log(1.3), +Inf), so its distance to its boundary
#' #   is simply x_{ij} - log(1.3)
#' max(abs(dist$dx[,3:p] - (x[,3:p] - log(1.3))))
#' # dist\'(xi2) is 1 if it is closer to log(1.3), or -1 if it is closer 1,
#' #   assuming x_{i2} %in% c(log(1.3), (1+log(1.3))/2, 1) with probability 0
#' all((dist$dpx[,2] == 2 * (1 - x[,2] > x[,2] - log(1.3)) - 1))
#' # x_{ij} for j != 2 is bounded from below but unbounded from above, so dist\'(xij) is always 1
#' all(dist$dpx[,-2] == 1)
#'
#' # log_log model on {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' x <- gen(n, setting="log_log", abs=FALSE, eta=eta, K=K, domain=domain,
#'        finite_infinity=100, xinit=NULL, seed=2, burn_in=1000, thinning=100)
#' dist <- get_dist(x, domain)
#' # Upper bound for j * xij so that sum_j j * xij <= 1
#' quota <- 1 - (rowSums(t(t(x) * 1:p)) - t(t(x) * 1:p))
#' # Distance of xij to its boundary is min(xij - 0, quota_{i,j} / j - xij)
#' max(abs(dist$dx - pmin((t(t(quota) / 1:p) - x), x)))
#'
#' # log_log model on the simplex with K having row and column sums 0 (Aitchison model)
#' domain <- make_domain("simplex", p=p)
#' K <- -cov_cons("band", p=p, spars=3, eig=1)
#' diag(K) <- diag(K) - rowSums(K) # So that rowSums(K) == colSums(K) == 0
#' eigen(K)$val[(p-1):p] # Make sure K has one 0 and p-1 positive eigenvalues
#' x <- gen(n, setting="log_log_sum0", abs=FALSE, eta=eta, K=K, domain=domain,
#'         xinit=NULL, seed=2, burn_in=1000, thinning=100, verbose=TRUE)
#' # Note that dist$dx and dist$dpx only has p-1 columns -- excluding the last coordinate in x
#' dist <- get_dist(x, domain)
#' # Upper bound for x_{i,j} so that x_{i,1} + ... + x_{i,p-1} <= 1
#' quota <- 1 - (rowSums(x[,-p]) - x[,-p])
#' # Distance of x_{i,j} to its boundary is min(xij - 0, quota_{i,j} - xij)
#' max(abs(dist$dx - pmin(quota - x[,-p], x[,-p])))
#' @export
#' @useDynLib genscore, .registration = TRUE
get_dist <- function(x, domain){
  # x must be a vector of length domain$p or a matrix of size n_sample x domain$p.
  # For simplex, the size can be p - 1 or p, where in the latter case only the first p-1 components are kept.
  if (!"checked" %in% names(domain))
    stop("domain must be an object returned by make_domain().")
  if (is.null(dim(x))) {
    if (length(x) != domain$p_deemed && length(x) != domain$p) # Assuming only simplex domains have different p_deemed and p!
      stop("If x is a dimensionless vector, its length must be domain$p = ", domain$p, " (or ", domain$p_deemed, " for simplex). Otherwise it must be a matrix with ", domain$p, " (or ", domain$p_deemed, " for simplex) columns.")
    x <- matrix(x, nrow=1)
  }

  if (domain$type == "simplex" && ncol(x) == domain$p)
    x <- x[, 1:domain$p_deemed, drop=FALSE]
  if (ncol(x) != domain$p_deemed)
    stop("x must have domain$p = ", domain$p, " (or ", domain$p_deemed, " for simplex) columns.")
  if (!"checked" %in% names(domain))
    stop("domain must be an object returned by make_domain().")

  inbound <- in_bound(x, domain)
  if (!all(inbound))
    stop("Row number(s) of x provided not in domain: ", paste(which(!inbound), collapse=", "))

  if (domain$type %in% c("uniform", "polynomial", "simplex")) {
    res <- do.call(".C",
                   c(list("dist", n=as.integer(nrow(x)),
                          p=as.integer(domain$p_deemed),
                          x=as.double(t(x)),
                          dists=numeric(length(x)),
                          dist_ps=integer(length(x))),
                     domain_for_C(domain),
                     list("errno"=as.integer(0))))
    if (res$errno)
      stop("Error occurred in C -> dist().")
    return (list("dx"=t(matrix(res$dists, domain$p_deemed, nrow(x))),
                 "dpx"=t(matrix(res$dist_ps, domain$p_deemed, nrow(x)))))
  } else if (domain$type == "R") {
    return (list("dx"=matrix(Inf, nrow=nrow(x), ncol=ncol(x)),
                 "dpx"=matrix(0, nrow=nrow(x), ncol=ncol(x))))
  } else if (domain$type == "R+") {
    return (list("dx"=x, "dpx"=matrix(1, nrow=nrow(x), ncol=ncol(x))))
  } else
    stop("Domain type ", domain$type, " not supported.")
}

#' Generator of h and hp (derivative of h) functions.
#'
#' Generator of \code{h} and \code{hp} (derivative of \eqn{h}) functions.
#'
#' @param mode A string, see details.
#' @param para May be optional. A number, the first parameter. Default to \code{NULL}.
#' @param para2 May be optional. A number, the second parameter. If \code{mode} is one of the adaptive mode below, this specifies the percentile (see details). Default to \code{NULL}.
#' @return A function that returns a list containing \code{hx=h(x)} (element-wise) and \code{hpx=hp(x)} (element-wise derivative of \eqn{h}) when applied to a vector (for mode names not ending with \code{"_ada"} only) or a matrix \code{x}, with both of the results having the same shape as \code{x}.
#' @details
#' The \code{mode} parameter can be chosen from the options listed below along with the corresponding definitions of \code{h} under appropriate choices of \code{para} and \code{para2} parameters. Unless otherwise noted, \code{para} and \code{para2}, must both be strictly positive if provided, and are set to 1 if not provided. Functions \code{h} and \code{hp} should only be applied to non-negative values \code{x} and this is not enforced or checked by the functions.
#' Internally calls \code{get_h_hp_vector}.
#' \describe{
#'     \item{\code{asinh}}{An asinh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\mathrm{asinh}(\mathrm{para}\cdot\boldsymbol{x})=\log\left(\mathrm{para}\cdot\boldsymbol{x}+\sqrt{(\mathrm{para}\cdot\boldsymbol{x})^2+1}\right)}{h(x)=asinh(para*x)=log(para*x+sqrt((para*x)^2+1))}. Unbounded and takes one parameter. Equivalent to \code{min_asinh(x, para, Inf)}.}
#'     \item{\code{cosh}}{A shifted cosh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\cosh(\mathrm{para}\cdot\boldsymbol{x})-1}{h(x)=cosh(para*x)-1}. Unbounded and takes one parameter. Equivalent to \code{min_cosh(x, para, Inf)}.}
#'     \item{\code{exp}}{A shifted exponential function \eqn{\boldsymbol{h}(\boldsymbol{x})=\exp(\mathrm{para}\cdot\boldsymbol{x})-1}{h(x)=exp(para*x)-1}. Unbounded and takes one parameter. Equivalent to \code{min_exp(x, para, Inf)}.}
#'     \item{\code{identity}}{The identity function \eqn{\boldsymbol{h}(\boldsymbol{x})=\boldsymbol{x}}{h(x)=x}. Unbounded and does not take any parameter. Equivalent to \code{pow(x, 1)} or \code{min_pow(x, 1, Inf)}.}
#'     \item{\code{log_pow}}{A power function on a log scale \eqn{\boldsymbol{h}(\boldsymbol{x})=\log(1+\boldsymbol{x})^{\mathrm{para}}}{log(1+x)^para}. Unbounded and takes one parameter. Equivalent to \code{min_log_pow(x, para, Inf)}.}
#'     \item{\code{mcp}}{Treating \eqn{\lambda}=para, \eqn{\gamma}=para2, the step-wise MCP function applied element-wise: \eqn{\lambda x-x^2/(2\gamma)}{\lambdax-x^2/(2\gamma)} if \eqn{x\leq\lambda\gamma}{x<=\lambda\gamma}, or \eqn{\gamma\lambda^2/2} otherwise. Bounded and takes two parameters.}
#'     \item{\code{min_asinh}}{A truncated asinh function applied element-wise: \eqn{\min(\mathrm{asinh}(\mathrm{para}\cdot\boldsymbol{x}),\mathrm{para}_2)}{pmin(asinh(para*x), para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_asinh_ada}}{Adaptive version of \code{min_asinh}.}
#'     \item{\code{min_cosh}}{A truncated shifted cosh function applied element-wise: \eqn{\min(\cosh(\mathrm{para}\cdot\boldsymbol{x})-1,\mathrm{para}_2)}{pmin(cosh(para*x)-1, para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_cosh_ada}}{Adaptive version of \code{min_cosh}.}
#'     \item{\code{min_exp}}{A truncated shifted exponential function applied element-wise: \eqn{\boldsymbol{h}(\boldsymbol{x})=\min(\exp(\mathrm{para}\cdot\boldsymbol{x})-1,\mathrm{para}_2)}{pmin(exp(para*x)-1, para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_exp_ada}}{Adaptive version of \code{min_exp}.}
#'     \item{\code{min_log_pow}}{A truncated power on a log scale applied element-wise: \eqn{\boldsymbol{h}(\boldsymbol{x})=\min(\log(1+\boldsymbol{x}),\mathrm{para}_2)^{\mathrm{para}}}{pmin(log(1+x), para2)^para}. Bounded and takes two parameters.}
#'     \item{\code{min_log_pow_ada}}{Adaptive version of \code{min_log_pow}.}
#'     \item{\code{min_pow}}{A truncated power function applied element-wise: \eqn{\boldsymbol{h}(\boldsymbol{x})=\min(\boldsymbol{x},\mathrm{para}_2)^{\mathrm{para}}}{pmin(x, para2)^para}. Bounded and takes two parameters.}
#'     \item{\code{min_pow_ada}}{Adaptive version of \code{min_pow}.}
#'     \item{\code{min_sinh}}{A truncated sinh function applied element-wise: \eqn{\min(\sinh(\mathrm{para}\cdot\boldsymbol{x}),\mathrm{para}_2)}{pmin(sinh(para*x), para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_sinh_ada}}{Adaptive version of \code{min_sinh}.}
#'     \item{\code{min_softplus}}{A truncated shifted softplus function applied element-wise: \eqn{\min(\log(1+\exp(\mathrm{para}\cdot\boldsymbol{x}))-\log(2),\mathrm{para}_2)}{pmin(log(1+exp(para*x))-log(2), para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_softplus_ada}}{Adaptive version of \code{min_softplus}.}
#'     \item{\code{pow}}{A power function \eqn{\boldsymbol{h}(\boldsymbol{x})=\boldsymbol{x}^{\mathrm{para}}}{h(x)=x^para}. Unbounded and takes two parameter. Equivalent to \code{min_pow(x, para, Inf)}.}
#'     \item{\code{scad}}{Treating \eqn{\lambda}=para, \eqn{\gamma}=para2, the step-wise SCAD function applied element-wise: \eqn{\lambda x}{\lambdax} if \eqn{x\leq\lambda}{x<=\lambda}, or \eqn{(2\gamma\lambda x-x^2-\lambda^2)/(2(\gamma-1))}{(2\gamma\lambdax-x^2-\lambda^2)/(2(\gamma-1))} if \eqn{\lambda<x<\gamma\lambda}, or \eqn{\lambda^2(\gamma+1)/2} otherwise. Bounded and takes two parameters, where \code{para2} must be larger than 1, and will be set to 2 by default if not provided.}
#'     \item{\code{sinh}}{A sinh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\sinh(\mathrm{para}\cdot\boldsymbol{x})}{h(x)=sinh(para*x)}. Unbounded and takes one parameter. Equivalent to \code{min_sinh(x, para, Inf)}.}
#'     \item{\code{softplus}}{A shifted softplus function \eqn{\boldsymbol{h}(\boldsymbol{x})=\log(1+\exp(\mathrm{para}\cdot\boldsymbol{x}))-\log(2)}{h(x)=log(1+exp(para*x))-log(2)}. Unbounded and takes one parameter. Equivalent to \code{min_softplus(x, para, Inf)}.}
#'     \item{\code{tanh}}{A tanh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\tanh(\mathrm{para}\cdot\boldsymbol{x})}{h(x)=tanh(para*x)}. Bounded and takes one parameter.}
#'     \item{\code{truncated_sin}}{A truncated sin function applied element-wise: \eqn{\sin(\mathrm{para}\cdot x)}{sin(para*x)} if \eqn{\mathrm{para}\cdot x\leq\pi/2}{para*x<=\pi/2}, or 1 otherwise. Bounded and takes one parameter.}
#'     \item{\code{truncated_tan}}{A truncated tan function applied element-wise: \eqn{\tan(\mathrm{para}\cdot x)}{tan(para*x)} if \eqn{\mathrm{para}\cdot x\leq\pi/4}{para*x<=\pi/4}, or 1 otherwise. Bounded and takes one parameter.}
#'  }
#' For the adaptive modes (names ending with \code{"_ada"}), \code{h} and \code{hp} are first applied to \code{x} without truncation. Then inside each column, values that are larger than the \code{para2}-th quantile will be truncated. The quantile is calculated using finite values only, and if no finite values exist the quantile is set to 1.
#' For example, if \code{mode == "min_pow_ada"}, \code{para == 2}, \code{para2 == 0.4}, the \code{j}-th column of the returned \code{hx} will be \code{pmin(x[,j]^2, stats::quantile(x[,j]^2, 0.4))}, and the \code{j}-th column of \code{hpx} will be \code{2*x[,j]*(x[,j] <= stats::quantile(x[,j]^2, 0.4))}.
#' @examples
#' get_h_hp("mcp", 2, 4)(0:10)
#' get_h_hp("min_log_pow", 1, log(1+3))(matrix(0:11, nrow=3))
#' get_h_hp("min_pow", 1.5, 3)(seq(0, 5, by=0.5))
#' get_h_hp("min_softplus")(matrix(seq(0, 2, by=0.1), nrow=7))
#'
#' get_h_hp("min_log_pow_ada", 1, 0.4)(matrix(0:49, nrow=10))
#' get_h_hp("min_pow_ada", 2, 0.3)(matrix(0:49, nrow=10))
#' get_h_hp("min_softplus_ada", 2, 0.6)(matrix(seq(0, 0.49, by=0.01), nrow=10))
#' @export
get_h_hp <- function(mode, para=NULL, para2=NULL){
  if (startsWith(mode, "min_") && endsWith(mode, "_ada"))
    return (get_h_hp_adaptive(gsub("_ada$", "", mode), para, para2))
  h_hp <- get_h_hp_vector(mode, para, para2)
  return (function(x) {
    if (is.matrix(x)) {
      hx_hpx <- t(apply(x, 1, h_hp))
      return (list(hx=hx_hpx[,1:ncol(x)], hpx=hx_hpx[,ncol(x)+1:ncol(x)]))
    } else {
      hx_hpx <- h_hp(x)
      return (list(hx=hx_hpx[,1], hpx=hx_hpx[,2]))
    }
  })
}

#' Generator of adaptive h and hp (derivative of h) functions.
#'
#' Generator of adaptive \code{h} and \code{hp} (derivative of \eqn{h}) functions.
#'
#' @param mode A string, the corresponding mode (with the suffix \code{"_ada"} removed from the input to \code{get_h_hp()}). Must be one of the modes starting with \code{"min_"} supported by \code{get_h_hp_vector()}.
#' @param para Must be provided, but can be \code{NULL}. A number, the first parameter; see \code{get_h_hp()} or \code{get_h_hp_vector()}.
#' @param percentile A number, the percentile for column-wise truncation on \code{hx} and \code{hpx}.
#' @return A function that returns a list containing \code{hx=h(x)} (element-wise) and \code{hpx=hp(x)} (element-wise derivative of \eqn{h}) when applied to a matrix \code{x}, with both of the results having the same shape as \code{x}.
#' @details
#' Helper function of \code{get_h_hp()}. Please refer to \code{get_hs_hp()}.
#' @examples
#' get_h_hp_adaptive("min_log_pow", 1, 0.4)(matrix(0:49, nrow=10))
#' get_h_hp_adaptive("min_pow", 2, 0.3)(matrix(0:49, nrow=10))
#' get_h_hp_adaptive("min_softplus", 2, 0.6)(matrix(seq(0, 0.49, by=0.01), nrow=10))
#'
#' hx_hpx <- get_h_hp_adaptive("min_log_pow", 1, 0.4)(matrix(0:49, nrow=10))
#' hx_hpx2 <- get_h_hp("min_log_pow_ada", 1, 0.4)(matrix(0:49, nrow=10))
#' c(max(abs(hx_hpx$hx - hx_hpx2$hx)), max(abs(hx_hpx$hpx - hx_hpx2$hpx)))
#' @export
get_h_hp_adaptive <- function(mode, para, percentile){
  if (percentile < 0 || percentile > 1) stop("para2 for adaptive modes must be between 0 and 1.")
  h_hp <- get_h_hp_vector(mode, para, Inf)
  return (function(x) {
    if (!is.matrix(x)) {stop("h_hp with adaptive modes can only be applied to matrices.")
    } else {
      hx_hpx <- t(apply(x, 1, h_hp))
      hx <- hx_hpx[,1:ncol(x)]; hpx <- hx_hpx[,ncol(x)+1:ncol(x)]
      for (j in 1:ncol(x)) {
        quant <- stats::quantile(hx[,j][is.finite(hx[,j])], percentile)
        if (is.na(quant)) quant <- 1
        truncated <- which(hx[,j] > quant)
        hx[truncated,j] <- quant
        hpx[truncated,j] <- 0
      }
      return (list(hx=hx, hpx=hpx))
    }
  })
}

#' Generator of h and hp (derivative of h) functions.
#'
#' Generator of \code{h} and \code{hp} (derivative of \eqn{h}) functions.
#'
#' @param mode A string, see details.
#' @param para May be optional. A number, the first parameter. Default to \code{NULL}.
#' @param para2 May be optional. A number, the second parameter. Default to \code{NULL}.
#' @return A function that returns a matrix with \code{hx=h(x)} (element-wise) and \code{hpx=hp(x)} (element-wise derivative of \eqn{h}) \code{cbind}ed when applied to a vector or a matrix \code{x}, where if \code{x} is a vector, the returned value will have two columns and number of rows equal to \code{length(x)}, otherwise it will have the same number of rows as \code{x} and number of columns doubled.
#' @details
#' Helper function of \code{get_h_hp()}. Please refer to \code{get_hs_hp()}.
#' @examples
#' get_h_hp_vector("mcp", 2, 4)
#' get_h_hp_vector("min_log_pow", 1, log(1+3))
#' get_h_hp_vector("min_pow", 1, 3)
#' get_h_hp_vector("min_softplus")
#' @export
get_h_hp_vector <- function(mode, para=NULL, para2=NULL){
  if (is.null(mode) || mode == "")
    stop ("Mode must be chosen from one of the following:",
          "1, asinh, cosh, exp, identity, log_pow, mcp, min_asinh, min_asinh_ada, ",
          "min_cosh, min_cosh_ada, min_exp, min_exp_ada, min_log_pow, min_log_pow_ada, ",
          "min_pow, min_pow_ada, min_sinh, min_sinh_ada, min_softplus, min_softplus_ada, ",
          "pow, scad, sinh, softplus, tanh, truncated_sin, truncated_tan.")
  number_of_params <- list("1"=0, "asinh"=1, "cosh"=1, "exp"=1, "identity"=0, "log_pow"=1, "mcp"=2,
                           "min_asinh"=2, "min_cosh"=2, "min_exp"=2, "min_log_pow"=2,
                           "min_pow"=2, "min_sinh"=2, "min_softplus"=2, "pow"=1,
                           "scad"=2, "sinh"=1, "softplus"=1, "tanh"=1, "truncated_sin"=1,
                           "truncated_tan"=1)
  if (!mode %in% names(number_of_params)){
    stop("Mode ", mode, " not supported.")
  }
  if (number_of_params[[mode]]){
    if (is.null(para)) {message("para not provided, default to 1.\n"); para <- 1
    } else if (para <= 0) {stop("para must be strictly positive.")}
  }
  if (number_of_params[[mode]] == 2){
    if (is.null(para2)) {
      if (mode == "scad") {message("para2 not provided, default to 2.\n"); para2 <- 2}
      else {message("para2 not provided, default to 1."); para2 <- 1}
    } else if (para2 <= 0) {stop("para2 must be strictly positive.")}
  }
  if (mode == "1")
    return (function(x){cbind(rep(1, length(x)), 0)})
  else if (mode == "asinh")
    return (function(x){cbind(asinh(para*x), para/sqrt((para*x)^2+1))})
  else if (mode == "cosh")
    return (function(x){cbind(cosh(para*x)-1, para*sinh(para*x))})
  else if (mode == "exp")
    return (function(x){tmp <- exp(para*x); cbind(tmp-1, para*tmp)})
  #else if (mode == "extpow") {if(para<=0)stop("para must be > 0"); return (list(h=function(a){abs(a)^para},hp=switch((para==0)+1, function(a){para*sign(a)*abs(a)^(para-1)}, function(a){rep(0,length(a))})))}
  else if (mode == "identity")
    return (function(x){cbind(x, 1)})
  else if (mode == "log_pow")
    return (function(x){tmp<-log(1+x); cbind(tmp^para, ifelse(is.finite(x), para*tmp^(para-1)/(1+x), 0))})
  else if (mode == "mcp") {
    lambda <- para; gamma <- para2
    return (function(x){cbind(ifelse(x<=gamma*lambda, lambda*x-x^2/2/gamma, gamma*lambda^2/2), pmax(lambda-x/gamma,0))}) # para,para2 > 0
  } else if (mode == "min_asinh")
    return (function(x){tmp<-asinh(para*x); cbind(pmin(tmp,para2), ifelse(tmp<=para2, para/sqrt((para*x)^2+1), 0))})
  else if (mode == "min_cosh")
    return (function(x){tmp<-cosh(para*x)-1; cbind(pmin(tmp,para2), ifelse(tmp<=para2, para*sinh(para*x), 0))})
  else if (mode == "min_exp")
    return (function(x){tmp<-exp(para*x)-1; cbind(pmin(tmp,para2), ifelse(tmp<=para2, para*exp(para*x), 0))})
  else if (mode == "min_log_pow")
    return (function(x){tmp<-log(1+x); cbind((pmin(tmp,para2))^para, ifelse(tmp<=para2, para*tmp^(para-1)/(1+x), 0))})
  else if (mode == "min_pow")
    return (function(x){cbind(pmin(x,para2)^para, ifelse(x<para2, para*x^(para-1), 0))})
  else if (mode == "min_sinh")
    return (function(x){tmp<-sinh(para*x); cbind(pmin(tmp,para2), ifelse(tmp<=para2, para*cosh(para*x), 0))})
  else if (mode == "min_softplus")
    return (function(x){tmp<-log(1+exp(para*x))-log(2); cbind(pmin(tmp,para2), ifelse(tmp<=para2, para/(1+exp(-para*x)), 0))})
  else if (mode == "pow")
    return (function(x){cbind(x^para, para*x^(para-1))})
  else if (mode == "scad") {
    if(para2 <= 1) stop("para2 must be > 1.")
    lambda <- para; gamma <- para2
    return (function(x){
      cbind(ifelse(x<=lambda, lambda*x, ifelse(x<gamma*lambda, (2*gamma*lambda*x-x^2-lambda^2)/(2*(gamma-1)), lambda^2*(gamma+1)/2)),
            ifelse(x<=lambda, lambda, ifelse(x<gamma*lambda, (gamma*lambda-x)/(gamma-1), 0)))})
  } else if (mode == "sinh")
    return (function(x){cbind(sinh(para*x), para*cosh(para*x))})
  else if (mode == "softplus")
    return (function(x){cbind(log(1+exp(para*x))-log(2), para/(1+exp(-para*x)))})
  else if (mode == "tanh")
    return (function(x){cbind(tanh(para*x), para/cosh(para*x)^2)}) # naturally bounded
  else if (mode == "truncated_sin")
    return (function(x){tmp<-para*x; res <- cbind(rep(1, length(x)), 0); idx <- tmp<=pi/2
    res[idx, ] <- cbind(sin(tmp[idx]), para*cos(tmp[idx])); res})
  else if (mode == "truncated_tan")
    return (function(x){tmp<-para*x; res <- cbind(rep(1, length(x)), 0); idx <- tmp<=pi/4
    res[idx, ] <- cbind(tan(tmp[idx]), para/cos(tmp[idx])^2); res})
  else
    stop("Mode not supported!")
}

#' Estimate \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta} using elts from \code{get_elts()} given one \eqn{\lambda_{\mathbf{K}}}{\lambda_K} (and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta} if non-profiled non-centered) and applying warm-start with strong screening rules.
#'
#' Estimate \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta} using elts from \code{get_elts()} given one \eqn{\lambda_{\mathbf{K}}}{\lambda_K} (and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta} if non-profiled non-centered) and applying warm-start with strong screening rules.
#'
#' @param elts A list, elements necessary for calculations returned by \code{get_elts()}.
#' @param symmetric A string. If equals \code{"symmetric"}, estimates the minimizer \eqn{\mathbf{K}}{K} over all symmetric matrices; if \code{"and"} or \code{"or"}, use the "and"/"or" rule to get the support.
#' @param lambda1 A number, the penalty parameter for \eqn{\mathbf{K}}{K}.
#' @param lambda2 A number, the penalty parameter for \eqn{\boldsymbol{\eta}}{\eta}. Default to \code{0}. Cannot be \code{Inf} if non-profiled non-centered.
#' @param tol Optional. A number, the tolerance parameter.
#' @param maxit Optional. A positive integer, the maximum number of iterations.
#' @param previous_res Optional. A list or \code{NULL}, the returned list by this function run previously with another lambda value.
#' @param is_refit A boolean, in the refit mode for BIC estimation if \code{TRUE}. If \code{TRUE}, \code{lambda1}, \code{previous_lambda1} and \code{lambda2} are all set to \code{0}, and estimation is restricted to entries in exclude that are \code{0}.
#' @return
#'     \item{converged}{A boolean indicating convergence.}
#'     \item{crit}{A number, the final penalized loss.}
#'     \item{edges}{A vector of the indices of entries in the \code{K} estimate that are non-zero.}
#'     \item{eta}{A p-vector, the \code{eta} estimate. Returned only if \code{elts$centered == FALSE}.}
#'     \item{eta_support}{A vector of the indices of entries in the \code{eta} estimate that are non-zero. Returned only if \code{elts$centered == FALSE && elts$profiled_if_noncenter == TRUE}.}
#'     \item{iters}{An integer, number of iterations run.}
#'     \item{K}{A p*p matrix, the \code{K} estimate.}
#'     \item{n}{An integer, the number of samples.}
#'     \item{p}{An integer, the dimension.}
#'     \item{is_refit,lambda1,maxit,previous_lambda1,symmetric,tol}{Same as in the input.}
#'     \item{lambda2}{Same as in the input, and returned only if \code{elts$centered == FALSE} and \cr
#'              \code{elts$profiled_if_noncenter == FALSE}.}
#' @details
#' If \code{elts$domain_type == "simplex"}, \code{symmetric != "symmetric"} or \code{elts$centered == FALSE && elts$profiled_if_noncenter} are currently not supported.
#' If \code{elts$domain_type == "simplex"} and \code{elts$setting} constains substring \code{"sum0"}, it is assumed that the column and row sums of \code{K} are all 0 and estimation will be done by profiling out the diagonal entries.
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()}, as the
#' #   way to call this function (\code{get_results()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=FALSE, scale="norm", diag=dm)
#' test_nc_np <- get_results(elts_gauss_np, symmetric="symmetric", lambda1=0.35,
#'                 lambda2=2, previous_res=NULL, is_refit=FALSE)
#' test_nc_np2 <- get_results(elts_gauss_np, symmetric="and", lambda1=0.25,
#'                  lambda2=2, previous_res=test_nc_np, is_refit=FALSE)
#'
#' elts_gauss_p <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                centered=FALSE, profiled=TRUE, scale="norm", diag=dm)
#' test_nc_p <- get_results(elts_gauss_p, symmetric="symmetric",
#'                lambda1=0.35, lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#'
#' elts_gauss_c <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                centered=TRUE, scale="norm", diag=dm)
#' test_c <- get_results(elts_gauss_c, symmetric="or", lambda1=0.35,
#'                lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#'
#' @export
#' @useDynLib genscore, .registration = TRUE
get_results <- function(elts, symmetric, lambda1, lambda2=0, tol=1e-6, maxit=10000,
                        previous_res=NULL, is_refit=FALSE){
  if (lambda1 < 0) stop("lambda1 must be non-negative.")
  if (!elts$centered && !elts$profiled_if_noncenter && lambda2 < 0) stop("lambda2 must be non-negative.")
  if (tol <= 0) stop("tol must be positive.")
  if (as.integer(maxit) < 1) stop("maxit must be a positive integer.")
  if (!symmetric %in% c("symmetric", "and", "or")){
    stop("Parameter symmetric must be one of \"symmetric\", \"and\", or \"or\".")
  }
  if (is_refit) {  ## If in refit mode, lambdas are set to 0.
    if (is.null(previous_res)) {stop("previous_res must be provided if is_refit == TRUE.")}  ## Otherwise do not know which edges to restrict to.
    if (symmetric != previous_res$symmetric){
      warning("symmetric changed to ", previous_res$symmetric, " according to previous_res$symmetric.")
      symmetric <- previous_res$symmetric
    }
    lambda1 <- lambda2 <- 0
  }
  if (is.null(previous_res)){ ## If previous result not given, fit from scratch.
    previous_res <- list()
    previous_res$lambda1 <- lambda1  # Set previous_lambda1 to current lambda1, and strong screening rules merely become KKT conditions
    previous_res$K <- diag(1,elts$p) ## Initialize K to the identity
    exclude <- matrix(0,elts$p,elts$p) ## Do not exclude any edge
    if ((!elts$centered) && (!elts$profiled_if_noncenter)){ ## If non-profiled non-centered
      previous_res$eta <- exclude_eta <- numeric(elts$p)  ## Initialize eta to the zero vector
    }
  } else {
    if (lambda1 != 0 && abs(previous_res$lambda1-lambda1)/lambda1 <= tol &&
        (elts$centered || elts$profiled_if_noncenter ||
         (is.infinite(lambda2) && is.infinite(previous_res$lambda2)) ||
         (!is.infinite(lambda2) && lambda2 != 0 && abs(previous_res$lambda2-lambda2)/lambda2 <= tol)) &&
        previous_res$symmetric == symmetric &&
        maxit == previous_res$maxit &&
        tol != 0 && abs(previous_res$tol-tol)/tol <= tol &&
        previous_res$is_refit == is_refit){ ## If same parameters as in previous_res, return that
      previous_res$iters <- 0
      return (previous_res)
    }
    if (is.null(previous_res$edges)) { ## If no edges given, do not exclude
      exclude <- matrix(0, elts$p, elts$p)
    } else { ## If edges given, exclude all non-edges
      exclude <- 1 - diag(elts$p)
      exclude[previous_res$edges] <- 0
      previous_res$K[exclude==1] <- 0
    }
    if ((!elts$centered) && (!elts$profiled_if_noncenter)){
      if (is.null(previous_res$eta_support)) { # If eta support not given, do not exclude
        exclude_eta <- numeric(elts$p)
      } else {
        exclude_eta <- rep(1, elts$p)
        exclude_eta[previous_res$eta_support] <- 0
        previous_res$eta[exclude_eta==1] <- 0
      }
    }
  }
  manual_ncnp_to_c <- !elts$centered && !elts$profiled_if_noncenter && is.infinite(lambda2)
  if (manual_ncnp_to_c) elts$centered <- TRUE
  if (elts$domain_type == "simplex") {
    sum_to_zero <- grepl("sum0", elts$setting)
    if (sum_to_zero)
      previous_res$K[diag(elts$p) == 1] <- 0
    if (symmetric != "symmetric") #### This should be updated in the future
      stop("Asymmetric estimation not supported for models on the simplex.")
    if (elts$centered) {
      test <- .C("simplex_centered", p=as.integer(elts$p), sum_to_zero=as.integer(sum_to_zero),
                 Gamma_K=as.double(elts$Gamma_K), Gamma_K_jp=as.double(elts$Gamma_K_jp),
                 g_K=as.double(elts$g_K), K=as.double(previous_res$K), lambda1=as.double(lambda1),
                 tol=as.double(tol), maxit=as.integer(maxit), iters=as.integer(0), converged=as.integer(0), crit=as.double(0),
                 exclude=as.integer(exclude), previous_lambda1=as.double(previous_res$lambda1),
                 is_refit=as.integer(is_refit), diagonals_with_multiplier=as.double(elts$diagonals_with_multiplier), PACKAGE="genscore")
      test$Gamma_K_jp <- NULL
    } else {
      if (elts$profiled_if_noncenter)
        stop("Profiled non-centered estimators not supported for models on the simplex.")
      test <- .C("simplex_full", p=as.integer(elts$p), sum_to_zero=as.integer(sum_to_zero),
                 Gamma_K=as.double(elts$Gamma_K), Gamma_K_eta=as.double(elts$Gamma_K_eta),
                 Gamma_K_jp=as.double(elts$Gamma_K_jp), Gamma_Kj_etap=as.double(elts$Gamma_Kj_etap),
                 Gamma_Kp_etaj=as.double(elts$Gamma_Kp_etaj), Gamma_eta=as.double(elts$Gamma_eta),
                 Gamma_eta_jp=as.double(elts$Gamma_eta_jp), g_K=as.double(elts$g_K),
                 g_eta=as.double(elts$g_eta), K=as.double(previous_res$K), eta=as.double(previous_res$eta),
                 lambda1=as.double(lambda1), lambda2=as.double(lambda2),
                 tol=as.double(tol), maxit=as.integer(maxit),
                 iters=as.integer(0), converged=as.integer(0), crit=as.double(0),
                 exclude=as.integer(exclude), exclude_eta=as.integer(exclude_eta),
                 previous_lambda1=as.double(previous_res$lambda1), is_refit=as.integer(is_refit),
                 diagonals_with_multiplier=as.double(elts$diagonals_with_multiplier))
      test$eta_support <- which(abs(test$eta) > tol)  ## For refit
      test$Gamma_K_eta <- test$Gamma_K_jp <- test$Gamma_Kj_etap <- test$Gamma_Kp_etaj <- test$Gamma_eta <- test$Gamma_eta_jp <- test$exclude_eta <- NULL
    }
    if (sum_to_zero) {
      test$K <- matrix(test$K, elts$p, elts$p)
      diag(test$K) <- diag(test$K) - rowSums(test$K) # Fill in the diagonals
    }
  } else {
    if (elts$centered || elts$profiled_if_noncenter){
      # Arguments are the same in the two cases. But cannot do call_name <- "..." and .C(call_name, ...) since R would complain.
      if (symmetric == "symmetric") {
        test <- .C("profiled", p = as.integer(elts$p), Gamma_K = as.double(elts$Gamma_K), g_K = as.double(elts$g_K),
                   K = as.double(previous_res$K), lambda1 = as.double(lambda1), tol = as.double(tol), maxit = as.integer(maxit),
                   iters = as.integer(0), converged = as.integer(0), crit = as.double(0), exclude = as.integer(exclude),
                   previous_lambda1 = as.double(previous_res$lambda1), is_refit = as.integer(is_refit),
                   diagonals_with_multiplier = as.double(elts$diagonals_with_multiplier),
                   gauss = as.integer(elts$setting=="gaussian" && elts$domain_type=="R"),
                   PACKAGE="genscore")
      } else {
        test <- .C("profiled_asymm", p = as.integer(elts$p), Gamma_K = as.double(elts$Gamma_K), g_K = as.double(elts$g_K),
                   K = as.double(previous_res$K), lambda1 = as.double(lambda1), tol = as.double(tol), maxit = as.integer(maxit),
                   iters = as.integer(0), converged = as.integer(0), crit = as.double(0), exclude = as.integer(exclude),
                   previous_lambda1 = as.double(previous_res$lambda1), is_refit = as.integer(is_refit),
                   diagonals_with_multiplier = as.double(elts$diagonals_with_multiplier),
                   gauss = as.integer(elts$setting=="gaussian" && elts$domain_type=="R"),
                   PACKAGE="genscore")
      }
    } else {
      # Arguments are the same in the two cases. But cannot do call_name <- "..." and .C(call_name, ...) since R would complain.
      if (symmetric == "symmetric") {
        test <- .C("full", p = as.integer(elts$p), Gamma_K = as.double(elts$Gamma_K), Gamma_K_eta = as.double(elts$Gamma_K_eta),
                   Gamma_eta = as.double(elts$Gamma_eta), g_K=as.double(elts$g_K), g_eta=as.double(elts$g_eta), K = as.double(previous_res$K),
                   eta=as.double(previous_res$eta), lambda1 = as.double(lambda1), lambda2 = as.double(lambda2), tol = as.double(tol),
                   maxit = as.integer(maxit), iters = as.integer(0), converged = as.integer(0), crit = as.double(0),
                   exclude = as.integer(exclude), exclude_eta = as.integer(exclude_eta),
                   previous_lambda1 = as.double(previous_res$lambda1), is_refit = as.integer(is_refit),
                   diagonals_with_multiplier = as.double(elts$diagonals_with_multiplier),
                   gauss = as.integer(elts$setting=="gaussian" && elts$domain_type=="R"), PACKAGE="genscore")
      } else {
        test <- .C("full_asymm", p = as.integer(elts$p), Gamma_K = as.double(elts$Gamma_K), Gamma_K_eta = as.double(elts$Gamma_K_eta),
                   Gamma_eta = as.double(elts$Gamma_eta), g_K=as.double(elts$g_K), g_eta=as.double(elts$g_eta), K = as.double(previous_res$K),
                   eta=as.double(previous_res$eta), lambda1 = as.double(lambda1), lambda2 = as.double(lambda2), tol = as.double(tol),
                   maxit = as.integer(maxit), iters = as.integer(0), converged = as.integer(0), crit = as.double(0),
                   exclude = as.integer(exclude), exclude_eta = as.integer(exclude_eta),
                   previous_lambda1 = as.double(previous_res$lambda1), is_refit = as.integer(is_refit),
                   diagonals_with_multiplier = as.double(elts$diagonals_with_multiplier),
                   gauss = as.integer(elts$setting=="gaussian" && elts$domain_type=="R"), PACKAGE="genscore")
      }
      test$eta_support <- which(abs(test$eta) > tol)  ## For refit
      test$Gamma_K_eta <- test$Gamma_eta <- test$exclude_eta <- NULL
    }
  }
  test$K <- matrix(test$K, elts$p, elts$p)
  if (manual_ncnp_to_c) {
    elts$centered <- FALSE
    test$eta <- numeric(elts$p)
    test$eta_support <- integer(0)
    test$lambda2 <- Inf
  }
  if ((!elts$centered) && elts$profiled_if_noncenter){
    if (elts$setting == "gaussian" && elts$domain_type == "R")
      test$eta <- elts$t1 - elts$t2 %*% test$K ## Note: elts$t1 = 0
    else
      test$eta <- elts$t1 - sapply(1:elts$p, function(k){crossprod(elts$t2[(k-1)*elts$p+1:elts$p], test$K[,k])})
  }
  if (symmetric == "symmetric")
    test$edges <- which(abs(test$K) > tol & diag(elts$p) == 0) ## all off-diagonal
  else if (symmetric == "and")
    test$edges <- which(abs(test$K) > tol & abs(t(test$K)) > tol & diag(elts$p) == 0)
  else ## "or"
    test$edges <- which((abs(test$K) > tol | abs(t(test$K)) > tol) & diag(elts$p) == 0)
  test$symmetric <- symmetric
  test$n <- elts$n
  test$diagonals_with_multiplier <- test$g_K <- test$Gamma_K <- test$g_eta <- test$exclude <- test$gauss <- NULL
  return (test)
}


#' Searches for a tight bound for \eqn{\lambda_{\boldsymbol{K}}}{\lambda_K} that gives the empty or complete graph starting from a given lambda with a given step size
#'
#' Searches for the smallest lambda that gives the empty graph (if \code{lower == FALSE}) or the largest that gives the complete graph (if \code{lower == TRUE}) starting from the given lambda, each time updating by multiplying or dividing by \code{step} depending on the search direction.
#'
#' @param elts A list, elements necessary for calculations returned by \code{get_elts()}.
#' @param symmetric A string. If equals \code{"symmetric"}, estimates the minimizer \eqn{\mathbf{K}}{K} over all symmetric matrices; if \code{"and"} or \code{"or"}, use the "and"/"or" rule to get the support
#' @param lambda A number, the initial searching point for \eqn{\lambda_{\mathbf{K}}}{\lambda_K}.
#' @param lambda_ratio A positive number (or \code{Inf}), the fixed ratio \eqn{\lambda_{\mathbf{K}}}{\lambda_K} and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta}, if \eqn{\lambda_{\boldsymbol{\eta}}\neq 0}{\lambda_\eta!=0} (non-profiled) in the non-centered setting.
#' @param step A number, the multiplicative constant applied to lambda at each iteration. Must be strictly larger than 1.
#' @param lower A boolean. If \code{TRUE}, finds the largest possible lambda that gives the complete graph (a \eqn{lower} bound). If \code{FALSE}, finds the smallest possible lambda that gives the empty graph (an \eqn{upper} bound).
#' @param verbose Optional. A boolean. If \code{TRUE}, prints out the lambda value at each iteration.
#' @param tol Optional. A number, the tolerance parameter.
#' @param maxit Optional. A positive integer, the maximum number of iterations in model fitting for each lambda.
#' @param cur_res Optional. A list, current results returned from a previous lambda. If provided, used as a warm start. Default to \code{NULL}.
#' @return
#'     \item{lambda}{A number, the best \code{lambda} that produces the desired number of edges. \code{1e-10} or \code{1e15} is returned if out of bound.}
#'     \item{cur_res}{A list, results for this \code{lambda}. May be \code{NULL} if \code{lambda} is out of bound.}
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()}, as the
#' #   way to call this function (\code{test_lambda_bounds()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                  centered=FALSE, profiled=FALSE, diag=dm)
#' lambda_cur_res <- test_lambda_bounds(elts_gauss_np, "symmetric", lambda=1,
#'                   lambda_ratio=1, step=1.5, lower=TRUE, cur_res=NULL)
#' lambda_cur_res2 <- test_lambda_bounds(elts_gauss_np, "symmetric", lambda=1,
#'                   lambda_ratio=1, step=1.5, lower=FALSE, cur_res=lambda_cur_res$cur_res)
#'@export
test_lambda_bounds <- function(elts, symmetric, lambda=1, lambda_ratio=1, step=2, lower = TRUE, verbose=TRUE, tol=1e-6, maxit=10000, cur_res=NULL){
  if (step <= 1)
    stop("step must be larger than 1.")
  if (verbose)
    print(paste("Testing lower bound for lambda: ", lambda))
  cur_res <- get_results(elts, symmetric, lambda1=lambda, lambda2=lambda/lambda_ratio, tol=tol, maxit=maxit, previous_res=cur_res, is_refit=FALSE)
  want_edges <- ifelse(lower, elts$p*(elts$p-1), 0)
  best_lambda <- best_res <- NULL
  if (length(cur_res$edges) == want_edges){
    best_lambda <- lambda; best_res <- cur_res
  }
  ## Search direction (+1 or -1) will never change, since otherwise we will go to the previous lambda we have just examined
  ## If finding a tight (i.e. highest possible) lower bound and the lambda given already gives the complete graph, then increase lambda; if the lambda does not give the complete graph (not small enough), then decrease lambda
  ## If finding a tight (i.e. loweest possible) higher bound and the lambda given already gives the empty graph, then lower lambda; if the lambda does not give the empty graph (not large enough), then increase lambda
  search_direction <- 2*(xor(lower, is.null(best_lambda))) - 1
  while (TRUE){
    lambda <- lambda * step ^ search_direction
    if (verbose)
      message("Testing lower bound for lambda:", lambda)
    cur_res <- get_results(elts, symmetric=symmetric, lambda1=lambda, lambda2=lambda/lambda_ratio, tol=tol, maxit=maxit, previous_res=cur_res)
    if (verbose)
      message(", number of edges=", length(cur_res$edges), "; want ", want_edges, ".\n", sep="")
    if (length(cur_res$edges) == want_edges){ ## If this lambda gives the desired graph (complete/empty), it is necessarily the best lambda so far
      if (is.null(best_lambda)) ## If no previous lambda is good (i.e. the initial lambda is bad), return this, as the next lambda will necessarily give the right graph but will not be tight
        return (list("lambda"=lambda, "cur_res"=cur_res))
      best_lambda <- lambda; best_res <- cur_res  ## Otherwise, just update the best lambda
    } else if (!is.null(best_lambda)) ## If this lambda does not give the desired graph, but we already found some (best) lambda that does, return that lambda; otherwise keep searching
      return (list("lambda"=best_lambda, "cur_res"=best_res))
    if (lambda <= 1e-10 || lambda >= 1e15){ ## If too small or too large, have to stop
      if (verbose)
        message("Stopped at ", max(1e-10, min(lambda, 1e15)), ".\n", sep = "")
      return (list("lambda"=max(1e-10, min(lambda, 1e15)))) ## cur_res=NULL
    }
  }
}

#' Searches for a tight bound for \eqn{\lambda_{\boldsymbol{K}}}{\lambda_K} that gives the empty or complete graph starting from a given lambda
#'
#' Searches for the smallest lambda that gives the empty graph (if \code{lower == FALSE}) or the largest that gives the complete graph (if \code{lower == TRUE}) starting from the given lambda.
#'
#' @param elts A list, elements necessary for calculations returned by get_elts().
#' @param symmetric A string. If equals \code{"symmetric"}, estimates the minimizer \eqn{\mathbf{K}}{K} over all symmetric matrices; if \code{"and"} or \code{"or"}, use the "and"/"or" rule to get the support
#' @param lambda_ratio A positive number (or \code{Inf}), the fixed ratio \eqn{\lambda_{\mathbf{K}}}{\lambda_K} and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta}, if \eqn{\lambda_{\boldsymbol{\eta}}\neq 0}{\lambda_\eta!=0} (non-profiled) in the non-centered setting.
#' @param lower A boolean. If \code{TRUE}, finds the largest possible lambda that gives the complete graph (a \eqn{lower} bound). If \code{FALSE}, finds the smallest possible lambda that gives the empty graph (an \eqn{upper} bound).
#' @param verbose Optional.  A boolean. If \code{TRUE}, prints out the lambda value at each iteration.
#' @param tol Optional. A number, the tolerance parameter.
#' @param maxit Optional. A positive integer, the maximum number of iterations in model fitting for each lambda.
#' @param lambda_start Optional. A number, the starting point for searching. If \code{NULL}, set to \code{1e-4} if \code{lower == TRUE}, or \code{1} if \code{lower == FALSE}.
#' @return A number, the best lambda that produces the desired number of edges. \code{1e-10} or \code{1e15} is returned if out of bound.
#' @details This function calls \code{test_lambda_bounds} three times with \code{step} set to \code{10}, \code{10^(1/4)}, \code{10^(1/16)}, respectively.
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()}, as the
#' #   way to call this function (\code{test_lambda_bounds2()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' test_lambda_bounds2(elts_gauss_np, "symmetric", lambda_ratio=2,
#'      lower=TRUE, verbose=TRUE, lambda_start=NULL)
#' test_lambda_bounds2(elts_gauss_np, "symmetric", lambda_ratio=2,
#'      lower=FALSE, verbose=TRUE, lambda_start=1.0)
#' @export
test_lambda_bounds2 <- function(elts, symmetric, lambda_ratio=Inf, lower = TRUE, verbose = TRUE, tol=1e-6, maxit=10000, lambda_start = NULL){
  ## To-do: boundary points tend to be evaluated twice
  ## To-do find the open-interval bounds, i.e. smallest discrete lambdas s.t. num_edges > 0 and < max_edges
  if (!is.null(lambda_start)) {
    if (lambda_start <= 0) stop ("lambda_start must be positive if provided.")
    lambda_cur_res <- list("lambda"=lambda_start)
  } else
    lambda_cur_res <- list("lambda"=ifelse(lower, 1e-4, 1))
  step <- 10
  while (step >= 1.5^(1/4)){
    if (verbose)
      message("step =", step, "\n")
    lambda_cur_res <- test_lambda_bounds(elts, symmetric=symmetric, lambda=lambda_cur_res$lambda, lambda_ratio=lambda_ratio, step=step, lower=lower, verbose=verbose, tol=tol, maxit=maxit, cur_res=lambda_cur_res$cur_res)
    if (lambda_cur_res$lambda <= 1e-10 || lambda_cur_res$lambda >= 1e15) {
      return (max(1e-10, min(lambda_cur_res$lambda, 1e15)))
    }
    step <- step ^ (1/4)
  }
  if (verbose)
    message("Final: ", lambda_cur_res$lambda, ", ", length(lambda_cur_res$cur_res$edges), " edges.\n", sep="")
  return (lambda_cur_res$lambda)
}

#' Analytic solution for the minimum \eqn{\lambda_{\mathbf{K}}}{\lambda_K} that gives the empty graph.
#'
#' Analytic solution for the minimum \eqn{\lambda_{\mathbf{K}}}{\lambda_K} that gives the empty graph. In the non-centered setting the bound is not tight, as it is such that both \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta} are empty. The bound is also not tight if \code{symmetric == "and"}.
#'
#' @param elts A list, elements necessary for calculations returned by \code{get_elts()}.
#' @param symmetric A string. If equals \code{"symmetric"}, estimates the minimizer \eqn{\mathbf{K}}{K} over all symmetric matrices; if \code{"and"} or \code{"or"}, use the "and"/"or" rule to get the support.
#' @param lambda_ratio A positive number (or \code{Inf}), the fixed ratio \eqn{\lambda_{\mathbf{K}}}{\lambda_K} and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta}, if \eqn{\lambda_{\boldsymbol{\eta}}\neq 0}{\lambda_\eta!=0} (non-profiled) in the non-centered setting.
#' @return A number, the smallest lambda that produces the empty graph in the centered case, or that gives zero solutions for \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta} in the non-centered case. If \code{symmetric == "and"}, it is not a tight bound for the empty graph.
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()},
#' #   as the way to call this function (\code{lambda_max()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#'
#' # Exact analytic solution for the smallest lambda such that K and eta are both zero,
#' #  but not a tight bound for K ONLY
#' lambda_max(elts_gauss_np, "symmetric", 2)
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_gauss_np, "symmetric", lambda_ratio=2, lower = FALSE,
#'      lambda_start = lambda_max(elts_gauss_np, "symmetric", 2))
#'
#' # Exact analytic solution for the smallest lambda such that K and eta are both zero,
#' #  but not a tight bound for K ONLY
#' lambda_max(elts_gauss_np, "or", 2)
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_gauss_np, "or", lambda_ratio=2, lower = FALSE,
#'      lambda_start = lambda_max(elts_gauss_np, "or", 2))
#'
#' # An upper bound, not tight.
#' lambda_max(elts_gauss_np, "and", 2)
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_gauss_np, "and", lambda_ratio=2, lower = FALSE,
#'      lambda_start = lambda_max(elts_gauss_np, "and", 2))
#'
#'
#' elts_gauss_p <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'               centered=FALSE, profiled=TRUE, diag=dm)
#' # Exact analytic solution
#' lambda_max(elts_gauss_p, "symmetric")
#' # Numerical solution, should be close to the analytic solution
#' test_lambda_bounds2(elts_gauss_p, "symmetric", lambda_ratio=Inf, lower = FALSE,
#'      lambda_start = NULL)
#'
#' # Exact analytic solution
#' lambda_max(elts_gauss_p, "or")
#' # Numerical solution, should be close to the analytic solution
#' test_lambda_bounds2(elts_gauss_p, "or", lambda_ratio=Inf, lower = FALSE,
#'      lambda_start = NULL)
#'
#' # An upper bound, not tight
#' lambda_max(elts_gauss_p, "and")
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_gauss_p, "and", lambda_ratio=Inf, lower = FALSE,
#'      lambda_start = lambda_max(elts_gauss_p, "and"))
#' @export
lambda_max <- function(elts, symmetric, lambda_ratio=Inf){
  p <- elts$p
  tol <- 1e-10
  if (elts$setting == "gaussian" && elts$domain_type == "R") {
    diag(elts$Gamma_K) <- elts$diagonals_with_multiplier
    K_init_diag <- 1/diag(elts$Gamma_K)
    if (symmetric == "symmetric"){
      K_init_off_max <- max(sapply(1:(p-1), function(i){max(abs(-elts$Gamma_K[(i+1):p,i]*K_init_diag[i] - elts$Gamma_K[i,(i+1):p]*K_init_diag[(i+1):p]))}))/2
    } else {
      K_init_off_max <- max(sapply(1:p, function(i){max(abs(-elts$Gamma_K[-i,i]*K_init_diag[i]))}))
    }
    if (elts$centered || elts$profiled_if_noncenter)
      return (max(K_init_off_max, tol))  ## If symmetric!="and", this lambda is exact when centered or profiled
    else {
      if (is.infinite(lambda_ratio))
        stop ("elts$profiled_if_noncenter should be TRUE if lambda_ratio=Inf.")
      eta_init <- max(abs(sapply(1:p, function(i){-K_init_diag[i]*elts$Gamma_K_eta[i]})))
      return (max(K_init_off_max, eta_init * lambda_ratio, tol)) ## This only serves an upper bound, since it is the max lambda where both K and eta are 0.
    }
  } else if (elts$domain_type == "simplex") {
    if (symmetric != "symmetric") #### This should be updated in the future
      stop("Asymmetric estimation not supported for models on the simplex.")
    if (!elts$centered && elts$profiled_if_noncenter)
      stop("Profiled non-centered estimators not supported for models on the simplex.")
    K_grad_from_eta <- function(elts, eta) { # \Gamma_{K,eta} %*% \eta in paper
      return (c(sapply(1:(p-1), function(i){elts$Gamma_K_eta[,i]*eta[i] + elts$Gamma_Kj_etap[,i]*eta[p]}),
                elts$Gamma_K_eta[,p]*eta[p] + elts$Gamma_Kp_etaj %*% eta[-p]))
    }
    K_grad_from_Kdiag <- function(elts, Kdiag) { # \Gamma_K %*% K in paper, assuming only the diagonals of K are non-zero
      return (c(sapply(1:(p-1), function(i){
        elts$Gamma_K[,(i-1)*p+i]*Kdiag[i] + elts$Gamma_K_jp[,i*p]*Kdiag[p]}),
        rowSums(sapply(1:(p-1), function(i){elts$Gamma_K_jp[i,(i-1)*p+1:p]*Kdiag[i]})) +
          elts$Gamma_K[,p*p]*Kdiag[p]))
    }
    if (grepl("sum0", elts$setting)) { # If K assumed to sum to 0, then lambda_max corresponds to K = 0.
      if (!elts$centered && is.infinite(lambda_ratio)) { # Non-centered with no penalty on eta (in which case we only care about when K is all 0)
        Gamma_eta <- diag(elts$Gamma_eta)
        Gamma_eta[p, -p] <- Gamma_eta[-p, p] <- elts$Gamma_eta_jp
        eta_init <- solve(Gamma_eta, elts$g_eta)
        K_grad <- matrix(K_grad_from_eta(elts, eta_init) - elts$g_K, p, p)
        K_grad[p,-p] <- K_grad[p,-p] / (p-1)
        K_grad[,p] <- K_grad[,p] / (p-1)
        return (max(max(abs(K_grad + t(K_grad))) / 2, tol))
      } else { # Centered or has penalty on eta
        elts$g_K <- matrix(elts$g_K, nrow=elts$p, ncol=elts$p)
        elts$g_K[p,-p] <- elts$g_K[p,-p] / (p-1)
        elts$g_K[,p] <- elts$g_K[,p] / (p-1)
        gK_max <- max(abs(elts$g_K + t(elts$g_K))) / 2
        if (elts$centered) # Centered
          return (max(gK_max, tol))
        else # Non-centered and has penalty on eta
          return (max(gK_max, max(abs(elts$g_eta[-p]), abs(elts$g_eta[p]) / (p-1)) * lambda_ratio, tol))
      }
    } else {
      elts$Gamma_K[(0:(p^2-1))*p+rep(1:p,p)] <- elts$diagonals_with_multiplier
      if (!elts$centered && is.infinite(lambda_ratio)) { # Non-centered with no penalty on eta (in which case we only care about when K is all 0)
        # Gamma_K_AND_eta: form the Gamma matrix restricted to K_diag and eta
        Gamma_K_AND_eta <- diag(c(elts$diagonals_with_multiplier[(0:(p-1))*p+1:p], elts$Gamma_eta)) # quadratic terms
        Gamma_K_AND_eta[p,1:(p-1)] <- Gamma_K_AND_eta[1:(p-1),p] <- elts$Gamma_K_jp[1:(p-1)*(p^2+1)-p] # elts$Gamma_K_jp[i,i*100] for i in 1:(p-1), interaction between Kjj and Kpp
        Gamma_K_AND_eta[2*p, (p+1):(2*p-1)] <- Gamma_K_AND_eta[(p+1):(2*p-1), 2*p] <- elts$Gamma_eta_jp # interaction between etaj and etap
        Gamma_K_AND_eta[(1:p)*(2*p+1)-2*p+p*p*2] <- Gamma_K_AND_eta[(1:p)*(2*p+1)-2*p+p] <- diag(elts$Gamma_K_eta) # interaction between Kjj and etaj
        Gamma_K_AND_eta[1:(p-1),2*p] <- Gamma_K_AND_eta[2*p,1:(p-1)] <- diag(elts$Gamma_Kj_etap) # interaction between Kjj and etap
        Gamma_K_AND_eta[p,(p+1):(2*p-1)] <- Gamma_K_AND_eta[(p+1):(2*p-1),p] <- elts$Gamma_Kp_etaj[p,]
        Kdiag_eta_init <- solve(Gamma_K_AND_eta, c(elts$g_K[(0:(p-1))*p+1:p], elts$g_eta))
        K_grad <- matrix(K_grad_from_Kdiag(elts, Kdiag_eta_init[1:p]) +
                           K_grad_from_eta(elts, Kdiag_eta_init[p+1:p]) - elts$g_K, p, p)
        K_grad[p,-p] <- K_grad[p,-p] / (p-1)
        K_grad[,p] <- K_grad[,p] / (p-1)
        return (max(max(abs(K_grad + t(K_grad))) / 2, tol))
      } else { # Centered or has penalty on eta
        Gamma_K <- diag(elts$diagonals_with_multiplier[(0:(p-1))*p+1:p])
        Gamma_K[p,-p] <- Gamma_K[-p,p] <- elts$Gamma_K_jp[1:(p-1)*(p^2+1)-p] # elts$Gamma_K_jp[i,i*100] for i in 1:(p-1)
        K_init_diag <- solve(Gamma_K, elts$g_K[(0:(p-1))*p+1:p])
        K_grad <- matrix(K_grad_from_Kdiag(elts, K_init_diag) - elts$g_K, p, p)
        K_grad[p,-p] <- K_grad[p,-p] / (p-1)
        K_grad[,p] <- K_grad[,p] / (p-1)
        gK_max <- max(abs(K_grad + t(K_grad))) / 2
        if (elts$centered) { # Centered
          return (max(gK_max, tol))
        } else { # Non-centered and has penalty on eta
          g_eta <- c(diag(elts$Gamma_K_eta)[-p] * K_init_diag[-p] + elts$Gamma_Kp_etaj[p,] * K_init_diag[p],
                     sum(diag(elts$Gamma_Kj_etap) * K_init_diag[-p]) + elts$Gamma_K_eta[p,p] * K_init_diag[p]) - elts$g_eta
          g_eta[p] <- g_eta[p] / (p-1)
          return (max(gK_max, max(abs(g_eta)) * lambda_ratio, tol))
        }
      }
    }
  } else {
    elts$Gamma_K[(0:(p^2-1))*p+rep(1:p,p)] <- elts$diagonals_with_multiplier
    K_init_diag <- elts$g_K[(0:(p-1))*p+1:p]/elts$Gamma_K[((0:(p-1))*p+0:(p-1))*p+1:p]
    if (symmetric == "symmetric"){
      K_init_off_max <- max(sapply(1:(p-1), function(i){max(abs(-elts$Gamma_K[(i-1)*p^2+(i-1)*p+(i+1):p]*K_init_diag[i] - elts$Gamma_K[(i:(p-1))*p^2+(i:(p-1))*p+i]*K_init_diag[(i+1):p] + elts$g_K[(i-1)*p+(i+1):p] + elts$g_K[(i:(p-1))*p+i]))}))/2
    } else {
      K_init_off_max <- max(sapply(1:p, function(i){noi<-(1:p)[-i]; max(abs(-elts$Gamma_K[(i-1)*p^2+(i-1)*p+noi]*K_init_diag[i]+elts$g_K[(i-1)*p+noi]))}))
    }
    if (elts$centered || elts$profiled_if_noncenter)
      return (K_init_off_max+1e-8)  ## If symmetric!="and", this lambda is exact when centered or profiled
    else {
      if (is.infinite(lambda_ratio))
        stop ("elts$profiled_if_noncenter should be TRUE if lambda_ratio=Inf.")
      eta_init <- max(abs(sapply(1:p, function(i){elts$g_eta[i]-K_init_diag[i]*elts$Gamma_K_eta[(i-1)*p+i]})))
      return (max(K_init_off_max, eta_init * lambda_ratio)+1e-8) ## This only serves an upper bound, since it is the max lambda where both K and eta are 0.
    }
  }
}

#' Helper function for outputting if verbose.
#' @param out Text string.
#' @param verbose Boolean.
#' @param verbosetext Text string.
#' @return If \code{verbose == TRUE}, outputs a string that concatenates \code{verbosetext} and \code{out}.
#' @export
s_output <- function(out, verbose, verbosetext){
  if (verbose)
    message(verbosetext, out, "\n")
}

#' Helper function for making fold IDs for cross validation.
#' @param nsamp Number of samples.
#' @param nfold Number of cross validation folds.
#' @param cv_fold_seed Seed for random shuffling.
#' @return A list of \code{nsamp} vectors, numbers 1 to \code{nsamp} shuffled and groupped into vectors of length \code{floor(nsamp/nfold)} followed by vectors of length \code{floor(nsamp/nfold)+1}.
#' @examples
#' make_folds(37, 5, NULL)
#' make_folds(100, 5, 2)
#' make_folds(100, 10, 3)
#' @export
make_folds <- function(nsamp, nfold, cv_fold_seed){
  if (!is.null(cv_fold_seed)) set.seed(cv_fold_seed)
  id <- sample(nsamp, nsamp)
  every <- floor(nsamp / nfold)
  ends <- cumsum(c(0,rep(every, nfold-(nsamp-nfold*every)), rep(every+1, nsamp-nfold*every)))
  return (lapply(1:nfold, function(i){id[(ends[i]+1):ends[i+1]]}))
}

#' The main function for the generalized score-matching estimator for graphical models.
#'
#' The main function for the generalized score-matching estimator for graphical models.
#'
#' @param x An \code{n} by \code{p} matrix, the data matrix, where \code{n} is the sample size and \code{p} the dimension.
#' @param setting A string that indicates the distribution type, must be one of \code{"exp"}, \code{"gamma"}, \code{"gaussian"}, \code{"log_log"}, \code{"log_log_sum0"}, or of the form \code{"ab_NUM1_NUM2"}, where \code{NUM1} is the \code{a} value and \code{NUM2} is the \code{b} value, and \code{NUM1} and \code{NUM2} must be integers or two integers separated by "/", e.g. "ab_2_2", "ab_2_5/4" or "ab_2/3_1/2".
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @param elts A list (optional), elements necessary for calculations returned by get_elts().
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param symmetric A string. If equals \code{"symmetric"}, estimates the minimizer \eqn{\mathbf{K}}{K} over all symmetric matrices; if \code{"and"} or \code{"or"}, use the "and"/"or" rule to get the support. Default to \code{"symmetric"}.
#' @param scale A string indicating the scaling method. If contains \code{"sd"}, columns are scaled by standard deviation; if contains \code{"norm"}, columns are scaled by l2 norm; if contains \code{"center"} and \code{setting == "gaussian" && domain$type == "R"}, columns are centered to have mean zero. Default to \code{"norm"}.
#' @param lambda1s A vector of lambdas, the penalty parameter for K.
#' @param lambda_length An integer >= 2, the number of lambda1s. Ignored if \code{lambda1s} is provided, otherwise a grid of lambdas is automatically chosen so that the results range from an empty graph to a complete graph. Default to \code{10} if neither \code{lambda1s} nor \code{lambda_length} is provided.
#' @param lambda_ratio A positive number, the fixed ratio between \eqn{\lambda_{\mathbf{K}}}{\lambda_K} and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta}, if \eqn{\lambda_{\boldsymbol{\eta}}\neq 0}{\lambda_\eta!=0} (non-profiled) in the non-centered setting.
#' @param mode A string, the class of the \code{h} function. Ignored if \code{elts}, or \code{h} and \code{hp} are provided, or if \code{setting == "gaussian" && domain$type == "R"}.
#' @param param1 A number, the first parameter to the \code{h} function. Ignored if \code{elts}, or \code{h} and \code{hp} are provided, or if \code{setting == "gaussian" && domain$type == "R"}.
#' @param param2 A number, the second parameter (may be optional depending on \code{mode}) to the \code{h} function. Ignored if \code{elts}, or \code{h} and \code{hp} are provided, or if \code{setting == "gaussian" && domain$type == "R"}.
#' @param h_hp A function that returns a list containing \code{hx=h(x)} (element-wise) and \code{hpx=hp(x)} (element-wise derivative of \eqn{h}) when applied to a vector or a matrix \code{x}, both of which has the same shape as \code{x}.
#' @param unif_dist Optional, defaults to \code{NULL}. If not \code{NULL}, \code{h_hp} must be \code{NULL} and \code{unif_dist(x)} must return a list containing \code{"g0"} of length \code{nrow(x)} and \code{"g0d"} of dimension \code{dim(x)}, representing the l2 distance and the gradient of the l2 distance to the boundary: the true l2 distance function to the boundary is used for all coordinates in place of h_of_dist; see "Estimating Density Models with Complex Truncation Boundaries" by Liu et al, 2019. That is, \eqn{(h_j\circ \phi_j)(x_i)}{(h_j\circ phi_j)(xi)} in the score-matching loss is replaced by \eqn{g_0(x_i)}{g0(xi)}, the l2 distance of xi to the boundary of the domain.
#' @param verbose Optional. A boolean, whether to output intermediate results.
#' @param verbosetext Optional. A string, text to be added to the end of each printout if \code{verbose == TRUE}.
#' @param tol Optional. A number, the tolerance parameter. Default to \code{1e-6}.
#' @param maxit Optional. A positive integer, the maximum number of iterations for each fit. Default to \code{1000}.
#' @param BIC_refit A boolean, whether to get the BIC scores by refitting an unpenalized model restricted to the estimated edges, with \code{lambda1=lambda2=0} and \code{diagonal_multiplier=1}. Default to \code{TRUE}.
#' @param warmstart Optional. A boolean, whether to use the results from a previous (larger) lambda as a warm start for each new lambda. Default to \code{TRUE}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier. Optional and ignored if elts is provided. If \code{ncol(x) > ncol(n)}, a value strictly larger than 1 is recommended. Default to \eqn{1+\left(1-\left(1+4e\max\left(6\log p/n, \sqrt{6\log p/n}\right)\right)^{-1}\right)}{1+(1-1/(1+4e*max(6*log(p)/n, sqrt(6*log(p)/n))))}.
#' @param eBIC_gammas Optional. A number of a vector of numbers. The \eqn{\gamma} parameter in eBIC. Default to \code{c(0,0.5,1)}.
#' @param cv_fold Optional. An integer larger than 1 if provided. The number of folds used for cross validation. If provided, losses will be calculated on each fold with model fitted on the other folds, and a \code{lambda_length x cv_fold} matrix \code{cv_losses} will be returned.
#' @param cv_fold_seed Optional. Seed for generating folds for cross validation.
#' @param return_raw A boolean, whether to return the raw estimates of \code{K}. Default to \code{FALSE}.
#' @param return_elts A boolean, whether to return the \code{elts} used for estimation. Default to \code{FALSE}.
#' @return
#'    \item{edgess}{A list of vectors of integers: indices of the non-zero edges.}
#'    \item{BICs}{A \code{lambda_length} by \code{length(eBIC_gammas)} matrix of raw eBIC scores (without refitting). If \code{return_raw == FALSE}, may contain \code{Inf}s for rows after the first lambda that gives the complete graph.}
#'    \item{lambda1s}{A vector of numbers of length \code{lambda_length}: the grid of \code{lambda1}s over which the estimates are obtained.}
#'    \item{converged}{A vector of booleans of length \code{lambda_length}: indicators of convergence for each fit. If \code{return_raw == FALSE}, may contain \code{0}s for all lambdas after the first lambda that gives the complete graph.}
#'    \item{iters}{A vector of integers of length \code{lambda_length}: the number of iterations run for each fit. If \code{return_raw == FALSE}, may contain \code{0}s for all lambdas after the first lambda that gives the complete graph.}
#'
#'    In addition,
#'    if \code{centered == FALSE},
#'    \item{etas}{A \code{lambda_length}*\code{p} matrix of \code{eta} estimates with the \eqn{i}-th row corresponding to the \eqn{i}-th \code{lambda1}. If \code{return_raw == FALSE},  may contain \code{NA}s after the first lambda that gives the complete graph.}
#'    if \code{centered == FALSE} and non-profiled,
#'    \item{lambda2s}{A vector of numbers of length \code{lambda_length}: the grid of \code{lambda2}s over which the estimates are obtained.}
#'    if \code{return_raw == TRUE},
#'    \item{raw_estimate}{A list that contains \code{lambda_length} estimates for \code{K} of size \code{ncol(x)}*\code{ncol(x)}.}
#'    if \code{BIC_refit == TRUE},
#'    \item{BIC_refits}{A \code{lambda_length} by \code{length(eBIC_gammas)} matrix of refitted eBIC scores, obtained by refitting unpenalized models restricted to the estimated edges. May contain \code{Inf}s for rows after the first lambda that gives the graph restricted to which an unpenalized model does not have a solution (loss unbounded from below).}
#'    if \code{cv_fold} is not \code{NULL},
#'    \item{cv_losses}{A \code{lambda_length x cv_fold} matrix of cross validation losses. If \code{return_raw == FALSE}, may contain \code{Inf}s for all lambdas after the first lambda that gives the complete graph.}
#'    if \code{return_elts == TRUE},
#'    \item{elts}{A list of elements returned from \code{get_elts()}.}
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()},
#' #   as the way to call this function (\code{estimate()}) is exactly the same in those cases.
#' n <- 30
#' p <- 20
#' domain <- make_domain("R+", p=p)
#' mu <- rep(0, p)
#' K <- diag(p)
#' lambda1s <- c(0.01,0.1,0.2,0.3,0.4,0.5)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' ## Centered estimates, no elts or h provided, mode and params provided
#' est1 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=TRUE,
#'           symmetric="symmetric", lambda1s=lambda1s, mode="min_pow",
#'           param1=1, param2=3, diag=dm, return_raw=TRUE)
#'
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' ## Centered estimates, no elts provided, h provided; equivalent to est1
#' est2 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=TRUE,
#'           symmetric="symmetric", lambda1s=lambda1s, h_hp=h_hp, diag=dm, return_raw=TRUE)
#' compare_two_results(est1, est2) ## Should be almost all 0
#'
#' elts_gauss_c <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'             centered=TRUE, diag=dm)
#' ## Centered estimates, elts provided; equivalent to est1 and est2
#' ## Here diagonal_multiplier will be set to the default value, equal to dm above
#' est3 <- estimate(x, "gaussian", domain=domain, elts=elts_gauss_c,
#'           symmetric="symmetric", lambda1s=lambda1s, diag=NULL,
#'           return_raw=TRUE)
#' compare_two_results(est1, est3) ## Should be almost all 0
#'
#' ## Noncentered estimates with Inf penalty on eta; equivalent to est1~3
#' est4 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=FALSE,
#'           lambda_ratio=0, symmetric="symmetric", lambda1s=lambda1s,
#'           h=h_hp, diag=dm, return_raw=TRUE)
#' sum(abs(est4$etas)) ## Should be 0 since non-centered with lambda ratio 0 is equivalent to centered
#' est4$etas <- NULL ## But different from est1 in that the zero etas are returned in est4
#' compare_two_results(est1, est4) ## Should be almost all 0
#'
#'
#' ## Profiled estimates, no elts or h provided, mode and params provided
#' est5 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=FALSE,
#'           lambda_ratio=Inf, symmetric="or", lambda1s=lambda1s,
#'           mode="min_pow", param1=1, param2=3, diag=dm, return_raw=TRUE)
#'
#' ## Profiled estimates, no elts provided, h provided; equivalent to est5
#' est6 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=FALSE,
#'           lambda_ratio=Inf, symmetric="or", lambda1s=lambda1s,
#'           h_hp=h_hp, diag=dm, return_raw=TRUE)
#' compare_two_results(est5, est6) ## Should be almost all 0
#'
#' elts_gauss_p <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=TRUE, diag=dm)
#' ## Profiled estimates, elts provided; equivalent to est5~6
#' est7 <- estimate(x, "gaussian", domain=domain, elts=elts_gauss_p, centered=FALSE,
#'           lambda_ratio=Inf, symmetric="or", lambda1s=lambda1s,
#'           diagonal_multiplier=NULL, return_raw=TRUE)
#' compare_two_results(est5, est7) ## Should be almost all 0
#'
#'
#' ## Non-centered estimates, no elts or h provided, mode and params provided
#' ## Using 5-fold cross validation and no BIC refit
#' est8 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=FALSE,
#'           lambda_ratio=2, symmetric="and", lambda_length=100,
#'           mode="min_pow", param1=1, param2=3, diag=dm, return_raw=TRUE,
#'           BIC_refit=FALSE, cv_fold=5, cv_fold_seed=2)
#'
#' ## Non-centered estimates, no elts provided, h provided; equivalent to est5
#' ## Using 5-fold cross validation and no BIC refit
#' est9 <- estimate(x, "gaussian", domain=domain, elts=NULL, centered=FALSE,
#'           lambda_ratio=2, symmetric="and", lambda_length=100, h_hp=h_hp,
#'           diag=dm, return_raw=TRUE, BIC_refit=FALSE, cv_fold=5, cv_fold_seed=2)
#' compare_two_results(est8, est9) ## Should be almost all 0
#'
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain, centered=FALSE,
#'                 profiled=FALSE, diag=dm)
#' ## Non-centered estimates, elts provided; equivalent to est8~9
#' ## Using 5-fold cross validation and no BIC refit
#' est10 <- estimate(x, "gaussian", domain, elts=elts_gauss_np, centered=FALSE,
#'            lambda_ratio=2, symmetric="and", lambda_length=100, diag=NULL,
#'            return_raw=TRUE, BIC_refit=FALSE, cv_fold=5, cv_fold_seed=2)
#' compare_two_results(est8, est10) ## Should be almost all 0
#'
#' @export
estimate <- function(x, setting, domain, elts=NULL, centered=TRUE, symmetric="symmetric", scale="",
                     lambda1s=NULL, lambda_length=NULL, lambda_ratio=Inf, mode=NULL, param1=NULL,
                     param2=NULL, h_hp=NULL, unif_dist=NULL, verbose=TRUE, verbosetext="", tol=1e-6,
                     maxit=1000, BIC_refit=TRUE, warmstart=TRUE, diagonal_multiplier=NULL,
                     eBIC_gammas=c(0,0.5,1), cv_fold=NULL, cv_fold_seed=NULL, return_raw=FALSE,
                     return_elts=FALSE){
  ## BIC_refit: calculate BIC (with refit) or not
  ## return_raw: return the raw estimates or not
  ## If elts is given, centered, scale, mode, param1, param2, h, hp are all ignored
  ## If both h and hp given, mode, param1, param2 are ignored; h must be continuous and positive almost thinningwhere and h(0)=0; hp must be the almost thinningwhere derivative of h
  ## Must provide at least one of: 1. elts, 2. mode & param1, 3. h & hp
  n <- nrow(x); p <- ncol(x)
  if (is.null(diagonal_multiplier))
    diagonal_multiplier <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
  if (!is.null(elts) && setting != elts$setting)
    stop(paste("The setting you chose was ", setting, ", but elts$setting was ", elts$setting, ".", sep=""))

  changed_from_nc_to_c <- FALSE
  if (lambda_ratio == 0) {
    centered <- changed_from_nc_to_c <- TRUE;
    if (!is.null(elts)){elts$centered <- TRUE; elts$Gamma_K_eta <- elts$Gamma_eta <- elts$g_eta <- NULL}
  } else if (!is.null(elts) && elts$centered != centered) {
    stop ("centered and elts$centered do not match.")
  } else if (is.infinite(lambda_ratio) && !is.null(elts) && !elts$centered && !elts$profiled_if_noncenter){
    stop("elts should be created with profiled_if_noncenter = TRUE if lambda_ratio is Inf and elts$centered is FALSE.")
  } else if (is.finite(lambda_ratio) && !is.null(elts) && !elts$centered && elts$profiled_if_noncenter){
    stop("elts should be created with profiled_if_noncenter = FALSE if lambda_ratio is finite and elts$centered is FALSE.")
  }
  if (is.null(diagonal_multiplier)){ ## If diagonal multiplier not provided
    if (is.null(elts) || is.null(elts$diagonal_multiplier))
      stop("diagonal_multiplier must be provided if elts is not provided or does not contain a diagonal multiplier.")
    diagonal_multiplier <- elts$diagonal_multiplier
  } else {
    if (!is.null(elts) && !is.null(elts$diagonal_multiplier) && abs(diagonal_multiplier - elts$diagonal_multiplier) > tol){
      warning("diagonal_multiplier and elts$diagonal_multiplier do not agree. Using elts$diagonal_multiplier instead.")
      diagonal_multiplier <- elts$diagonal_multiplier
    } else {
      if (diagonal_multiplier < 1 + tol && p > n){ ## ill-behaved
        warning("p > n and diagonal_multiplier should be larger than 1.")
      } else if (diagonal_multiplier > 1 - tol){ ## Allows numerical error; if diagonal_multiplier slightly below 1, set it to 1.
        diagonal_multiplier <- max(1, diagonal_multiplier)
      } else
        stop("diagonal_multiplier must be at least 1.")
    }
  }
  if (!symmetric %in% c("symmetric", "and", "or"))
    stop("Parameter symmetric must be one of \"symmetric\", \"and\", or \"or\".")
  if (scale != "") {
    if (domain$type %in% c("simplex", "uniform"))
      stop("Simplex- and uniform-type domains are not invariate to scaling. Please set \"scale\" to \"\".")
    if (domain$type == "polynomial")
      warning("Polynomial-type domains may not be invariate to scaling. Please set \"scale\" to \"\" if this is the case, otherwise you may ignore this warning.")
  }

  if (!is.null(lambda1s)){
    if (!is.null(lambda_length) && length(lambda1s) != lambda_length)
      warning("lambda1s should have be a vector with the same length as lambda_length. lambda_length ignored.\n")
    if (any(lambda1s < 0))  ## If lambda1s given but some are negative
      stop("All lambda1s must be non-negative.")
  } else if (is.null(lambda_length)){  ## If neither lambda1s or lambda_length given
    lambda_length = 10
    warning("No lambda1s or lambda_length given. Automatically choosing 10 lambda1s.\n")
  } else if (as.integer(lambda_length) <= 2){
    stop("If lambda1s are not provided, lambda_length must be at least 2.")
  } else
    lambda_length <- as.integer(lambda_length)
  if ((setting != "gaussian" || domain$type != "R") && is.null(elts) && (is.null(mode) || is.null(param1)) && (is.null(h_hp)) && is.null(unif_dist))
    stop("At least one of 1. elts, 2. mode & param1, 3. h_hp, or 4. unif_dist has to be provided.")
  if ((setting != "gaussian" || domain$type != "R") && is.null(elts) && !is.null(mode)){
    if (!is.null(h_hp)){
      warning("mode and h_hp should not be provided at the same time. Using h_hp instead.\n")
      #if (abs(h(0)) > tol)
      #  stop("h(0) must be equal to 0.")
    } else{
      if (is.null(param1))
        stop("param1 (and param2 optionally) must be provided with mode.")
      h_hp <- get_h_hp(mode=mode, para=param1, para2=param2);
      if (length(h_hp(matrix(1:4, nrow=2)))==0)
        stop("Error occurred in generating h_hp. Possibly due to invalid param1 and/or param2.")
      #if (abs(h(0)) > tol)
      #  stop(paste("h(0)=", h(0), ", larger than 0. Stopped.", sep=""))
    }
  }
  if (!is.null(cv_fold) && (cv_fold %% 1 != 0 || cv_fold <= 1))
    stop("cv_fold must be an integer larger than 1 if provided.")
  s_output("Calculating elements necessary for estimation.", verbose, verbosetext)
  if (is.null(elts))
    elts <- get_elts(h_hp, x, setting, domain, centered=centered,
                     profiled_if_noncenter = is.infinite(lambda_ratio) && (domain$type != "simplex"), # profiled not supported for models on the simplex
                     scale=scale, diagonal_multiplier=diagonal_multiplier,
                     unif_dist=unif_dist)
  if (is.null(lambda1s)){
    s_output("Calculating lower bound for lambda.", verbose, verbosetext)
    lambda_lo <- test_lambda_bounds2(elts, symmetric, lambda_ratio, lower=TRUE, verbose=verbose, tol=tol, maxit=maxit)
    s_output("Calculating upper bound for lambda.", verbose, verbosetext)
    lambda_hi <- lambda_max(elts = elts, symmetric=symmetric, lambda_ratio = lambda_ratio)
    s_output(paste("Upper bound for lambda: ", lambda_hi, ".", sep=""), verbose, verbosetext)
    if (lambda_hi > 1e10) # If unreasonable lambda_hi detected, use search; 1e10 chosen since it is smaller than 1e15 in test_lambda_bounds()
      lambda_hi <- test_lambda_bounds2(elts, symmetric, lambda_ratio, lower=FALSE, verbose=verbose, tol=tol, maxit=maxit, lambda_start = 1)
    ## In the non-profiled non-centered case, the analytic solution for profiled/centered does not hold; use that result as a starting point to search
    if ((!elts$centered && !elts$profiled_if_noncenter) || symmetric == "and")
      lambda_hi <- test_lambda_bounds2(elts, symmetric, lambda_ratio, lower=FALSE, verbose=verbose, tol=tol, maxit=maxit, lambda_start = lambda_hi)
    lambda1s <- exp(seq(log(lambda_hi), log(lambda_lo), length.out=lambda_length))
    if (warmstart)
      res <- get_results(elts, symmetric, lambda_hi, lambda2=lambda_hi/lambda_ratio, tol=tol, maxit=maxit, previous_res = NULL, is_refit=FALSE)
    res <- get_results(elts, symmetric, lambda_lo, lambda2=lambda_hi/lambda_ratio, tol=tol, maxit=maxit, previous_res = NULL, is_refit=FALSE)
  } else {
    lambda_length <- length(lambda1s)
    lambda1s <- sort(lambda1s, decreasing = TRUE) # Sort lambdas
    upper_lambda_calculated <- FALSE
    res <- NULL
  }
  edgess <- list()
  BICs <- matrix(Inf, ncol=2*length(eBIC_gammas), nrow=lambda_length)
  if (!elts$centered) {etas <- matrix(NA, nrow=lambda_length, ncol=p)
  } else {etas <- NULL}
  convergeds <- numeric(lambda_length)
  iters <- numeric(lambda_length)
  s_output("Calculating estimates.", verbose, verbosetext)
  if (verbose)
    checkpoints <- ceiling(c(0.1, 0.2, 0.5, 1:10) * lambda_length / 10)
  raw_estimates <- list()
  for (lambda_index in 1:lambda_length){
    if (!warmstart)
      res <- NULL
    res <- get_results(elts, symmetric, lambda1=lambda1s[lambda_index], lambda2=lambda1s[lambda_index]/lambda_ratio, tol=tol, maxit=maxit, previous_res=res, is_refit=FALSE)
    convergeds[lambda_index] <- res$converged
    iters[lambda_index] <- res$iters
    edgess[[lambda_index]] <- res$edges
    if (return_raw) {raw_estimates[[lambda_index]] <- res$K}
    if (!is.null(res$eta)){
      etas[lambda_index, ] <- res$eta;# exclude_eta <- res$exclude_eta
    }
    BICs[lambda_index, ] <- eBIC(res, elts, BIC_refit=BIC_refit, gammas=eBIC_gammas)
    #if (BIC_refit && any(is.infinite(BICs[lambda_index, (length(eBIC_gammas)+1):(2*length(eBIC_gammas))]))) ## If asked for BIC with refit but the BIC is infinity
    #  warning("Design sub-matrix not invertible, returning Inf eBIC.\n")
    if (verbose && lambda_length > 10 && lambda_index %in% checkpoints) {
      pretext <- ifelse (is.null(cv_fold), "", "Entire dataset: ")
      s_output(paste(pretext, floor(lambda_index/lambda_length*100), "% done", sep=""), verbose, verbosetext)
    }
    if (!return_raw && length(res$edges) == elts$p*(elts$p-1) && lambda_index < lambda_length){ ## If all edges are selected for some lambda, end early
      BICs[(lambda_index+1):lambda_length, ] <- Inf
      convergeds[(lambda_index+1):lambda_length] <- 0#convergeds[lambda_index]
      iters[(lambda_index+1):lambda_length] <- 0
      for (li in (lambda_index+1):lambda_length)
        edgess[[li]] <- edgess[[li-1]]
      break
    }
  }
  lambda_index_stopped <- lambda_index # In case ended early
  if (!is.null(cv_fold)) {
    ids <- make_folds(n, cv_fold, cv_fold_seed)
    cv_losses <- matrix(Inf, nrow=lambda_length, ncol=cv_fold)
    res <- NULL
    for (fold in 1:cv_fold) {
      this_ids <- ids[[fold]]; rest_ids <- Reduce("c", (ids[c(-fold)]))
      elts_this <- get_elts(h_hp, x[this_ids,], setting, domain, centered=centered,
                            profiled_if_noncenter = FALSE, scale=scale, diagonal_multiplier=1,
                            unif_dist=unif_dist)
      elts_rest <- get_elts(h_hp, x[rest_ids,], setting, domain, centered=centered,
                            profiled_if_noncenter = is.infinite(lambda_ratio) && (domain$type != "simplex"), # profiled not supported for models on the simplex
                            scale=scale, diagonal_multiplier=diagonal_multiplier, unif_dist=unif_dist)
      for (lambda_index in 1:lambda_index_stopped){
        if (!warmstart)
          res <- NULL
        res <- get_results(elts_rest, symmetric, lambda1=lambda1s[lambda_index], lambda2=lambda1s[lambda_index]/lambda_ratio, tol=tol, maxit=maxit, previous_res=res, is_refit=FALSE)
        cv_losses[lambda_index, fold] <- calc_crit(elts_this, res, penalty=FALSE)
        if (verbose && lambda_index_stopped > 10 && lambda_index %in% checkpoints)
          s_output(paste("CV fold ", fold, ": ", floor(lambda_index/lambda_index_stopped*100), "% done", sep=""), verbose, verbosetext)
      }
    }
  }

  s_output("Done.", verbose, verbosetext)
  return (list("edgess"=edgess,
               "etas"=switch(changed_from_nc_to_c+1, etas, matrix(0, nrow=lambda_length, ncol=p)),
               "BICs"=BICs[,1:length(eBIC_gammas)],
               "BIC_refits"=switch(BIC_refit+1, NULL, BICs[,(length(eBIC_gammas)+1):(2*length(eBIC_gammas))]),
               "lambda1s"=lambda1s,
               "lambda2s"=switch(1+elts$centered+2*(!elts$centered && elts$profiled), lambda1s/lambda_ratio, NULL, rep(0, length(lambda1s))),
               "cv_losses"=switch(1+is.null(cv_fold), cv_losses, NULL),
               "converged"=convergeds, "iters"=iters,
               "raw_estimates"=raw_estimates, "symmetric"=symmetric,
               "elts"=switch(1+return_elts, NULL, elts)))
}

#' Loss for a refitted (restricted) unpenalized model
#'
#' Refits an unpenalized model restricted to the estimated edges, with \code{lambda1=0}, \code{lambda2=0} and \code{diagonal_multiplier=1}. Returns \code{Inf} if no unique solution exists (when the loss is unbounded from below or has infinitely many minimizers).
#'
#' @param res A list of results returned by \code{get_results()}.
#' @param elts A list, elements necessary for calculations returned by \code{get_elts()}.
#' @return A number, the loss of the refitted model.
#' @details Currently the function only returns \code{Inf} when the maximum node degree is >= the sample size, which is a sufficient and necessary condition for inexistence of a unique solution with probability 1 if \code{symmetric != "symmetric"}. In practice it is also a sufficient and necessary condition for most cases and \code{symmetric == "symmetric"}.
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()},
#' #   as the way to call this function (\code{refit()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' res_nc_np <- get_results(elts_gauss_np, symmetric="symmetric",
#'                lambda1=0.35, lambda2=2, previous_res=NULL, is_refit=FALSE)
#' refit(res_nc_np, elts_gauss_np)
#'
#' elts_gauss_p <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                centered=FALSE, profiled=TRUE, diag=dm)
#' res_nc_p <- get_results(elts_gauss_p, symmetric="symmetric",
#'               lambda1=0.35, lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#' refit(res_nc_p, elts_gauss_p)
#'
#' elts_gauss_c <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'             centered=TRUE, diag=dm)
#' res_c <- get_results(elts_gauss_c, symmetric="or", lambda1=0.35,
#'            lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#' refit(res_c, elts_gauss_c)
#'
#'@export
refit <- function(res, elts){
  if (max(table(res$edges %% elts$p), 0) >= elts$n-2) # If max degree > n, (with multiplier 1) Gamma sub-matrix not invertible, so do not refit; n-1 because res$edges does not contain diagonals
    return (Inf)
  if (res$symmetric != "symmetric"){
    return (get_crit_nopenalty(elts, previous_res=res))
  }
  res$K[-res$edges] <- 0
  if (elts$centered || elts$profiled_if_noncenter){
    test <- get_results(elts, res$symmetric, lambda1=0, tol=res$tol, maxit=res$maxit, previous_res=res, is_refit=TRUE)
  }
  else {
    res$eta[-res$eta_support] <- 0
    test <- get_results(elts, res$symmetric, lambda1=0, lambda2=0, tol=res$tol, maxit=res$maxit, previous_res=res, is_refit=TRUE)
  }
  #if (test$converged)
  #  return (test$crit)
  #else # Temporary fix for non-positive definiteness problem for log_log models on simplex not assuming K1=0 !!!!
  #  return (Inf)
  return (test$crit)
}

#' eBIC score with or without refitting.
#'
#' Calculates the eBIC score both with and without refitting an unpenalized model restricted to the estimated support.
#'
#' @param res A list of results returned by get_results().
#' @param elts A list, elements necessary for calculations returned by get_elts().
#' @param BIC_refit A boolean, whether to get the BIC scores by refitting an unpenalized model restricted to the estimated edges, with \code{lambda1=0}, \code{lambda2=0} and \code{diagonal_multiplier=1}. Default to \code{TRUE}.
#' @param gammas Optional. A number of a vector of numbers. The \eqn{\gamma} parameter in eBIC. Default to \code{c(0,0.5,1)}.
#' @return A vector of length \code{2*length(gammas)}. The first \code{length(gammas)} numbers are the eBIC scores without refitting for each \code{gamma} value, and the rest are those with refitting if \code{BIC_refit == TRUE}, or \code{Inf} if \code{BIC_refit == FALSE}.
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()},
#' #   as the way to call this function (\code{eBIC()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' res_nc_np <- get_results(elts_gauss_np, symmetric="symmetric",
#'                lambda1=0.35, lambda2=2, previous_res=NULL,
#'                is_refit=FALSE)
#' eBIC(res_nc_np, elts_gauss_np, BIC_refit=TRUE, gammas=c(0,0.5,1))
#' @export
eBIC <- function(res, elts, BIC_refit=TRUE, gammas=c(0,0.5,1)){
  ## BIC_refit: calculate BIC (with refit) or not
  n <- elts$n; p <- elts$p
  loss <- 2*n*res$crit
  eBIC_refit <- rep(Inf, length(gammas))
  if (elts$centered || elts$profiled_if_noncenter){
    if (res$symmetric == "symmetric"){
      nedges <- length(res$edges)/2
      eBIC_p2 <- nedges*log(n)
      eBIC_p3s <- 2 * gammas * lchoose(p*(p-1)/2, nedges)
    } else{
      nedges <- length(res$edges)
      eBIC_p2 <- nedges*log(n)
      eBIC_p3s <- 2 * gammas * lchoose(p*(p-1), nedges)
    }
  } else{
    if (res$symmetric == "symmetric")
      nedges1 <- length(res$edges)/2
    else
      nedges1 <- length(res$edges)
    nedges2 <- length(res$eta_support)
    eBIC_p2 <- (nedges1+nedges2)*log(n)
    if (res$symmetric == "symmetric")
      eBIC_p3s <- 2 * gammas * (lchoose(p*(p-1)/2, nedges1) + lchoose(p, nedges2))
    else
      eBIC_p3s <- 2 * gammas * (lchoose(p*(p-1), nedges1) + lchoose(p, nedges2))
  }
  eBIC <- loss + eBIC_p2 + eBIC_p3s
  if (BIC_refit){
    loss_refit <- 2*n*refit(res, elts)
    eBIC_refit <- loss_refit + eBIC_p2 + eBIC_p3s ###
  }
  return (c(eBIC, eBIC_refit))
}


#' Calculates penalized or unpenalized loss in K and eta given arbitrary data
#'
#' Calculates penalized or unpenalized loss in K and eta given arbitrary data
#'
#' @param elts An element list returned from \code{get_elts()}. Need not be the same as the elements used to estimate \code{res}, but they must be both centered or both non-centered, and their dimension \code{p} must match. \code{elts} cannot be profiled as this is supposed to be elements for a new data unseen by \code{res}, in which case the loss must be explicitly written in \code{K} and \code{eta} with \code{Gamma} and \code{g} from a new dataset \code{x}.
#' @param res A result list returned from \code{get_results()}. Must be centered if \code{elts} is centered, and must be non-centered otherwise. Can be profiled. \code{res$p} must be equal to \code{elts$p}.
#' @param penalty A boolean, indicates whether the loss should be penalized (using \code{elts$diagonals_with_multiplier}, \code{res$lambda1} and \code{res$lambda2}).
#' @return A number, the loss.
#' @details This function calculates the loss in some estimated \code{K} and \code{eta} given an \code{elts} generated using \code{get_elts()} with a new dataset \code{x}. This is helpful for cross-validation.
#' @examples
#' # In the following examples, all printed numbers should be close to 0.
#' # In practice, \code{res} need not be estimates fit to \code{elts},
#' # but in the examples we use \code{res <- get_results(elts)} just to
#' # demonstrate that the loss this function returns matches that returned
#' # by the C code during estimation using \code{get_results}.
#'
#' n <- 6
#' p <- 3
#' eta <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#'
#' domains <- list(make_domain("R", p=p),
#'                 make_domain("R+", p=p),
#'                 make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3)),
#'                 make_domain("polynomial", p=p,
#'                   ineqs=list(list("expression"="sum(x^2)<=1", nonnegative=FALSE, abs=FALSE))))
#' \donttest{
#' domains <- c(domains,
#'              list(make_domain("polynomial", p=p,
#'                     ineqs=list(list("expression"="sum(x^2)<=1", nonnegative=TRUE, abs=FALSE))),
#'                   make_domain("polynomial", p=p,
#'                     ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                       function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                       abs=FALSE, nonnegative=TRUE))),
#'                   make_domain("simplex", p=p)))
#' }
#' for (domain in domains) {
#'   if (domain$type == "R" ||
#'        (domain$type == "uniform" && any(domain$lefts < 0)) ||
#'        (domain$type == "polynomial" && !domain$ineqs[[1]]$nonnegative))
#'     settings <- c("gaussian")
#'   else if (domain$type == "simplex")
#'     settings <- c("log_log", "log_log_sum0")
#'   else
#'     settings <- c("gaussian", "exp", "gamma", "log_log", "ab_3/4_2/3")
#'
#'   if (domain$type == "simplex")
#'     symms <- c("symmetric")
#'   else
#'     symms <- c("symmetric", "and", "or")
#'
#'   for (setting in settings) {
#'     x <- gen(n, setting=setting, abs=FALSE, eta=eta, K=K, domain=domain,
#'          finite_infinity=100, xinit=NULL, burn_in=1000, thinning=100, verbose=FALSE)
#'     h_hp <- get_h_hp("min_pow", 1, 3)
#'
#'     for (symm in symms) {
#'
#'        # Centered, penalized loss
#'        elts <- get_elts(h_hp, x, setting, domain, centered=TRUE, scale="", diag=dm)
#'        res <- get_results(elts, symm, 0.1)
#'        print(calc_crit(elts, res, penalty=TRUE) - res$crit) # Close to 0
#'
#'        # Non-centered, unpenalized loss
#'        elts_nopen <- get_elts(h_hp, x, setting, domain, centered=TRUE, scale="", diag=1)
#'        res_nopen <- get_results(elts_nopen, symm, 0)
#'        print(calc_crit(elts_nopen, res_nopen, penalty=FALSE) - res_nopen$crit) # Close to 0
#'
#'        # Non-centered, non-profiled, penalized loss
#'        elts_nc_np <- get_elts(h_hp, x, setting, domain, centered=FALSE,
#'          profiled_if_noncenter=FALSE, scale="", diag=dm)
#'        res_nc_np <- get_results(elts_nc_np, symm, lambda1=0.1, lambda2=0.05)
#'        print(calc_crit(elts_nc_np, res_nc_np, penalty=TRUE) - res_nc_np$crit) # Close to 0
#'
#'        # Non-centered, non-profiled, unpenalized loss
#'        elts_nc_np_nopen <- get_elts(h_hp, x, setting, domain, centered=FALSE,
#'          profiled_if_noncenter=FALSE, scale="", diag=1)
#'        res_nc_np_nopen <- get_results(elts_nc_np_nopen, symm, lambda1=0, lambda2=0)
#'        print(calc_crit(elts_nc_np_nopen, res_nc_np_nopen, penalty=FALSE) -
#'          res_nc_np_nopen$crit) # Close to 0
#'
#'        if (domain$type != "simplex") {
#'          # Non-centered, profiled, penalized loss
#'          elts_nc_p <- get_elts(h_hp, x, setting, domain, centered=FALSE,
#'            profiled_if_noncenter=TRUE, scale="", diag=dm)
#'          res_nc_p <- get_results(elts_nc_p, symm, lambda1=0.1)
#'          if (elts_nc_np$setting != setting || elts_nc_np$domain_type != "R")
#'            res_nc_p$crit <- res_nc_p$crit - sum(elts_nc_np$g_eta ^ 2 / elts_nc_np$Gamma_eta) / 2
#'          print(calc_crit(elts_nc_np, res_nc_p, penalty=TRUE) - res_nc_p$crit)  # Close to 0
#'          # Note that the elts argument cannot be profiled, so
#'          # calc_crit(elts_nc_p, res_nc_p, penalty=TRUE) is not allowed
#'
#'          # Non-centered, profiled, unpenalized loss
#'          elts_nc_p_nopen <- get_elts(h_hp, x, setting, domain, centered=FALSE,
#'            profiled_if_noncenter=TRUE, scale="", diag=1)
#'          res_nc_p_nopen <- get_results(elts_nc_p_nopen, symm, lambda1=0)
#'          if (elts_nc_np_nopen$setting != setting || elts_nc_np_nopen$domain_type != "R")
#'            res_nc_p_nopen$crit <- (res_nc_p_nopen$crit -
#'               sum(elts_nc_np_nopen$g_eta ^ 2 / elts_nc_np_nopen$Gamma_eta) / 2)
#'          print(calc_crit(elts_nc_np_nopen, res_nc_p_nopen, penalty=TRUE) -
#'            res_nc_p_nopen$crit) # Close to 0
#'           # Again, calc_crit(elts_nc_p_nopen, res_nc_p, penalty=TRUE) is not allowed
#'        } # if domain$type != "simplex"
#'
#'     } # for symm in symms
#'   } # for setting in settings
#' } # for domain in domains
#' @export
calc_crit <- function(elts, res, penalty) {
  if (!elts$centered && elts$profiled)
    stop("In calc_crit(): elts must not be profiled if noncentered.")
  if (elts$centered != is.null(res$eta))
    stop("elts and res must be both centered or both noncentered.")
  if (elts$p != res$p)
    stop("elts$p and res$p must be equal.")
  if (penalty) {
    if (elts$setting == "gaussian" && elts$domain_type == "R")
      diag(elts$Gamma_K) <- elts$diagonals_with_multiplier
    else
      elts$Gamma_K[(1:(elts$p*elts$p)-1)*elts$p + 1:elts$p] <- elts$diagonals_with_multiplier
  }
  if (elts$setting == "gaussian" && elts$domain_type == "R") {
    crit <- sum(sapply(1:elts$p, function(i){
      crossprod(res$K[,i], elts$Gamma_K) %*% res$K[,i]})) / 2 - sum(diag(res$K))
    if (!elts$centered) {
      crit <- crit + sum(sapply(1:elts$p, function(i){
        sum(res$K[,i] * elts$Gamma_K_eta) * res$eta[i]
      })) + sum(res$eta ^ 2) / 2
    }
  } else {
    crit <- sum(sapply(1:elts$p, function(i){
      crossprod(res$K[,i], elts$Gamma_K[, (i-1)*elts$p+1:elts$p] %*% res$K[,i] / 2 -
                  elts$g_K[(i-1)*elts$p+1:elts$p])
    }))
    if (!elts$centered) {
      crit <- crit + sum(sapply(1:elts$p, function(i){
        crossprod(res$K[,i], elts$Gamma_K_eta[,i]) * res$eta[i]
      })) - sum(res$eta * elts$g_eta) + sum(res$eta ^ 2 * elts$Gamma_eta) / 2
    }
    if (elts$domain_type == "simplex") {
      crit <- crit + sum(sapply(1:(elts$p-1), function(i){
        res$K[,i] %*% elts$Gamma_K_jp[, (i-1)*elts$p+1:elts$p] %*% res$K[,elts$p]
      }))
      if (!elts$centered) {
        crit <- crit + sum(sapply(1:(elts$p-1), function(i){
          res$K[,i] %*% elts$Gamma_Kj_etap[,i] * res$eta[elts$p] +
            res$K[elts$p,] %*% elts$Gamma_Kp_etaj[,i] * res$eta[i] +
            res$eta[i] * elts$Gamma_eta_jp[i] * res$eta[elts$p]
        }))
      }
    }
  }
  if (penalty) {
    crit <- crit + res$lambda1 * sum(abs(res$K[diag(elts$p) == 0]))
    if (!is.null(res$lambda2)) # res$lambda2 is NULL if res is centered or profiled
      crit <- crit + res$lambda2 * sum(abs(res$eta))
    if (elts$domain_type == "simplex") {
      crit <- crit + (elts$p-2) * res$lambda1 * sum(abs(res$K[-elts$p,elts$p])+abs(res$K[elts$p,-elts$p]))
      if (!is.null(res$lambda2))
        crit <- crit + (elts$p-2) * res$lambda2 * abs(res$eta[elts$p])
    }
  }
  return (crit)
}


#' Minimized loss for unpenalized restricted asymmetric models.
#'
#' Analytic solution of the minimized loss for an unpenalized asymmetric model restricted to a given support. Does not work if \code{symmetric == "symmetric"}.
#'
#' @param elts A list, elements necessary for calculations returned by get_elts().
#' @param exclude Optional. A p*p binary matrix or a p^2 binary vector, where \code{1} indicates the entry in K was estimated to 0 in the previous estimate. Default to \code{NULL}.
#' @param exclude_eta Optional. A p-binary vector, similar to \code{exclude}. Default to \code{NULL}.
#' @param previous_res Optional. A list, the returned list by \code{get_results()} run previously with another lambda value. Default to \code{NULL}.
#' @return A number, the refitted loss.
#' @details If \code{previous_res} is provided, \code{exclude} and \code{exclude_eta} must be \code{NULL} or be consistent with the estimated support in \code{previous_res}. If \code{previous_res} and \code{exclude} are both \code{NULL}, assume all edges are present. The same applies to the non-profiled non-centered case when \code{previous_res} and \code{exclude_eta} are both \code{NULL}.
#' @examples
#' # Examples are shown for Gaussian truncated to R+^p only. For other distributions
#' #   on other types of domains, please refer to \code{gen()} or \code{get_elts()}, as the
#' #   way to call this function (\code{get_crit_nopenalty()}) is exactly the same in those cases.
#' n <- 50
#' p <- 30
#' domain <- make_domain("R+", p=p)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_gauss_np <- get_elts(h_hp, x, setting="gaussian", domain=domain,
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' res_nc_np <- get_results(elts_gauss_np, symmetric="symmetric", lambda1=0.35,
#'                lambda2=2, previous_res=NULL, is_refit=FALSE)
#' get_crit_nopenalty(elts_gauss_np, previous_res=res_nc_np)
#' @export
get_crit_nopenalty <- function(elts, exclude=NULL, exclude_eta=NULL, previous_res=NULL){
  if (elts$setting %in% c("log_log", "log_log_sum0") || elts$domain_type == "simplex")
    stop("get_crit_nopenalty() not supported for log_log and log_log_sum0 settings and simplex domains.")
  tryCatch({
    if (!is.null(previous_res)){
      exclude_given <- !is.null(exclude)
      if (exclude_given) exclude_old <- exclude
      exclude <- matrix(0, elts$p, elts$p)
      exclude <- 1 - diag(elts$p)
      exclude[previous_res$edges] <- 0
      if (!elts$centered && !elts$profiled_if_noncenter){
        exclude_eta_given <- !is.null(exclude_eta)
        if (exclude_eta_given) exclude_eta_old <- exclude_eta
        exclude_eta <- rep(1, elts$p)
        exclude_eta[previous_res$eta_support] <- 0
        if (exclude_eta_given && !elts$centered && !elts$profiled_if_noncenter && any(exclude_eta_old != exclude_eta))
          stop("previous_res inconsistent with exclude_eta provided. Stopped.")
      }
      if (exclude_given && any(exclude != exclude_old))
        stop("previous_res inconsistent with exclude provided. Stopped.")
    }
    if (!is.null(exclude)) diag(exclude) <- 0
    if (elts$setting != "gaussian" || elts$domain_type != "R"){
      sum(sapply(1:elts$p, function(j){
        this_j_edges <- switch(is.null(exclude)+1, which(exclude[,j]==0), 1:elts$p)
        GamKj <- elts$Gamma_K[this_j_edges, (j-1)*elts$p+this_j_edges]
        gKj <- elts$g_K[(j-1)*elts$p+this_j_edges]
        if (!elts$centered && !elts$profiled_if_noncenter && (is.null(exclude_eta) || exclude_eta[j] == 0)){
          gKj <- gKj-elts$Gamma_K_eta[this_j_edges,j]*(elts$g_eta[j]/elts$Gamma_eta[j])
          GamKj <- GamKj-tcrossprod(elts$Gamma_K_eta[this_j_edges,j]/sqrt(elts$Gamma_eta[j]))
          return (-gKj%*%solve(GamKj,gKj)/2-elts$g_eta[j]^2/elts$Gamma_eta[j]/2)
        } else {
          -gKj%*%solve(GamKj,gKj)/2
        }
      }))
    } else {
      sum(sapply(1:elts$p, function(j){
        this_j_edges <- switch(is.null(exclude)+1, which(exclude[,j]==0), 1:elts$p)
        GamKj <- elts$Gamma_K[this_j_edges, this_j_edges]
        ind_j <- match(j, this_j_edges) ## index of j in this_j_edges
        gKj <- replace(numeric(length(this_j_edges)), ind_j, 1)
        if (!elts$centered && !elts$profiled_if_noncenter && (is.null(exclude_eta) || exclude_eta[j] == 0)){
          GamKj <- GamKj-tcrossprod(elts$Gamma_K_eta[this_j_edges])
        }
        return (-solve(GamKj,gKj)[ind_j]/2)
      }))
    }
  }, error=function(s){
    if (s$call == "solve.default(GamKj, gKj)"){
      return (Inf) ## Actually the critical value should be -Inf, but for the purpose of getting a refitted BIC, we return Inf.
    } else {stop(s)}
  })
}


############### 1-D ###############
#' The Cram\'er-Rao lower bound (times \code{n}) for estimating the mean parameter from a univariate truncated normal sample with known variance parameter.
#'
#' The Cram\'er-Rao lower bound (times \code{n}) on the variance for estimating the mean parameter \code{mu} from a univariate truncated normal sample, assuming the true variance parameter \code{sigmasq} is known.
#'
#' @param mu The mean parameter.
#' @param sigmasq The variance parameter.
#' @return A number, the Cram\'er-Rao lower bound.
#' @details The Cram\'er-Rao lower bound in this case is defined as \eqn{\sigma^4/var(X-\mu)}.
#' @export
crbound_mu <- function(mu,sigmasq){
  sigma <- sqrt(sigmasq)
  #return (-1/(exp(-mu^2/sigma^2)*(2*(exp(mu^2/2/sigma^2)*mu*sqrt(2*pi)+sigma-2*exp(mu^2/sigma^2)*pi*sigma)+(-exp(mu^2/2/sigma^2)*mu*sqrt(2*pi)+4*exp(mu^2/sigma^2)*pi*sigma)*(2-2*stats::pnorm(mu/sigma))-exp(mu^2/sigma^2)*pi*sigma*(2-2*stats::pnorm(mu/sigma))^2)/(pi*sigma^3*(2*stats::pnorm(mu/sigma,0,1))^2)))
  EZc2 <- -exp(-mu^2/2/sigma^2)*mu*sigma/sqrt(2*pi)+sigma^2/2+sigma^2*(1-stats::pnorm(0,mu,sigma)*2)/2
  E1 <- 1-stats::pnorm(0,mu,sigma)
  EZc <- exp(-mu^2/2/sigma^2)*sigma/sqrt(2*pi)
  return (sigma^4 / (EZc2/E1 - (EZc/E1)^2))
}

#' The Cram\'er-Rao lower bound (times \code{n}) for estimating the variance parameter from a univariate truncated normal sample with known mean parameter.
#'
#' The Cram\'er-Rao lower bound (times \code{n}) on the variance for estimating the variance parameter \code{sigmasq} from a univariate truncated normal sample, assuming the true mean parameter \code{mu} is known.
#'
#' @name crbound_sigma
#' @param mu The mean parameter.
#' @param sigmasq The variance parameter.
#' @return A number, the Cram\'er-Rao lower bound .
#' @details The Cram\'er-Rao lower bound in this case is defined as \eqn{4\sigma^8/var((X-\mu)^2)}.
#' @export
crbound_sigma <- function(mu,sigmasq){
  sigma <- sqrt(sigmasq)
  EZc2 <- -exp(-mu^2/2/sigma^2)*mu*sigma/sqrt(2*pi)+sigma^2/2+sigma^2*(1-stats::pnorm(0,mu,sigma)*2)/2
  EZc4 <- sigma/2*(-exp(-mu^2/2/sigma^2)*(mu^3*sqrt(2/pi)+3*mu*sqrt(2/pi)*sigma^2)+3*sigma^3*(2-stats::pnorm(0,mu,sigma)*2))
  E1 <- 1-stats::pnorm(0,mu,sigma)
  return (4*sigma^8 / (EZc4/E1 - (EZc2/E1)^2))
}

#' Estimates the mu and sigma squared parameters from a univariate truncated normal sample.
#'
#' Estimates the mu and sigma squared parameters from a univariate truncated normal sample.
#'
#' @param x A vector, the data.
#' @param mode A string, the class of the \code{h} function.
#' @param param1 A number, the first parameter to the \code{h} function.
#' @param param2 A number, the second parameter (may be optional depending on \code{mode}) to the \code{h} function.
#' @param mu A number, may be \code{NULL}. If \code{NULL}, an estimate will be given; otherwise, the value will be treated as the known true \code{mu} parameter and is used to calculate an estimate for \code{sigmasq}, if \code{sigmasq} is \code{NULL}.
#' @param sigmasq A number, may be \code{NULL}. If \code{NULL}, an estimate will be given; otherwise, the value will be treated as the known true \code{sigmasq} parameter and is used to calculate an estimate for \code{mu}, if \code{mu} is \code{NULL}.
#' @return A vector that contains the \code{mu} and the \code{sigmasq} estimates.
#' @details If both \code{mu} and \code{sigmasq} are provided, they are returned immediately. If neither is provided, the estimates are given as \deqn{[1/\sigma^2,\mu/\sigma^2]=\left\{\sum_{i=1}^nh(X_i)[X_i,-1][X_i,-1]^{\top}\right\}^{-1}\left\{\sum_{i=1}^n\left[h(X_i)+h'(X_i)X_i,-h'(X_i)\right]\right\}.}{[1/\sigma^2,\mu/\sigma^2]=(sum(h(Xi)*[Xi,-1][Xi,-1]'))^(-1) (sum([h(Xi)+h'(Xi)Xi, -h'(Xi)])).} If only \code{sigmasq} is provided, the estimate for \code{mu} is given as \deqn{\sum_{i=1}^n[h(X_i)X_i-\sigma^2 h'(X_i)]/\sum_{i=1}^nh(X_i).}{sum(h(Xi)Xi-\sigma^2*h'(Xi))/sum(h(Xi)).} If only \code{mu} is given, the estimate for \code{sigmasq} is given as \deqn{\sum_{i=1}^n h(X_i)(X_i-\mu)^2/\sum_{i=1}^n[h(X_i)+h'(X_i)(X_i-\mu)].}{sum(h(Xi)(Xi-\mu)^2)/(sum(h(Xi)+h'(Xi)(Xi-\mu))).}
#' @export
mu_sigmasqhat <- function(x, mode, param1, param2, mu=NULL, sigmasq=NULL){
  if (!is.null(mu) && !is.null(sigmasq)){
    return (c(mu, sigmasq))
  }
  hx_hpx <- get_h_hp(mode,param1,param2)(x)
  hx <- hx_hpx$hx
  hpx <- hx_hpx$hpx
  if (!is.null(mu)){
    return (c(mu, sum(hx*(x-mu)^2)/sum(hx+hpx*(x-mu))))
  } else if (!is.null(sigmasq)) {
    return (c(-sigmasq*(sum(-hx*x/sigmasq+hpx))/sum(hx), sigmasq))
  }
  hxx <- hx*x
  mean_hx_x <- mean(hxx)
  mean_hx_x2 <- mean(hxx*x)
  mean_hx <- mean(hx)
  mean_hpx_x <- mean(hpx*x)
  mean_hpx <- mean(hpx)
  k_eta <- solve(matrix(c(mean_hx_x2, -mean_hx_x, -mean_hx_x, mean_hx), nrow=2, ncol=2),
                 c(mean_hx+mean_hpx_x, -mean_hpx))
  return (c(k_eta[2]/k_eta[1], 1/k_eta[1]))
}

#' Asymptotic log of \code{h} and \code{hp} functions for large \code{x} for some modes.
#'
#' Asymptotic log of \code{h} and \code{hp} functions for large \code{x} for some modes.
#'
#' @param mode A string, the class of the \code{h} function.
#' @param para A number, the first parameter to the \code{h} function.
#' @return A list of two functions, \code{logh} and \code{loghp}.
#' @export
get_safe_log_h_hp <- function(mode, para){
  # works for e.g. when para*a >= 100, except for log_pow and pow which are exact for all para and a
  if (mode == "asinh") return (list(logh=function(a){log(log(para*a*2))},loghp=function(a){-log(a)}))
  else if (mode == "cosh") return (list(logh=function(a){para*a-log(2)},loghp=function(a){log(para)+para*a-log(2)}))
  else if (mode == "exp") {return (list(logh=function(a){para*a},loghp=function(a){log(para)+para*a}))}
  else if (mode == "identity") {return (list(logh=function(a){log(a)},loghp=function(a){numeric(length(a))}))}
  else if (mode == "log_pow") {return (list(logh=function(a){para*log(log(1+a))},loghp=function(a){log(para)+log(log(1+a))*(para-1)-log(1+a)}))}
  else if (mode == "pow") {return (list(logh=function(a){para*log(a)},loghp=function(a){log(para)+log(a)*(para-1)}))}
  else if (mode == "sinh") return (list(logh=function(a){para*a-log(2)},loghp=function(a){log(para)+para*a-log(2)}))
  else if (mode == "softplus") return (list(logh=function(a){log(para*a-log(2))},loghp=function(a){log(para)}))
  else if (mode == "tanh") {return (list(logh=function(a){rep(1,length(a))},loghp=function(a){numeric(length(a))}))}
  else {warning("not supported"); return(list(logh=NULL, loghp=NULL))}
}

#' The truncation point for \code{h}.
#'
#' The truncation point for \code{h}.
#'
#' @param mode A string, the class of the \code{h} function.
#' @param param1 A number, the first parameter to the \code{h} function.
#' @param param2 A number, the second parameter (may be optional depending on \code{mode}) to the \code{h} function.
#'
#' @return Returns the truncation point (the point where \code{h} becomes constant and \code{hp} becomes 0) for some selected modes.
#' @export
get_trun <- function(mode, param1, param2){
  trun <- Inf
  if (mode %in% c("mcp", "scad")) {trun <- param1*param2
  } else if (mode == "min_asinh") {trun <- sinh(param2)/param1
  } else if (mode == "min_cosh") {trun <- acosh(param2+1)/param1
  } else if (mode == "min_exp") {trun <- log(param2+1)/param1
  } else if (mode == "min_log_pow") {trun <- exp(param2)-1
  } else if (mode == "min_pow") {trun <- param2
  } else if (mode == "min_sinh") {trun <- asinh(param2)/param1
  } else if (mode == "min_softplus") {trun <- log(exp(param2)*2-1)/param1
  } else if (mode == "truncated_sin") {trun <- pi/2/param1
  } else if (mode == "truncated_tan") {trun <- pi/4/param1}
  return (trun)
}

#' Asymptotic variance (times \code{n}) of the estimator for \code{mu} or \code{sigmasq} for the univariate truncated normal assuming the other parameter is known.
#'
#' Asymptotic variance (times \code{n}) of the estimator for \code{mu} or \code{sigmasq} for the univariate truncated normal assuming the other parameter is known.
#'
#' @param mu A number, the true \code{mu} parameter.
#' @param sigmasq A number, the true \code{sigmasq} parameter.
#' @param mode A string, the class of the \code{h} function.
#' @param param1 A number, the first parameter to the \code{h} function.
#' @param param2 A number, the second parameter (may be optional depending on \code{mode}) to the \code{h} function.
#' @param est_mu A boolean. If \code{TRUE}, returns the asymptotic variance of muhat assuming sigmasq is known; if \code{FALSE}, returns the asymptotic variance of sigmasqhat assuming mu is known.
#' @return A number, the asymptotic variance.
#' @details The estimates may be off from the empirical variance, or may even be \code{Inf} or \code{NaN} if \code{"mode"} is one of \code{"cosh"}, \code{"exp"}, and \code{"sinh")} as the functions grow too fast.
#' If \code{est_mu == TRUE}, the function numerically calculates
#' \deqn{E\left[\sigma^2 h^2(X)+\sigma^4 {h'}^2(X)\right]/E^2[h(X)],}{E[\sigma^2h(X)^2+\sigma^4hp(X)^2]/E[h(X)]^2,}
#' and if \code{est_mu == FALSE}, the function numerically calculates
#' \deqn{E\left[\left(2\sigma^6h^2(X)+\sigma^8{h'}^2(X)\right)(X-\mu)^2\right]/E^2\left[h(X)(X-\mu)^2\right],}{E[(2\sigma^6h(X)^2+\sigma^8hp(X)^2)(X-\mu)^2]/E[h(X)(X-\mu)^2]^2,}
#' where \eqn{E} is the expectation over the true distribution \eqn{TN(\mu,\sigma)} of \eqn{X}.
#' @examples
#' varhat(0, 1, "min_log_pow", 1, 1, TRUE)
#' varhat(0.5, 4, "min_pow", 1, 1, TRUE)
#' @export
varhat <- function(mu, sigmasq, mode, param1, param2, est_mu){
  if (!requireNamespace("cubature", quietly = TRUE))
    stop("Please install package \"cubature\".")
  sigma <- sqrt(sigmasq)
  inte <- function(f,from=0,to=Inf){stats::integrate(f,lower=from,upper=to,rel.tol=1e-10)$value}
  adaptinte <- function(f,from=0,to=Inf){cubature::adaptIntegrate(function(t){x<-t/(1-t);1/(1-t)^2*f(x)},lowerLimit=from/(from+1),upperLimit=ifelse(is.infinite(to),1,to/(to+1)),tol=1e-10)$integral}
  all_inte <- function(f,from=0,to=Inf){tryCatch(inte(f,from=from,to=to), error=function(e){adaptinte(f,from,to)})}
  h_hp <- get_h_hp(mode, param1, param2)
  h <- function(x){h_hp(x)$hx}; hp <- function(x){h_hp(x)$hpx}
  right_limit <- sqrt(-(log(1e-20)+log(sigma*sqrt(2*pi)))*sigma^2*2)+mu # such that dnorm(right_limit, mu, sigma) = 1e-20
  if (!mode %in% c("asinh", "cosh", "exp", "identity", "log_pow", "pow", "sinh", "softplus", "tanh")) {
    trun <- get_trun(mode, param1, param2)
    if (est_mu) {
      if (trun < right_limit) {
        Es2h2xANDs4hp2x <- all_inte(function(x){(sigma^2*h(x)^2+sigma^4*hp(x)^2)*stats::dnorm(x,mu,sigma)},to=trun) + sigma^2*h(trun)^2*(1-stats::pnorm(trun,mu,sigma))
        Eh <- all_inte(function(x){h(x)*stats::dnorm(x,mu,sigma)},to=trun) + h(trun)*(1-stats::pnorm(trun,mu,sigma))
      } else {
        Es2h2xANDs4hp2x <- all_inte(function(x){(sigma^2*h(x)^2+sigma^4*hp(x)^2)*stats::dnorm(x,mu,sigma)},to=right_limit)
        Eh <- all_inte(function(x){h(x)*stats::dnorm(x,mu,sigma)},to=right_limit)
      }
      return (Es2h2xANDs4hp2x / Eh^2 * (1-stats::pnorm(0,mu,sigma)))
    } else {
      if (trun < right_limit) {
        EXmu2_from_trun <- all_inte(function(x){exp(2*log(abs(x-mu))+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun)
        E2s6h2Xmu2ANDs8hp2Xmu2 <- all_inte(function(x){sigma^6*(2*h(x)^2+sigma^2*hp(x)^2)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=trun) + 2*sigma^6*h(trun)^2*EXmu2_from_trun
        EhXmu2 <- all_inte(function(x){h(x)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=trun) + h(trun)*EXmu2_from_trun
      } else {
        E2s6h2Xmu2ANDs8hp2Xmu2 <- all_inte(function(x){sigma^6*(2*h(x)^2+sigma^2*hp(x)^2)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=right_limit)
        EhXmu2 <- all_inte(function(x){h(x)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=right_limit)
      }
      return (E2s6h2Xmu2ANDs8hp2Xmu2/EhXmu2^2 * (1-stats::pnorm(0,mu,sigma)))
    }
  } else{
    trun <- 100/param1
    log_h_hp <- get_safe_log_h_hp(mode, param1); logh <- log_h_hp$logh; loghp <- log_h_hp$loghp
    if (est_mu) {
      if (trun < right_limit) {
        Es2h2xANDs4hp2x_to100 <- all_inte(function(x){(sigma^2*h(x)^2+sigma^4*hp(x)^2)*stats::dnorm(x,mu,sigma)},to=trun)
        Eh2_from100 <- all_inte(function(x){exp(2*logh(x)+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun,to=right_limit)
        Ehp2_from100 <- all_inte(function(x){exp(2*loghp(x)+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun,to=right_limit)
        Eh_to100 <- all_inte(function(x){h(x)*stats::dnorm(x,mu,sigma)},to=trun)
        Eh_from100 <- all_inte(function(x){exp(logh(x)+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun,to=right_limit)
      } else {
        Es2h2xANDs4hp2x_to100 <- all_inte(function(x){(sigma^2*h(x)^2+sigma^4*hp(x)^2)*stats::dnorm(x,mu,sigma)},to=right_limit)
        Eh_to100 <- all_inte(function(x){h(x)*stats::dnorm(x,mu,sigma)},to=right_limit)
        Eh2_from100 <- Ehp2_from100 <- Eh_from100 <- 0
      }
      return ((Es2h2xANDs4hp2x_to100+sigma^2*Eh2_from100+sigma^4*Ehp2_from100) / (Eh_from100+Eh_to100)^2 * (1-stats::pnorm(0,mu,sigma)))
    } else {
      if (trun < right_limit) {
        E2s6h2Xmu2ANDs8hp2Xmu2_to100 <- all_inte(function(x){sigma^6*(2*h(x)^2+sigma^2*hp(x)^2)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=trun)
        Eh2Xmu2_from100 <- all_inte(function(x){exp(2*logh(x)+2*log(x-mu)+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun,to=right_limit)
        EhpXmu2_from100 <- all_inte(function(x){exp(2*loghp(x)+2*log(x-mu)+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun,to=right_limit)
        EhXmu2_to100 <- all_inte(function(x){h(x)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=trun)
        EhXmu2_from100 <- all_inte(function(x){exp(logh(x)+2*log(x-mu)+stats::dnorm(x,mu,sigma,log=TRUE))},from=trun,to=right_limit)
      } else {
        E2s6h2Xmu2ANDs8hp2Xmu2_to100 <- all_inte(function(x){sigma^6*(2*h(x)^2+sigma^2*hp(x)^2)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=right_limit)
        EhXmu2_to100 <- all_inte(function(x){h(x)*(x-mu)^2*stats::dnorm(x,mu,sigma)},to=right_limit)
        Eh2Xmu2_from100 <- EhpXmu2_from100 <- EhXmu2_from100 <- 0
      }
      return ((E2s6h2Xmu2ANDs8hp2Xmu2_to100+2*sigma^6*Eh2Xmu2_from100+sigma^8*EhpXmu2_from100)/(EhXmu2_to100+EhXmu2_from100)^2 * (1-stats::pnorm(0,mu,sigma)))
    }
  }
}



