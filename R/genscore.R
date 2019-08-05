### Now includes implementations for sqrt exp, gamma
### Update in 0928 from 0925: implemented centered and profiled for sqrtexp and gamma
### Update in 1002 from 0928: implemented sampler for general ab; !!!
### For exp, 1002.gibbs_c(...,Theta,...) = 0928.gibbs_c(...,Theta*2,...): IN genscore.so, linear part in def of exp density changed from Theta%*%sqrt(x) to (Theta%*%((x^1/2-1)/(1/2))) as in a/b !!!
### 1011: flipped the sign for K in non-trun_gaussian settings to conform with trun_gaussian
### 310322: Added implementation for asymmetric optimization.
### Change from 310322 in package:
###    1. Renamed Theta as eta and Phi as K.
###    2. The random generators now takes a PSD K as the input and passes -K to C.
###    3. Renamed every as thinning.
###    4. Renamed get_elts_all as get_elts.
###    5. Renamed Sd20, Sd21, Sd22, g1, g2 to their appropriate names.
###    6. In crossprod for Gamma for ab and gamma, changed cbind(sqrt(hx[,j]/n)*.., sqrt(hx[,j]/n)*...) to sqrt(hx[,j]/n)*cbind(..,...)
###    7. Added get_elts_trun_gauss: R implementation of elts for truncated gaussian.
###    8. Renamed rsqrt_gibbs_c to rexp_gamma_reject_R
###    9. Renamed rsqrt_gibbs_ab_c to rab_arms_R
###    10. Renamed all min_..._trun as  min_....
###    11. In get_h_hp, changed the first line to stop execution when mode not given.
###    12. In get_h_hp, removed functions related to ext and slope, and generalized functions min_... with truncation point 1 to truncation point para2 (merged if both functions existed).
###    13. In get_h_hp, changed (and also generalized) the definition of min_exp.
###    14. In get_h_hp, added asinh, min_asinh, softplus and min_softplus.
###    15. Sanity checks added to get_results.
###    16. Added variable upper_lambda_calculated: TRUE if the lambda sequence is calculated and the first lambda ensures an empty graph, so initialize exclude to 1-diag(p); otherwise, initialize exclude to matrix(0,p,p).
###    17. Restructured test_lambda_bound, and removed the want_edges parameter as it can be inferred from \code{lower}. Now correctly returns the result for the best lambda instead of that for the previous lambda as in the previous  version.
###    18. (Also changed in 310322.R) lambda_max accommodated for asymmetric.
###    19. Changed default value of eBIC_gammas in estimate() from c(1) to c(0,0.5,1)
###    20. (Also in 310322.R) estimate(): changed the default return values for lambdas after complete graph
###    21. (Also changed in 310322.R) get_results(): reset previous_res$K[exclude==1] and previous_res$eta[exclude_eta==1] to 0 before calling .C. Otherwise for symmetric=="and" the non-zero values for entries corresponding to absent edges (i.e. Kij!=0 but ij is absent because Kji=0) are ignored and not updated by .C (because they are excluded) but the critical value and the estimates will be wrong (since they are based on the old values).
###    22. Added get_crit_nopenalty() that analytically calculates the refitted loss. This can be faster than estimation using gradient descent and can be more numerically stable: #' n <- 50
###       n <- 20; p <- 30; h_hp <- get_h_hp("min_pow", 1, 3); mu <- rep(0, p); K <- diag(p); diagonal_multiplier <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n)))); x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K), lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs", burn.in.samples = 100, thinning = 10); elts_C <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian", centered=TRUE, profiled_if_noncenter=TRUE, scale="norm", diagonal_multiplier=diagonal_multiplier, use_C=TRUE); previous_res <- get_results(elts, "or", 0.359917, 0, previous_res=NULL, is_refit=FALSE)
###       ## Different:
###       get_results(elts, "or", 0, 0, previous_res=previous_res, is_refit=TRUE)$crit; get_crit_nopenalty(elts, previous_res=previous_res)
###       ## Then look into the details:
###       exclude <- matrix(0, elts$p, elts$p); exclude <- 1 - diag(elts$p); exclude[previous_res$edges] <- 0; sapply(1:elts$p, function(j){this_j_edges <- switch(is.null(exclude)+1, which(exclude[,j]==0), 1:elts$p); GamKj <- elts$Gamma_K[this_j_edges, (j-1)*elts$p+this_j_edges]; gKj <- elts$g_K[(j-1)*elts$p+this_j_edges];  tmp$K[,j][this_j_edges] - solve(GamKj, gKj)})
###
###  In Score_init.R:
###    1. Moved helper functions for the competing methods to helper_funs_other_methods.R
###    2. Renamed Posdef as ran_mat.
###    3. Renamed mode "p15" in cov_cons to "sub".
###    4. In cov_cons: now returns as.matrix(K) to ensure the result (especially for mode == "sub") is a numeric matrix.
###    5. In avgrocs: Fixed error: inside tpr_for_fpr(), if pointers[curve_index] already reaches the end, interpolate between roc[pointers[curve_index], ] and c(1,1), instead of roc[pointers[curve_index]+1, ].



### TO DO: 1. Gaussian, 2. check asymmetric cases: Loss should be in terms of both K%*%mu and t(K)%*%mu instead of K%*%mu due to asymmetric??? Also check loss in the .c file.

require(utils) ## txtProgressBar, setTxtProgressBar
require(cubature) ## adaptIntegrate
#dyn.load("src/Score_gen.so")
#dyn.load("src/sampler.so")

## Package:
#library(knitr); library(rmarkdown); library(devtools); library(roxygen2)
#document(); build(); install(); check()

#' The R implementation to get the elements necessary for calculations for general a and b.
#'
#' The R implementation to get the elements necessary for calculations for general \eqn{a} and \eqn{b}.
#'
#' @param hx A matrix, \eqn{h(\mathbf{x})}{h(x)}, should be of the same dimension as \code{x}.
#' @param hpx A matrix, \eqn{h'(\mathbf{x})}{h'(x)}, should be of the same dimension as \code{x}.
#' @param x A matrix, the data matrix.
#' @param a A number, must be strictly larger than \eqn{b/2}.
#' @param b A number, must be >= 0.
#' @param setting A string that indicates the setting. Returned without being checked or used in the function body.
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' hx <- t(apply(x, 1, h_hp$h))
#' hpx <- t(apply(x, 1, h_hp$hp))
#' get_elts_ab(hx, hpx, x, a=0.5, b=0.5, setting="trun_gaussian",
#'             centered=TRUE, scale="norm", diag=1.5)
#' get_elts_ab(hx, hpx, x, a=0.7, b=1.2, setting="ab_0.7_1.2",
#'             centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#' @export
get_elts_ab <- function(hx, hpx, x, a, b, setting,
                        centered=TRUE, profiled_if_noncenter=TRUE,
                        scale="norm", diagonal_multiplier=1){
  ### centered and profiled_if_noncenter IGNORED, just for compatibility
  if (b < 0 || 2*a <= b) {stop("a and b must satisfy 2a > b >= 0.")}
  n <- dim(x)[1]; p <- dim(x)[2]
  xa <- x^a
  if (a == 1/2){xa_1 <- 1/xa
  } else {xa_1 <- xa/x}
  if (b-a == 1) {xb_1 <- xa
  } else {xb_1 <- x^(b-1)}
  g_K <- c(crossprod(xa, hpx*xa_1+(a-1)*hx*xa_1/x)/n + diag(a*colMeans(hx*xa_1*xa_1)))
  if (centered){
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hx[,j]/n)*xa_1[,j]*xa)}))
    diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting=setting))
  }
  #x2 <- x^2
  #Gamma_K <- Reduce(cbind, lapply(1:p, function(j){sweep(t(xa), 2, hx[,j]*xa[,j]*xa[,j]/x2[,j]/n, "*") %*% xa}))
  #Gamma_K_eta <- crossprod(xa, hx*xa*xb_1/x)/n
  #Gamma_eta <- colMeans(hx*xb_1*xb_1)
  Gamma0 <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hx[,j]/n)*cbind(-xa_1[,j]*xa, xb_1[,j]))}))
  Gamma_K <- Gamma0[-p-1, -c(1:p)*(p+1)]
  Gamma_K_eta <- Gamma0[-p-1, c(1:p)*(p+1)]
  Gamma_eta <- Gamma0[p+1, c(1:p)*(p+1)]
  remove(Gamma0)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  g_eta <- -colMeans(hpx*xb_1+hx*(b-1)*xb_1/x)
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
#' @param hx A matrix, \eqn{h(\mathbf{x})}{h(x)}, should be of the same dimension as \code{x}.
#' @param hpx A matrix, \eqn{h'(\mathbf{x})}{h'(x)}, should be of the same dimension as \code{x}.
#' @param x A matrix, the data matrix.
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' hx <- t(apply(x, 1, h_hp$h))
#' hpx <- t(apply(x, 1, h_hp$hp))
#' get_elts_exp(hx, hpx, x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_exp(hx, hpx, x, centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#' @export
get_elts_exp <- function(hx, hpx, x, centered=TRUE, profiled_if_noncenter=TRUE, scale="norm", diagonal_multiplier=1){
  ## scale only used in the returned list
  ## Note that elts_exp with centered=TRUE is the same as elts_gamma with centered=TRUE
  n <- dim(x)[1]; p <- dim(x)[2]
  xsqrt <- sqrt(x)
  tmp <- (hpx-0.5*hx/x)/xsqrt
  g_K <- c(crossprod(xsqrt, tmp) / n + diag(colMeans(hx/x)/2))
  if (centered){
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hx[,j]/n)/xsqrt[,j]*xsqrt)}))
    diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="exp"))
  }
  xsqrt <- sqrt(x)
  Gamma0 <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hx[,j]/n)/xsqrt[,j]*cbind(-xsqrt,1))}))
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
#' @param hx A matrix, \eqn{h(\mathbf{x})}{h(x)}, should be of the same dimension as \code{x}.
#' @param hpx A matrix, \eqn{h'(\mathbf{x})}{h'(x)}, should be of the same dimension as \code{x}.
#' @param x A matrix, the data matrix.
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' hx <- t(apply(x, 1, h_hp$h))
#' hpx <- t(apply(x, 1, h_hp$hp))
#' get_elts_gamma(hx, hpx, x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_gamma(hx, hpx, x, centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#' @export
get_elts_gamma <- function(hx, hpx, x, centered=TRUE, profiled_if_noncenter=TRUE, scale="norm", diagonal_multiplier=1){
  ## scale only used in the returned list
  ## Note that elts_exp with centered=TRUE is the same as elts_gamma with centered=TRUE
  n <- dim(x)[1]; p <- dim(x)[2]
  xsqrt <- sqrt(x)
  g_K <- c(crossprod(xsqrt, (hpx-0.5*hx/x)/xsqrt)/n + diag(0.5*colMeans(hx/x)))
  if (centered){
    Gamma_K <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hx[,j]/n)/xsqrt[,j]*xsqrt)}))
    diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gamma"))
  }
  Gamma0 <- Reduce(cbind, lapply(1:p, function(j){crossprod(sqrt(hx[,j]/n)*cbind(-xsqrt/xsqrt[,j], 1/x[,j]))}))
  Gamma_K <- Gamma0[-p-1, -c(1:p)*(p+1)]
  Gamma_K_eta <- Gamma0[-p-1, c(1:p)*(p+1)]
  Gamma_eta <- Gamma0[p+1, c(1:p)*(p+1)]
  remove(Gamma0)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  g_eta <- -colMeans((hpx-hx/x)/x)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gamma"))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="gamma"))
}

#' The R implementation to get the elements necessary for calculations for the truncated gaussian setting (a=1, b=1).
#'
#' @param hx A matrix, \eqn{h(\mathbf{x})}{h(x)}, should be of the same dimension as \code{x}.
#' @param hpx A matrix, \eqn{h'(\mathbf{x})}{h'(x)}, should be of the same dimension as \code{x}.
#' @param x A matrix, the data matrix.
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
#'   \item{setting}{The setting \code{"trun_gaussian"}.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{\Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}.}
#'   \item{t1,t2}{Returned in the profiled non-centered setting, where the \eqn{\boldsymbol{\eta}}{\eta} estimate can be retrieved from \eqn{\boldsymbol{t_1}-\boldsymbol{t_2}\hat{\mathbf{K}}}{t1-t2*\hat{K}} after appropriate resizing.}
#' @details For details on the returned values, please refer to \code{get_elts_ab} or \code{get_elts}.
#' @examples
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' hx <- t(apply(x, 1, h_hp$h))
#' hpx <- t(apply(x, 1, h_hp$hp))
#' get_elts_trun_gauss(hx, hpx, x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_trun_gauss(hx, hpx, x, centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#' @export
get_elts_trun_gauss <- function(hx, hpx, x, centered=TRUE, profiled_if_noncenter=TRUE, scale="norm", diagonal_multiplier=1){
  n <- dim(x)[1]; p <- dim(x)[2]
  g_K <- crossprod(hpx, x)/n
  Gamma_K <- do.call("cbind", lapply(1:p, function(k){1/n * crossprod(x, hx[,k]*x)}))
  Gamma_eta <- colMeans(hx)
  diag(g_K) <- diag(g_K) + Gamma_eta
  g_K <- c(g_K)
  diagonals_with_multiplier <- Gamma_K[(1:p^2-1)*p+rep(1:p,p)]*diagonal_multiplier
  if (centered){
    return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="trun_gaussian"))
  }
  Gamma_K_eta <- -crossprod(x, hx)/n
  g_eta <- -colMeans(hpx)
  if (!profiled_if_noncenter)
    return (list("n"=n, "p"=p, "g_K"=g_K, "g_eta"=g_eta, "Gamma_K"=Gamma_K, "Gamma_K_eta"=Gamma_K_eta, "Gamma_eta"=Gamma_eta, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="trun_gaussian"))
  Gamma12Gamma22inv <- sweep(Gamma_K_eta, MARGIN=2, Gamma_eta, `/`)
  subtmp <- do.call("cbind", lapply(1:p, function(k){tcrossprod(Gamma12Gamma22inv[,k], Gamma_K_eta[,k])})) ## Gamma1flat
  Gamma_K <- Gamma_K - subtmp
  diagonals_with_multiplier <- diagonals_with_multiplier - subtmp[(1:p^2-1)*p+rep(1:p,p)]
  g_K <- g_K - sweep(Gamma12Gamma22inv, MARGIN=2, g_eta, `*`)
  return (list("n"=n, "p"=p, "g_K"=g_K, "Gamma_K"=Gamma_K, "t1"=g_eta/Gamma_eta, "t2"=Gamma12Gamma22inv, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=diagonals_with_multiplier, setting="trun_gaussian"))
}

#' The R implementation to get the elements necessary for calculations for the untruncated gaussian setting.
#'
#' @param x A matrix, the data matrix.
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
#' if (!requireNamespace("mvtnorm", quietly = TRUE)){
#'   stop("Please install package \"mvtnorm\" first.", call. = FALSE)
#' }
#' require(mvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' x <- mvtnorm::rmvnorm(n, mean=mu, sigma=solve(K))
#' get_elts_gauss(x, centered=TRUE, scale="norm", diag=1.5)
#' get_elts_gauss(x, centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#' @export
get_elts_gauss <- function(x, centered=TRUE, profiled_if_noncenter=TRUE, scale="norm", diagonal_multiplier=1){
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


#' The function wrapper to get the elements necessary for calculations for all settings.
#'
#' @param h A function, the \eqn{h} function. Must evaluate to 0 at 0. Ignored if \code{elts} is provided.
#' @param hp A function, the derivative of the \eqn{h} function. Must be provided if \code{h} is provided.
#' @param x A matrix, the data matrix.
#' @param setting A string that indicates the setting, must be one of \code{"exp"}, \code{"gamma"}, \code{"gaussian"}, \code{"trun_gaussian"}, or of the form \code{"ab_NUM1_NUM2"}, where \code{NUM1} is the \code{a} value and \code{NUM2} is the \code{b} value.
#' @param centered A boolean, whether in the centered setting(assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param profiled_if_noncenter A boolean, whether in the profiled setting (\eqn{\lambda_{\boldsymbol{\eta}}=0}{\lambda_\eta=0}) if noncentered. Parameter ignored if \code{centered=TRUE}. Default to \code{TRUE}.
#' @param scale A string indicating the scaling method. If contains \code{"sd"}, columns are scaled by standard deviation; if contains \code{"norm"}, columns are scaled by l2 norm; if contains \code{"center"} and \code{setting == "gaussian"}, columns are centered to have mean zero. Default to \code{"norm"}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier.
#' @param use_C Optional. A boolean, use C (\code{TRUE}) or R (\code{FALSE}) functions for computation. Default to \code{TRUE}. Ignored if \code{setting == "gaussian"}.
#' @param tol Optional. A positive number. If \code{setting != "gaussian"}, function stops if any entry if smaller than -tol, and all entries between -tol and 0 are set to tol, for numerical stability and to avoid violating the assumption that \eqn{h(\mathbf{x})>0}{h(x)>0} almost surely.
#' @return A list that contains the elements necessary for estimation.
#'   \item{n}{The sample size.}
#'   \item{p}{The dimension.}
#'   \item{centered}{The centered setting or not. Same as input.}
#'   \item{scale}{The scaling method. Same as input.}
#'   \item{diagonal_multiplier}{The diagonal multiplier. Same as input.}
#'   \item{diagonals_with_multiplier}{A vector that contains the diagonal entries of \eqn{\boldsymbol{\Gamma}}{\Gamma} after applying the multiplier.}
#'   \item{setting}{The setting. Same as input.}
#'   \item{g_K}{The \eqn{\boldsymbol{g}}{g} vector. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\mathbf{K}}{K}. A \eqn{p^2}-vector. Not returned if \code{setting == "gaussian"} since it is just \eqn{diag(p)}.}
#'   \item{Gamma_K}{The \eqn{\boldsymbol{\Gamma}}{Gamma} matrix with no diagonal multiplier. In the non-profiled non-centered setting, this is the \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\mathbf{K}}{K}. A vector of length \eqn{p^2} if \code{setting == "gaussian"} or \eqn{p^3} otherwise.}
#'   \item{g_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{g}}{g} sub-vector corresponding to \eqn{\boldsymbol{\eta}}{\eta}. A \eqn{p}-vector. Not returned if \code{setting == "gaussian"} since it is just \eqn{numeric(p)}.}
#'   \item{Gamma_K_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to interaction between \eqn{\mathbf{K}}{K} and \eqn{\boldsymbol{\eta}}{\eta}. If \code{setting == "gaussian"}, returns a vector of length \eqn{p}, or \eqn{p^2} otherwise.}
#'   \item{Gamma_eta}{Returned in the non-profiled non-centered setting. The \eqn{\boldsymbol{\Gamma}}{\Gamma} sub-matrix corresponding to \eqn{\boldsymbol{\eta}}{\eta}. A \eqn{p}-vector. Not returned if \code{setting == "gaussian"} since it is just \code{rep(1,p)}.}
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
#' if (!requireNamespace("mvtnorm", quietly = TRUE)){
#'   stop("Please install package \"mvtnorm\" first.", call. = FALSE)
#' }
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(mvtnorm)
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' diagonal_multiplier <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'   centered=TRUE, scale="norm", diag=1.5)
#' get_elts(h_hp$h, h_hp$hp, x, setting="ab_0.7_1.2",
#'   centered=FALSE, profiled=FALSE, scale="sd", diag=1.9)
#'
#' x <- mvtnorm::rmvnorm(n, mean=mu, sigma=solve(K))
#' get_elts(NULL, NULL, x, setting="gaussian", centered=FALSE,
#'   profiled=FALSE, scale="center_norm", diag=1.3)
#' @export
#' @useDynLib genscore elts_exp_c elts_exp_np elts_exp_p elts_gamma_np elts_gamma_p elts_c elts_nc_np elts_nc_p elts_ab_c elts_ab_np elts_ab_p
get_elts <- function(h, hp, x, setting, centered=TRUE, profiled_if_noncenter=TRUE, scale="norm", diagonal_multiplier=1, use_C=TRUE, tol=.Machine$double.eps^0.5){
  ## Note that in the result the diagonals of elts$Gamma_K are without multipliers.
  ## The diagonal entries with multipliers are stored in elts$diagonals_with_multiplier
  if (!(setting %in% c("exp", "gamma", "trun_gaussian", "gaussian") || startsWith(setting, "ab_"))){
    stop("\"setting\" parameter must be one of exp, gamma, trun_gaussian, or gaussian, or ab_A_B.")
  }
  n <- dim(x)[1]; p <- dim(x)[2]
  if (tol <= 0) {stop("tol must be >= 0.")}
  if (setting != "gaussian") {
    if (any(x < -tol)) {stop("All entries in x should be >= 0.")}
    if (any(x <= 0)) {
      cat("Entries in x between -tol and 0 are set to tol.")
      x <- pmax(x, tol)
    }
  }
  ### Violates the assumption that h(x)>0 almost surely since each column will have at least one 0
  #if (substr(scale, 1,3) == "min")
  #  x <- t(t(x)-apply(x, 2, min))
  if (grepl("center", scale)){
    if (setting == "gaussian") {x <- scale(x, center=T, scale=F)
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
  if (setting != "gaussian"){
    hx <- t(apply(x, 1, h))
    hpx <- t(apply(x, 1, hp))
  }
  if (setting == "exp"){
    if (use_C){
      if (centered) {
        res <- .C("elts_exp_c", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_exp_np", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else {
        res <- .C("elts_exp_p", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      }
    } else {return (get_elts_exp(hx, hpx, x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier))}
  } else if (setting == "gamma") {
    if (use_C){
      if (centered) {
        res <- .C("elts_exp_c", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_gamma_np", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else {
        res <- .C("elts_gamma_p", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      }
    } else {return (get_elts_gamma(hx, hpx, x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier))}
  } else if (setting == "trun_gaussian") {
    if (use_C){
      if (centered) {
        res <- .C("elts_c", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_nc_np", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else {
        res <- .C("elts_nc_p", nIn=as.integer(n), pIn=as.integer(p), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      }
    } else {return (get_elts_trun_gauss(hx, hpx, x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier))}
  } else if (setting == "gaussian") {
    #use_C ignored
    return (get_elts_gauss(x, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier))
  } else {
    s <- strsplit(setting,"_")[[1]]
    if (s[1] != "ab" || is.na(as.numeric(s[2])) || is.na(as.numeric(s[3])))
      stop("Setting should be of the form ab_Aval_Bval.")
    a <- as.numeric(s[2]); b <- as.numeric(s[3])
    if (use_C){
      if (centered) {
        res <- .C("elts_ab_c", nIn=as.integer(n), pIn=as.integer(p), a=as.double(a), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), Gamma=as.double(numeric(p^3)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "centered"=TRUE, "scale"=scale, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else if (!profiled_if_noncenter) {
        res <- .C("elts_ab_np", nIn=as.integer(n), pIn=as.integer(p), a=as.double(a), b=as.double(b), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "g_eta"=res$g2, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "Gamma_K_eta"=matrix(res$Gamma12,p,p), "Gamma_eta"=res$d, "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=FALSE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      } else {
        res <- .C("elts_ab_p", nIn=as.integer(n), pIn=as.integer(p), a=as.double(a), b=as.double(b), hx=as.double(hx), hpx=as.double(hpx), x=as.double(x), g1=as.double(numeric(p^2)), g2=as.double(numeric(p)), d=as.double(numeric(p)), Gamma=as.double(numeric(p^3)), Gamma12=as.double(numeric(p^2)), diagonal_multiplier=as.double(diagonal_multiplier), diagonals_with_multiplier=as.double(numeric(p^2)), PACKAGE="genscore")
        return (list("n"=n, "p"=p, "g_K"=res$g1, "Gamma_K"=matrix(res$Gamma,nrow=p,ncol=p^2), "t1"=res$g2, "t2"=matrix(res$Gamma12,p,p), "centered"=FALSE, "scale"=scale, "profiled_if_noncenter"=TRUE, "diagonal_multiplier"=diagonal_multiplier, "diagonals_with_multiplier"=res$diagonals_with_multiplier, setting=setting))
      }
    } else {return (get_elts_ab(hx, hpx, x, a, b, setting, centered=centered, profiled_if_noncenter=profiled_if_noncenter, scale=scale, diagonal_multiplier=diagonal_multiplier))}
  }
}

#' Random data generator from exp or gamma graphs.
#'
#' Random data generator from exponential square-root (\eqn{a=0.5}, \eqn{b=0.5}) or gamma (\eqn{a=0.5}, \eqn{b=0}) graphs using basic rejection sampling.
#'
#' @param n A number, number of observations.
#' @param gamm A boolean, generate from gamma (\code{TRUE}) or exponential square-root (\code{FALSE}) graphs.
#' @param eta A vector, the linear part in the distribution.
#' @param K A strictly co-positive square matrix, the interaction matrix.
#' @param seed Optional. A number, the seed for the random generator.
#' @param burn_in Optional. A positive integer, the number of burn-in samples in Gibbs sampling to be discarded.
#' @param thinning Optional. A positive integer, thinning factor in Gibbs sampling. Samples are taken at iteration steps \eqn{\mathrm{burn\_in}+1}{burn_in+1}, \eqn{\mathrm{burn\_in}+1+\mathrm{thinning}}{burn_in+1+thinning}, ..., \eqn{\mathrm{burn\_in}+1+(n-1)*\mathrm{thinning}}{burn_in+1+(n-1)*thinning}. Default to \code{1000}.
#' @param max_iter Optional. A positive integer, the maximum number of proposals for each sample. Default to \code{100000}.
#' @param verbose Optional. A boolean. If \code{TRUE}, prints a progress bar showing the progress. Defaults to \code{TRUE}.
#' @return An \eqn{n\times p}{n*p} matrix of samples, where \eqn{p} is the length of \code{eta}.
#' @details
#' Note: rab_arms_R() with a=0.5 and b=0.5 (exp) or b=0 (gamma) may be a more stable generator.
#'
#' Randomly generates \code{n} samples from the \code{p}-variate exponential square-root or gamma distributions with parameters \code{eta} and \code{K}, where \code{p} is the length of \code{eta} or the dimension of the square matrix \code{K}.
#'
#' The exponential square-root distribution is proportional to
#' \deqn{\exp\left(-\sqrt{\boldsymbol{x}}^{\top}\mathbf{K}\sqrt{\boldsymbol{x}}+2\boldsymbol{\eta}^{\top}\sqrt{\boldsymbol{x}}\right)}{exp(-sqrt(x)'K*sqrt(x)+2*\eta'sqrt(x))}
#' on the non-negative orthant of \eqn{\mathbf{R}^p}{R^p}.
#'
#' The gamma distribution is proportional to
#' \deqn{\exp\left(-\sqrt{\boldsymbol{x}}^{\top}\mathbf{K}\sqrt{\boldsymbol{x}}+\boldsymbol{\eta}^{\top}\log(\boldsymbol{x})\right)}{exp(-sqrt(x)'K*sqrt(x)+\eta'log(x))}
#' on the non-negative orthant of \eqn{\mathbf{R}^p}{R^p}.
#' @examples
#' p <- 30
#' set.seed(1)
#' eta <- stats::rnorm(p)*0.5
#' eta[eta <= -1] <- abs(eta[eta <= -1])
#' K <- cov_cons("er", p, seed = 3, spars = 0.05, eig = 0.1)
#' rexp_gamma_reject_R(n=1000, gamm=TRUE, eta=eta, K=K, seed=1,
#'   burn_in=100, thinning=1000, max_iter=1e+05, verbose=TRUE)
#' rexp_gamma_reject_R(n=1000, gamm=FALSE, eta=eta, K=K, seed=2,
#'   burn_in=200, thinning=2000, max_iter=1e+05, verbose=TRUE)
#' @export
#' @useDynLib genscore rexp_gamma_reject
rexp_gamma_reject_R <- function(n, gamm, eta, K, seed=NULL, burn_in=1000, thinning=1000, max_iter=100000, verbose=TRUE){
  ### If gamm=TRUE, generates from the gamma graphs
  ### If gamm=FALSE, generates from the sqrt exp graphs
  if (verbose && !requireNamespace("utils", quietly = TRUE)){
    stop("Please install package \"utils\" first.", call. = FALSE)
  }
  m <- length(eta)
  if (nrow(K) != m || ncol(K) != m)
    stop("Dimensions of eta and K do not match.")
  if (any(diag(K) <= 0))
    stop("All diagonal entries of K must be positive.")
  if (!is.null(seed)){set.seed(seed)
  } else {seed <- -1}
  xprev <- abs(stats::rnorm(m))
  sqrtx <- sqrt(xprev)
  res <- matrix(0, n, m)
  num_samples <- burn_in+(n-1)*thinning+1
  if (verbose)
    pb <- utils::txtProgressBar(min=0, max=burn_in+1+(n-1)*thinning, style = 3)
  for (i in 1:n){
    #####cat(sqrtx[1:5], " | ", sum(sqrtx), " | ", min(sqrtx),  " | ", min(xprev))
    ## If fails to generate within the max number of iterations specified, error may occur from NA values in sqrtx (c.f. one() in C code where x[i] is initialized to -1 and sqrt(-1)=NA)
    if (i == 1){
      tmp <- .C("rexp_gamma_reject", gamm=as.integer(gamm), xinit=as.double(xprev), sqrtx=as.double(sqrtx),
                steps=as.integer(burn_in+1), m=as.integer(m), Theta=as.double(eta), Phi=as.double(-K),
                max_iter=as.integer(max_iter), seed=as.integer(seed), PACKAGE="genscore")
    } else {
      tmp <- .C("rexp_gamma_reject", gamm=as.integer(gamm), xinit=as.double(xprev), sqrtx=as.double(sqrtx),
                steps=as.integer(thinning), m=as.integer(m), Theta=as.double(eta), Phi=as.double(-K),
                max_iter=as.integer(max_iter), seed=as.integer(seed), PACKAGE="genscore")
    }
    if (tmp$xinit[1] < 0){
      stop("Failed to generate samples within max number of iterations specified.")
    }
    seed <- tmp$seed
    xprev <- tmp$xinit; sqrtx <- tmp$sqrtx
    res[i,] <- xprev
    if (verbose)
      utils::setTxtProgressBar(pb, burn_in+1+(i-1)*thinning)
  }
  return (res)
}


#' Random data generator from general \code{a}-\code{b} distributions.
#'
#' Random data generator from general \code{a}-\code{b} graphs using adaptive rejection metropolis sampling (ARMS).
#'
#' @param n A number, number of observations.
#' @param a A number, must be strictly larger than \eqn{b/2}.
#' @param b A number, must be >= 0.
#' @param eta A vector, the linear part in the distribution.
#' @param K A strictly co-positive square matrix, the interaction matrix.
#' @param seed Optional. A number, the seed for the random generator.
#' @param burn_in Optional. A positive integer, the number of burn-in samples in ARMS to be discarded.
#' @param thinning Optional. A positive integer, thinning factor in ARMS. Samples are taken at iteration steps \eqn{\mathrm{burn\_in}+1}{burn_in+1}, \eqn{\mathrm{burn\_in}+1+\mathrm{thinning}}{burn_in+1+thinning}, ..., \eqn{\mathrm{burn\_in}+1+(n-1)*\mathrm{thinning}}{burn_in+1+(n-1)*thinning}. Default to \code{1000}.
#' @param verbose Optional. A boolean. If \code{TRUE}, prints a progress bar showing the progress. Defaults to \code{TRUE}.
#' @return An \eqn{n\times p}{n*p} matrix of samples, where \eqn{p} is the length of \code{eta}.
#' @details
#' Randomly generates \code{n} samples from the \code{p}-variate \code{a}-\code{b} distributions with parameters \eqn{\boldsymbol{\eta}}{\eta} and \eqn{\mathbf{K}}{K}, where \code{p} is the length of \eqn{\boldsymbol{\eta}}{\eta} or the dimension of the square matrix \eqn{\mathbf{K}}{K}.
#'
#' The \code{a}-\code{b} distribution is proportional to
#' \deqn{\exp\left(-\frac{1}{2a}{\boldsymbol{x}^a}^{\top}\mathbf{K}{\boldsymbol{x}}^a+\boldsymbol{\eta}^{\top}\frac{\boldsymbol{x}^b-\mathbf{1}_p}{b}\right)}{exp(-x^a'Kx^a/(2a)+eta'(x^b-rep(1,p))/b)}
#' on the non-negative orthant of \eqn{\mathbf{R}^p}{R^p}.
#'
#' @examples
#' p <- 30
#' set.seed(1)
#' eta <- stats::rnorm(p)*0.5
#' K <- cov_cons("er", p, seed = 3, spars = 0.05, eig = 0.1)
#' rab_arms_R(n=100, a=1.5, b=2.3, eta=eta, K=K,
#'   seed=1, burn_in=100, thinning=1000, verbose=TRUE)
#' rab_arms_R(n=100, a=1.2, b=0.9, eta=eta, K=K,
#'   seed=2, burn_in=200, thinning=2000, verbose=TRUE)
#' @export
#' @useDynLib genscore rab_arms
rab_arms_R <- function(n, a, b, eta, K, seed=NULL, burn_in=1000, thinning=1000, verbose=TRUE){
  ### If gamm=TRUE, generates from the gamma graphs
  ### If gamm=FALSE, generates from the sqrt exp graphs
  if (verbose && !requireNamespace("utils", quietly = TRUE)){
    stop("Please install package \"utils\" first.", call. = FALSE)
  }
  m <- length(eta)
  if (nrow(K) != m || ncol(K) != m)
    stop("Dimensions of eta and K do not match.")
  if (any(diag(K) <= 0))
    stop("All diagonal entries of K must be positive.")
  if (!is.null(seed)){set.seed(seed)
  } else {seed <- -1}
  xprev <- abs(stats::rnorm(m))
  xa <- xprev ^ a
  if (b == 0) {xb <- log(xprev)
  } else {xb <- xprev ^ b}
  res <- matrix(0, n, m)
  num_samples <- burn_in+(n-1)*thinning+1
  if (verbose)
    pb <- utils::txtProgressBar(min=0, max=burn_in+1+(n-1)*thinning, style = 3)
  for (i in 1:n){
    if (i == 1){
      tmp <- .C("rab_arms", a=as.double(a), b=as.double(b), xinit=as.double(xprev), xa=as.double(xa), xb=as.double(xb),
                steps=as.integer(burn_in+1), m=as.integer(m),
                Theta=as.double(eta), Phi=as.double(-K), seed=as.integer(seed), PACKAGE="genscore")
    } else {
      tmp <- .C("rab_arms", a=as.double(a), b=as.double(b), xinit=as.double(xprev), xa=as.double(xa), xb=as.double(xb),
                steps=as.integer(thinning), m=as.integer(m),
                Theta=as.double(eta), Phi=as.double(-K), seed=as.integer(seed), PACKAGE="genscore")
    }
    seed <- tmp$seed
    xprev <- tmp$xinit; xa <- tmp$xa; xb <- tmp$xb
    res[i,] <- xprev
    if (verbose)
      utils::setTxtProgressBar(pb, burn_in+1+(i-1)*thinning)
  }
  return (res)
}


#' Generator of h and hp functions.
#'
#' Generator of \code{h} and \code{hp} (\eqn{h'}) functions.
#'
#' @param mode A string, see details.
#' @param para May be optional. A number, the first parameter. Default to \code{NULL}.
#' @param para2 May be optional. A number, the second parameter. Default to \code{NULL}.
#' @return A list containing two functions \code{h} and \code{hp} (\eqn{h'}).
#' @details
#' The \code{mode} parameter can be chosen from the options listed below along with the corresponding definitions of \code{h} under appropriate choices of \code{para} and \code{para2} parameters. Unless otherwise noted, \code{para} and \code{para2}, must both be strictly positive if provided, and are set to 1 if not provided. Functions \code{h} and \code{hp} should only be applied to non-negative values \code{x} and this is not enforced or checked by the functions.
#' \describe{
#'     \item{\code{asinh}}{An asinh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\mathrm{asinh}(\mathrm{para}\cdot\boldsymbol{x})=\log\left(\mathrm{para}\cdot\boldsymbol{x}+\sqrt{(\mathrm{para}\cdot\boldsymbol{x})^2+1}\right)}{h(x)=asinh(para*x)=log(para*x+sqrt((para*x)^2+1))}. Unbounded and takes one parameter. Equivalent to \code{min_asinh(x, para, Inf)}.}
#'     \item{\code{cosh}}{A shifted cosh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\cosh(\mathrm{para}\cdot\boldsymbol{x})-1}{h(x)=cosh(para*x)-1}. Unbounded and takes one parameter. Equivalent to \code{min_cosh(x, para, Inf)}.}
#'     \item{\code{exp}}{A shifted exponential function \eqn{\boldsymbol{h}(\boldsymbol{x})=\exp(\mathrm{para}\cdot\boldsymbol{x})-1}{h(x)=exp(para*x)-1}. Unbounded and takes one parameter. Equivalent to \code{min_exp(x, para, Inf)}.}
#'     \item{\code{identity}}{The identity function \eqn{\boldsymbol{h}(\boldsymbol{x})=\boldsymbol{x}}{h(x)=x}. Unbounded and does not take any parameter. Equivalent to \code{pow(x, 1)} or \code{min_pow(x, 1, Inf)}.}
#'     \item{\code{log_pow}}{A power function on a log scale \eqn{\boldsymbol{h}(\boldsymbol{x})=\log(1+\boldsymbol{x})^{\mathrm{para}}}{log(1+x)^para}. Unbounded and takes one parameter. Equivalent to \code{min_log_pow(x, para, Inf)}.}
#'     \item{\code{mcp}}{Treating \eqn{\lambda}=para, \eqn{\gamma}=para2, the step-wise MCP function applied element-wise: \eqn{\lambda x-x^2/(2\gamma)}{\lambdax-x^2/(2\gamma)} if \eqn{x\leq\lambda\gamma}{x<=\lambda\gamma}, or \eqn{\gamma\lambda^2/2} otherwise. Bounded and takes two parameters.}
#'     \item{\code{min_asinh}}{A truncated asinh function applied element-wise: \eqn{\min(\mathrm{asinh}(\mathrm{para}\cdot\boldsymbol{x}),\mathrm{para}_2)}{pmin(asinh(para*x), para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_cosh}}{A truncated shifted cosh function applied element-wise: \eqn{\min(\cosh(\mathrm{para}\cdot\boldsymbol{x})-1,\mathrm{para}_2)}{pmin(cosh(para*x)-1, para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_exp}}{A truncated shifted exponential function applied element-wise: \eqn{\boldsymbol{h}(\boldsymbol{x})=\min(\exp(\mathrm{para}\cdot\boldsymbol{x})-1,\mathrm{para}_2)}{pmin(exp(para*x)-1, para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_log_pow}}{A truncated power on a log scale applied element-wise: \eqn{\boldsymbol{h}(\boldsymbol{x})=\min(\log(1+\boldsymbol{x}),\mathrm{para}_2)^{\mathrm{para}}}{pmin(log(1+x), para2)^para}. Bounded and takes two parameters.}
#'     \item{\code{min_pow}}{A truncated power function applied element-wise: \eqn{\boldsymbol{h}(\boldsymbol{x})=\min(\boldsymbol{x},\mathrm{para}_2)^{\mathrm{para}}}{pmin(x, para2)^para}. Bounded and takes two parameters.}
#'     \item{\code{min_sinh}}{A truncated sinh function applied element-wise: \eqn{\min(\sinh(\mathrm{para}\cdot\boldsymbol{x}),\mathrm{para}_2)}{pmin(sinh(para*x), para2)}. Bounded and takes two parameters.}
#'     \item{\code{min_softplus}}{A truncated shifted softplus function applied element-wise: \eqn{\min(\log(1+\exp(\mathrm{para}\cdot\boldsymbol{x}))-\log(2),\mathrm{para}_2)}{pmin(log(1+exp(para*x))-log(2), para2)}. Bounded and takes two parameters.}
#'     \item{\code{pow}}{A power function \eqn{\boldsymbol{h}(\boldsymbol{x})=\boldsymbol{x}^{\mathrm{para}}}{h(x)=x^para}. Unbounded and takes two parameter. Equivalent to \code{min_pow(x, para, Inf)}.}
#'     \item{\code{scad}}{Treating \eqn{\lambda}=para, \eqn{\gamma}=para2, the step-wise SCAD function applied element-wise: \eqn{\lambda x}{\lambdax} if \eqn{x\leq\lambda}{x<=\lambda}, or \eqn{(2\gamma\lambda x-x^2-\lambda^2)/(2(\gamma-1))}{(2\gamma\lambdax-x^2-\lambda^2)/(2(\gamma-1))} if \eqn{\lambda<x<\gamma\lambda}, or \eqn{\lambda^2(\gamma+1)/2} otherwise. Bounded and takes two parameters, where \code{para2} must be larger than 1, and will be set to 2 by default if not provided.}
#'     \item{\code{sinh}}{A sinh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\sinh(\mathrm{para}\cdot\boldsymbol{x})}{h(x)=sinh(para*x)}. Unbounded and takes one parameter. Equivalent to \code{min_sinh(x, para, Inf)}.}
#'     \item{\code{softplus}}{A shifted softplus function \eqn{\boldsymbol{h}(\boldsymbol{x})=\log(1+\exp(\mathrm{para}\cdot\boldsymbol{x}))-\log(2)}{h(x)=log(1+exp(para*x))-log(2)}. Unbounded and takes one parameter. Equivalent to \code{min_softplus(x, para, Inf)}.}
#'     \item{\code{tanh}}{A tanh function \eqn{\boldsymbol{h}(\boldsymbol{x})=\tanh(\mathrm{para}\cdot\boldsymbol{x})}{h(x)=tanh(para*x)}. Bounded and takes one parameter.}
#'     \item{\code{truncated_sin}}{A truncated sin function applied element-wise: \eqn{\sin(\mathrm{para}\cdot x)}{sin(para*x)} if \eqn{\mathrm{para}\cdot x\leq\pi/2}{para*x<=\pi/2}, or 1 otherwise. Bounded and takes one parameter.}
#'     \item{\code{truncated_tan}}{A truncated tan function applied element-wise: \eqn{\tan(\mathrm{para}\cdot x)}{tan(para*x)} if \eqn{\mathrm{para}\cdot x\leq\pi/4}{para*x<=\pi/4}, or 1 otherwise. Bounded and takes one parameter.}
#'  }
#' @examples
#' get_h_hp("mcp", 2, 4)
#' get_h_hp("min_log_pow", 1, log(1+3))
#' get_h_hp("min_pow", 1, 3)
#' get_h_hp("min_softplus")
#' @export
get_h_hp <- function(mode, para=NULL, para2=NULL){
  if (is.null(mode) || mode == "") stop ("Mode must be chosen from one of the following:
                           asinh, cosh, exp, identity, log_pow, mcp, min_asinh, min_cosh,
                           min_exp, min_log_pow,
                           min_pow, min_sinh, min_softplus, pow, scad,
                           sinh, softplus, tanh, truncated_sin, truncated_tan.")
  number_of_params <- list("asinh"=1, "cosh"=1, "exp"=1, "identity"=0, "log_pow"=1, "mcp"=2,
                           "min_asinh"=2, "min_cosh"=2, "min_exp"=2, "min_log_pow"=2,
                           "min_pow"=2, "min_sinh"=2, "min_softplus"=2, "pow"=1,
                           "scad"=2, "sinh"=1, "softplus"=1, "tanh"=1, "truncated_sin"=1,
                           "truncated_tan"=1)
  if (!mode %in% names(number_of_params)){
    stop("Mode ", mode, " not supported.")
  }
  if (number_of_params[[mode]]){
    if (is.null(para)) {cat("para not provided, default to 1."); para <- 1
    } else if (para <= 0) {stop("para must be strictly positive.")}
  }
  if (number_of_params[[mode]] == 2){
    if (is.null(para2)) {
      if (mode == "scad") {cat("para2 not provided, default to 2."); para2 <- 2}
      else {cat("para2 not provided, default to 1."); para2 <- 1}
    } else if (para2 <= 0) {stop("para2 must be strictly positive.")}
  }
  if (mode == "asinh") return (list(h=function(a){asinh(para*a)},hp=function(a){para/sqrt((para*a)^2+1)}))
  else if (mode == "cosh") return (list(h=function(a){cosh(para*a)-1},hp=function(a){para*sinh(para*a)}))
  else if (mode == "exp") {return (list(h=function(a){exp(para*a)-1},hp=function(a){para*exp(para*a)}))}
  #else if (mode == "extpow") {if(para<=0)stop("para must be > 0"); return (list(h=function(a){abs(a)^para},hp=switch((para==0)+1, function(a){para*sign(a)*abs(a)^(para-1)}, function(a){rep(0,length(a))})))}
  else if (mode == "identity") return (list(h=identity,hp=function(a){rep(1,length(a))}))
  else if (mode == "log_pow") {return (list(h=function(a){(log(1+a))^para},hp=function(a){para*(log(1+a))^(para-1)/(1+a)}))}
  #else if (mode == "log_slope") return (list(h=function(a){para*log(1+a)},hp=function(a){para/(1+a)}))
  else if (mode == "mcp") {lambda<-para; gamma<-para2; return (list(h=function(a){(lambda*a-a^2/2/gamma)*(a<=gamma*lambda)+gamma*lambda^2/2*(a>gamma*lambda)},hp=function(a){pmax(lambda-a/gamma,0)}))} # para,para2 > 0
  else if (mode == "min_asinh") {return (list(h=function(a){pmin(asinh(para*a),para2)},hp=function(a){(para<=sinh(para2)/a)*para/sqrt((para*a)^2+1)}))}
  else if (mode == "min_cosh") {return (list(h=function(a){pmin(cosh(para*a)-1,para2)},hp=function(a){(a<=acosh(para2+1)/para)*para*sinh(para*a)}))}
  else if (mode == "min_exp") {return (list(h=function(a){pmin(exp(para*a)-1,para2)},hp=function(a){(a<=log(para2+1)/para)*para*exp(para*a)}))}
  #else if (mode == "min_extpow") {if(para<=0)stop("para must be > 0"); return (list(h=function(a){pmin(abs(a),para2)^para},hp=function(a){(abs(a)<para2)*para*abs(a)^(para-1)*sign(a)}))}  #else if (mode == "exp") return (list(h=function(a){exp(a)-1},hp=exp))
  #else if (mode == "min_log_pow") return (list(h=function(a){(pmin(log(1+a),1))^para},hp=function(a){(a<exp(1)-1)*para*(log(1+a))^(para-1)/(1+a)}))
  else if (mode == "min_log_pow") {return (list(h=function(a){(pmin(log(1+a),para2))^para},hp=function(a){(a<exp(para2)-1)*para*(log(1+a))^(para-1)/(1+a)}))}
  #else if (mode == "min_log_slope") return (list(h=function(a){para*pmin(log(1+a),1)},hp=function(a){(a<exp(1)-1)*para/(1+a)}))
  #else if (mode == "min_pow") return (list(h=function(a){pmin(a,1)^para},hp=function(a){(a<1)*para*a^(para-1)}))
  else if (mode == "min_pow") {return (list(h=function(a){pmin(a,para2)^para},hp=function(a){(a<para2)*para*a^(para-1)}))}  #else if (mode == "exp") return (list(h=function(a){exp(a)-1},hp=exp))
  else if (mode == "min_sinh") return (list(h=function(a){pmin(sinh(para*a),para2)},hp=function(a){(para<=asinh(para2)/a)*para*cosh(para*a)}))
  #else if (mode == "min_slope") return (list(h=function(a){para*pmin(a,1)},hp=function(a){(a<1)*para}))
  else if (mode == "min_softplus") return (list(h=function(a){pmin(log(1+exp(para*a))-log(2),para2)},hp=function(a){(a<=log(exp(para2)*2-1)/para)*para/(1+exp(-para*a))}))
  else if (mode == "pow") {return (list(h=function(a){a^para},hp=function(a){para*a^(para-1)}))}
  #else if (mode == "quad_exp") return (list(h=function(a){a*(para*2-a)*(a<=para) + (a>para)*para^2*exp(para2*(para-a))}, hp=function(a){(a<=para)*2*(para-a) - para2*(a>para)*para^2*exp(para2*(para-a))}))
  else if (mode == "scad") {if(para2<=1){stop("para2 must be > 1.")};lambda=para;gamma=para2;return (list(h=function(a){(a<=lambda)*lambda*a+(lambda<a&a<gamma*lambda)*(2*gamma*lambda*a-a^2-lambda^2)/(2*(gamma-1))+(a>=gamma*lambda)*lambda^2*(gamma+1)/2},hp=function(a){(a<=lambda)*lambda+(a>lambda)*pmax(gamma*lambda-a,0)/(gamma-1)}))} # para>0, para2>1
  else if (mode == "sinh") return (list(h=function(a){sinh(para*a)},hp=function(a){para*cosh(para*a)}))
  #else if (mode == "slope") return (list(h=function(a){para*a},hp=function(a){para*rep(1,length(a))}))
  else if (mode == "softplus") return (list(h=function(a){log(1+exp(para*a))-log(2)},hp=function(a){para/(1+exp(-para*a))}))
  else if (mode == "tanh") return (list(h=function(a){tanh(para*a)},hp=function(a){para/cosh(para*a)^2})) # naturally bounded
  else if (mode == "truncated_sin") {return (list(h=function(a){sin(para*a)*(para*a<=pi/2)+(para*a>pi/2)},hp=function(a){(a<=pi/2/para)*cos(para*a)*para}))}
  else if (mode == "truncated_tan") {return (list(h=function(a){tan(para*a)*(para*a<=pi/4)+(para*a>pi/4)},hp=function(a){(a<=pi/4/para)/cos(para*a)^2*para}))}
  else {warning("Mode not supported!");
  return (list(h=NULL, hp=NULL))}
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
#'     \item{lambda2}{Same as in the input, and returned only if \code{elts$centered == FALSE} and \code{elts$profiled_if_noncenter == FALSE}.}
#' @examples
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=FALSE, scale="norm", diag=dm)
#' test_nc_np <- get_results(elts_NC_NP, symmetric="symmetric", lambda1=0.35,
#'                 lambda2=2, previous_res=NULL, is_refit=FALSE)
#' test_nc_np2 <- get_results(elts_NC_NP, symmetric="and", lambda1=0.25,
#'                  lambda2=2, previous_res=test_nc_np, is_refit=FALSE)
#'
#' elts_NC_P <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                centered=FALSE, profiled=TRUE, scale="norm", diag=dm)
#' test_nc_p <- get_results(elts_NC_P, symmetric="symmetric",
#'                lambda1=0.35, lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#'
#' elts_C <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                centered=TRUE, scale="norm", diag=dm)
#' test_c <- get_results(elts_C, symmetric="or", lambda1=0.35,
#'                lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#'
#' @export
#' @useDynLib genscore profiled profiled_asymm full full_asymm
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
  if (elts$centered || elts$profiled_if_noncenter){
    if (symmetric == "symmetric") {call_name <- "profiled"
    } else {call_name <- "profiled_asymm"}
    test <- .C(call_name, p = as.integer(elts$p), Gamma_K = as.double(elts$Gamma_K), g_K = as.double(elts$g_K),
               K = as.double(previous_res$K), lambda1 = as.double(lambda1), tol = as.double(tol), maxit = as.integer(maxit),
               iters = as.integer(0), converged = as.integer(0), crit = as.double(0), exclude = as.integer(exclude),
               previous_lambda1 = as.double(previous_res$lambda1), is_refit = as.integer(is_refit),
               diagonals_with_multiplier = as.double(elts$diagonals_with_multiplier), gauss = as.integer(elts$setting=="gaussian"), PACKAGE="genscore")
  } else{
    if (symmetric == "symmetric") {call_name <- "full"}
    else {call_name <- "full_asymm"}
    test <- .C(call_name, p = as.integer(elts$p), Gamma_K = as.double(elts$Gamma_K), Gamma_K_eta = as.double(elts$Gamma_K_eta),
               Gamma_eta = as.double(elts$Gamma_eta), g_K=as.double(elts$g_K), g_eta=as.double(elts$g_eta), K = as.double(previous_res$K),
               eta=as.double(previous_res$eta), lambda1 = as.double(lambda1), lambda2 = as.double(lambda2), tol = as.double(tol),
               maxit = as.integer(maxit), iters = as.integer(0), converged = as.integer(0), crit = as.double(0),
               exclude = as.integer(exclude), exclude_eta = as.integer(exclude_eta),
               previous_lambda1 = as.double(previous_res$lambda1), is_refit = as.integer(is_refit),
               diagonals_with_multiplier = as.double(elts$diagonals_with_multiplier), gauss = as.integer(elts$setting=="gaussian"), PACKAGE="genscore")
    test$eta_support <- which(abs(test$eta) > tol)  ## For refit
    test$Gamma_K_eta <- test$Gamma_eta <- test$exclude_eta <- NULL
  }
  test$K <- matrix(test$K, elts$p, elts$p)
  if (manual_ncnp_to_c) {
    elts$centered <- FALSE
    test$eta <- numeric(elts$p)
    test$eta_support <- integer(0)
    test$lambda2 <- Inf
  }
  if ((!elts$centered) && elts$profiled_if_noncenter){
    if (elts$setting == "gaussian")
      test$eta <- elts$t1 - test$K %*% elts$t2 ## Note: elts$t1 = 0
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' #require(tmvtnorm)
#' #n <- 50
#' #p <- 30
#' #h_hp <- get_h_hp("min_pow", 1, 3)
#' #mu <- rep(0, p)
#' #K <- diag(p)
#' #dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' #x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#' #       lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#' #       burn.in.samples = 100, thinning = 10)
#'
#' #elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#' #                centered=FALSE, profiled=FALSE, diag=dm)
#' #lambda_cur_res <- test_lambda_bounds(elts_NC_NP, "symmetric", lambda=1,
#' #                       lambda_ratio=1, step=1.5, lower=TRUE, cur_res=NULL)
#' #lambda_cur_res2 <- test_lambda_bounds(elts_NC_NP, "symmetric", lambda=1,
#' #                       lambda_ratio=1, step=1.5, lower=FALSE, cur_res=lambda_cur_res$cur_res)
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
      cat("Testing lower bound for lambda:", lambda)
    cur_res <- get_results(elts, symmetric=symmetric, lambda1=lambda, lambda2=lambda/lambda_ratio, tol=tol, maxit=maxit, previous_res=cur_res)
    if (verbose)
      cat(", number of edges=", length(cur_res$edges), "; want ", want_edges, ".\n", sep="")
    if (length(cur_res$edges) == want_edges){ ## If this lambda gives the desired graph (complete/empty), it is necessarily the best lambda so far
      if (is.null(best_lambda)) ## If no previous lambda is good (i.e. the initial lambda is bad), return this, as the next lambda will necessarily give the right graph but will not be tight
        return (list("lambda"=lambda, "cur_res"=cur_res))
      best_lambda <- lambda; best_res <- cur_res  ## Otherwise, just update the best lambda
    } else if (!is.null(best_lambda)) ## If this lambda does not give the desired graph, but we already found some (best) lambda that does, return that lambda; otherwise keep searching
      return (list("lambda"=best_lambda, "cur_res"=best_res))
    if (lambda <= 1e-10 || lambda >= 1e15){ ## If too small or too large, have to stop
      if (verbose)
        cat("Stopped at ", max(1e-10, min(lambda, 1e15)), ".\n", sep = "")
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' test_lambda_bounds2(elts_NC_NP, "symmetric", lambda_ratio=2,
#'      lower=TRUE, lambda_start=NULL)
#' test_lambda_bounds2(elts_NC_NP, "symmetric", lambda_ratio=2,
#'      lower=FALSE, lambda_start=1.0)
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
      cat("step =", step, "\n")
    lambda_cur_res <- test_lambda_bounds(elts, symmetric=symmetric, lambda=lambda_cur_res$lambda, lambda_ratio=lambda_ratio, step=step, lower=lower, verbose=verbose, tol=tol, maxit=maxit, cur_res=lambda_cur_res$cur_res)
    if (lambda_cur_res$lambda <= 1e-10 || lambda_cur_res$lambda >= 1e15) {
      return (max(1e-10, min(lambda_cur_res$lambda, 1e15)))
    }
    step <- step ^ (1/4)
  }
  if (verbose)
    cat("Final: ", lambda_cur_res$lambda, ", ", length(lambda_cur_res$cur_res$edges), " edges.\n", sep="")
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#'
#' # Exact analytic solution for the smallest lambda such that K and eta are both zero,
#' #  but not a tight bound for K only
#' lambda_max(elts_NC_NP, "symmetric", 2)
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_NC_NP, "symmetric", lambda_ratio=2, lower = FALSE,
#'      lambda_start = lambda_max(elts_NC_NP, "symmetric", 2))
#'
#' # Exact analytic solution for the smallest lambda such that K and eta are both zero,
#' #  but not a tight bound for K only
#' lambda_max(elts_NC_NP, "or", 2)
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_NC_NP, "or", lambda_ratio=2, lower = FALSE,
#'      lambda_start = lambda_max(elts_NC_NP, "or", 2))
#'
#' # An upper bound, not tight.
#' lambda_max(elts_NC_NP, "and", 2)
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_NC_NP, "and", lambda_ratio=2, lower = FALSE,
#'      lambda_start = lambda_max(elts_NC_NP, "and", 2))
#'
#'
#' elts_NC_P <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'               centered=FALSE, profiled=TRUE, diag=dm)
#' # Exact analytic solution
#' lambda_max(elts_NC_P, "symmetric")
#' # Numerical solution, should be close to the analytic solution
#' test_lambda_bounds2(elts_NC_P, "symmetric", lambda_ratio=Inf, lower = FALSE,
#'      lambda_start = NULL)
#'
#' # Exact analytic solution
#' lambda_max(elts_NC_P, "or")
#' # Numerical solution, should be close to the analytic solution
#' test_lambda_bounds2(elts_NC_P, "or", lambda_ratio=Inf, lower = FALSE,
#'      lambda_start = NULL)
#'
#' # An upper bound, not tight
#' lambda_max(elts_NC_P, "and")
#' # Use the upper bound as a starting point for numerical search
#' test_lambda_bounds2(elts_NC_P, "and", lambda_ratio=Inf, lower = FALSE,
#'      lambda_start = lambda_max(elts_NC_P, "and"))
#' @export
lambda_max <- function(elts, symmetric, lambda_ratio=Inf){
  p <- elts$p
  if (elts$setting != "gaussian"){
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
  } else {
    diag(elts$Gamma_K) <- elts$diagonals_with_multiplier
    K_init_diag <- 1/diag(elts$Gamma_K)
    if (symmetric == "symmetric"){
      K_init_off_max <- max(sapply(1:(p-1), function(i){max(abs(-elts$Gamma_K[(i+1):p,i]*K_init_diag[i] - elts$Gamma_K[i,(i+1):p]*K_init_diag[(i+1):p]))}))/2
    } else {
      K_init_off_max <- max(sapply(1:p, function(i){max(abs(-elts$Gamma_K[-i,i]*K_init_diag[i]))}))
    }
    if (elts$centered || elts$profiled_if_noncenter)
      return (K_init_off_max+1e-8)  ## If symmetric!="and", this lambda is exact when centered or profiled
    else {
      if (is.infinite(lambda_ratio))
        stop ("elts$profiled_if_noncenter should be TRUE if lambda_ratio=Inf.")
      eta_init <- max(abs(sapply(1:p, function(i){-K_init_diag[i]*elts$Gamma_K_eta[i]})))
      return (max(K_init_off_max, eta_init * lambda_ratio)+1e-8) ## This only serves an upper bound, since it is the max lambda where both K and eta are 0.
    }
  }
}

#' Helper function for outputting if verbose.
#' @param out Text string.
#' @param verbose Boolean.
#' @param verbosetext Text string.
#' @return If \code{verbose == TRUE}, outputs a string that concatenates \code{verbosetext} and \code{out}.
output <- function(out, verbose, verbosetext){
  if (verbose)
    cat(verbosetext, out, "\n")
}

#' The main function for the generalized score-matching estimator for graphical models.
#'
#' The main function for the generalized score-matching estimator for graphical models.
#'
#' @param x A matrix, the data.
#' @param setting A string that indicates the setting, must be one of \code{"exp"}, \code{"gamma"}, \code{"gaussian"}, \code{"trun_gaussian"}, or of the form \code{"ab_NUM1_NUM2"}, where \code{NUM1} is the \code{a} value and \code{NUM2} is the \code{b} value.
#' @param elts A list (optional), elements necessary for calculations returned by get_elts().
#' @param centered A boolean, whether in the centered setting (assume \eqn{\boldsymbol{\mu}=\boldsymbol{\eta}=0}{\mu=\eta=0}) or not. Default to \code{TRUE}.
#' @param symmetric A string. If equals \code{"symmetric"}, estimates the minimizer \eqn{\mathbf{K}}{K} over all symmetric matrices; if \code{"and"} or \code{"or"}, use the "and"/"or" rule to get the support. Default to \code{"symmetric"}.
#' @param scale A string indicating the scaling method. If contains \code{"sd"}, columns are scaled by standard deviation; if contains \code{"norm"}, columns are scaled by l2 norm; if contains \code{"center"} and \code{setting == "gaussian"}, columns are centered to have mean zero. Default to \code{"norm"}.
#' @param lambda1s A vector of lambdas, the penalty parameter for K.
#' @param lambda_length An integer >= 2, the number of lambda1s. Ignored if \code{lambda1s} is provided, otherwise a grid of lambdas is automatically chosen so that the results range from an empty graph to a complete graph. Default to \code{10} if neither \code{lambda1s} nor \code{lambda_length} is provided.
#' @param lambda_ratio A positive number, the fixed ratio between \eqn{\lambda_{\mathbf{K}}}{\lambda_K} and \eqn{\lambda_{\boldsymbol{\eta}}}{\lambda_\eta}, if \eqn{\lambda_{\boldsymbol{\eta}}\neq 0}{\lambda_\eta!=0} (non-profiled) in the non-centered setting.
#' @param mode A string, the class of the \code{h} function. Ignored if \code{elts}, or \code{h} and \code{hp} are provided, or if \code{setting == "gaussian"}.
#' @param param1 A number, the first parameter to the \code{h} function. Ignored if \code{elts}, or \code{h} and \code{hp} are provided, or if \code{setting == "gaussian"}.
#' @param param2 A number, the second parameter (may be optional depending on \code{mode}) to the \code{h} function. Ignored if \code{elts}, or \code{h} and \code{hp} are provided, or if \code{setting == "gaussian"}.
#' @param h A function, the \eqn{h} function. Must evaluate to 0 at 0. Ignored if \code{elts} is provided, or if \code{setting == "gaussian"}.
#' @param hp A function, the derivative of the \eqn{h} function. Must be provided if \code{h} is provided, or if \code{setting == "gaussian"}.
#' @param verbose Optional. A boolean, whether to output intermediate results.
#' @param verbosetext Optional. A string, text to be added to the end of each printout if \code{verbose == TRUE}.
#' @param tol Optional. A number, the tolerance parameter. Default to \code{1e-6}.
#' @param maxit Optional. A positive integer, the maximum number of iterations for each fit. Default to \code{1000}.
#' @param BIC_refit A boolean, whether to get the BIC scores by refitting an unpenalized model restricted to the estimated edges, with \code{lambda1=lambda2=0} and \code{diagonal_multiplier=1}. Default to \code{TRUE}.
#' @param warmstart Optional. A boolean, whether to use the results from a previous (larger) lambda as a warm start for each new lambda. Default to \code{TRUE}.
#' @param diagonal_multiplier A number >= 1, the diagonal multiplier. Optional and ignored if elts is provided. If \code{ncol(x) > ncol(n)}, a value strictly larger than 1 is recommended. Default to \eqn{1+\left(1-\left(1+4e\max\left(6\log p/n, \sqrt{6\log p/n}\right)\right)^{-1}\right)}{1+(1-1/(1+4e*max(6*log(p)/n, sqrt(6*log(p)/n))))}.
#' @param eBIC_gammas Optional. A number of a vector of numbers. The \eqn{\gamma} parameter in eBIC. Default to \code{c(0,0.5,1)}.
#' @param return_raw A boolean, whether to return the raw estimates of \code{K}. Default to \code{FALSE}.
#' @return
#'    \item{edgess}{A list of vectors of integers: indices of the non-zero edges.}
#'    \item{etas}{If applicable, a \code{lambda_length}*\code{p} matrix of \code{eta} estimates with the \eqn{i}-th row corresponding to the \eqn{i}-th \code{lambda1}, and may contain \code{NA}s after the first lambda that gives the complete graph. Otherwise \code{NULL}.}
#'    \item{BICs}{A \code{lambda_length} by \code{length(eBIC_gammas)} matrix of raw eBIC scores (without refitting). May contain \code{Inf}s for rows after the first lambda that gives the complete graph.}
#'    \item{BIC_refits}{\code{NULL} if \code{BIC_refit == FALSE}, otherwise a \code{lambda_length} by \code{length(eBIC_gammas)} matrix of refitted eBIC scores, obtained by refitting unpenalized models restricted to the estimated edges. May contain \code{Inf}s for rows after the first lambda that gives the graph restricted to which an unpenalized model does not have a solution (loss unbounded from below).}
#'    \item{lambda1s}{A vector of numbers of length \code{lambda_length}: the grid of \code{lambda1}s over which the estimates are obtained.}
#'    \item{lambda2s}{A vector of numbers of length \code{lambda_length}: the grid of \code{lambda2}s over which the estimates are obtained, if applicable, otherwise \code{NULL}.}
#'    \item{converged}{A vector of booleans of length \code{lambda_length}: indicators of convergence for each fit. May contain \code{0}s for all lambdas after the first lambda that gives the complete graph.}
#'    \item{iters}{A vector of integers of length \code{lambda_length}: the number of iterations run for each fit. May contain \code{0}s for all lambdas after the first lambda that gives the complete graph.}
#'    \item{raw_estimate}{An empty list if \code{return_raw == FALSE}, otherwise a list that contains \code{lambda_length} estimates for \code{K} of size \code{ncol(x)}*\code{ncol(x)}. May contain \code{ncol(x)}*\code{ncol(x)} matrices of \code{NA}s for all lambdas after the first lambda that gives the complete graph.}
#' @examples
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' mu <- rep(0, p)
#' K <- diag(p)
#' lambda1s <- c(0.01,0.1,0.2,0.3,0.4,0.5)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' ## Centered estimates, no elts or h provided, mode and params provided
#' est1 <- estimate(x, "trun_gaussian", elts=NULL, centered=TRUE,
#'           symmetric="symmetric", lambda1s=lambda1s, mode="min_pow",
#'           param1=1, param2=3, diag=dm, return_raw=TRUE)
#'
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' ## Centered estimates, no elts provided, h provided; equivalent to est1
#' est2 <- estimate(x, "trun_gaussian", elts=NULL, centered=TRUE,
#'           symmetric="symmetric", lambda1s=lambda1s, h=h_hp$h,
#'           hp=h_hp$hp, diag=dm, return_raw=TRUE)
#'
#' elts_C <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'             centered=TRUE, diag=dm)
#' ## Centered estimates, elts provided; equivalent to est1 and est2
#' est3 <- estimate(x, "trun_gaussian", elts=elts_C,
#'           symmetric="symmetric", lambda1s=lambda1s, diag=NULL,
#'           return_raw=TRUE)
#'
#' ## Noncentered estimates with Inf penalty on eta; equivalent to est1~3
#' est4 <- estimate(x, "trun_gaussian", elts=NULL, centered=FALSE,
#'           lambda_ratio=0, symmetric="symmetric", lambda1s=lambda1s,
#'           h=h_hp$h, hp=h_hp$hp, diag=dm, return_raw=TRUE)
#' compare_two_results(est1, est2) ## Should be almost all 0
#' compare_two_results(est1, est3) ## Should be almost all 0
#' sum(abs(est4$etas)) ## Should be 0 since non-centered with lambda ratio 0 is equivalent to centered
#' est4$etas <- NULL ## But different from est1 in that the zero etas are returned in est4
#' compare_two_results(est1, est4) ## Should be almost all 0
#'
#'
#' ## Profiled estimates, no elts or h provided, mode and params provided
#' est5 <- estimate(x, "trun_gaussian", elts=NULL, centered=FALSE,
#'           lambda_ratio=Inf, symmetric="or", lambda1s=lambda1s,
#'           mode="min_pow", param1=1, param2=3, diag=dm, return_raw=TRUE)
#'
#' ## Profiled estimates, no elts provided, h provided; equivalent to est5
#' est6 <- estimate(x, "trun_gaussian", elts=NULL, centered=FALSE,
#'           lambda_ratio=Inf, symmetric="or", lambda1s=lambda1s,
#'           h=h_hp$h, hp=h_hp$hp, diag=dm, return_raw=TRUE)
#'
#' elts_NC_P <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=TRUE, diag=dm)
#' ## Profiled estimates, elts provided; equivalent to est5~6
#' est7 <- estimate(x, "trun_gaussian", elts=elts_NC_P, centered=FALSE,
#'           lambda_ratio=Inf, symmetric="or", lambda1s=lambda1s,
#'           diagonal_multiplier=NULL, return_raw=TRUE)
#' compare_two_results(est5, est6) ## Should be almost all 0
#' compare_two_results(est5, est7) ## Should be almost all 0
#'
#'
#' ## Non-centered estimates, no elts or h provided, mode and params provided
#' est8 <- estimate(x, "trun_gaussian", elts=NULL, centered=FALSE,
#'           lambda_ratio=2, symmetric="and", lambda_length=100,
#'           mode="min_pow", param1=1, param2=3, diag=dm, return_raw=TRUE)
#'
#' ## Non-centered estimates, no elts provided, h provided; equivalent to est5
#' est9 <- estimate(x, "trun_gaussian", elts=NULL, centered=FALSE,
#'           lambda_ratio=2, symmetric="and", lambda_length=100,
#'           h=h_hp$h, hp=h_hp$hp, diag=dm, return_raw=TRUE)
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian", centered=FALSE,
#'                 profiled=FALSE, diag=dm)
#' ## Non-centered estimates, elts provided; equivalent to est8~9
#' est10 <- estimate(x, "trun_gaussian", elts=elts_NC_NP,
#'            centered=FALSE, lambda_ratio=2, symmetric="and",
#'            lambda_length=100, diag=NULL, return_raw=TRUE)
#' compare_two_results(est8, est9) ## Should be almost all 0
#' compare_two_results(est8, est10) ## Should be almost all 0
#'
#' @export
estimate <- function(x, setting, elts=NULL, centered=TRUE, symmetric="symmetric", scale="norm", lambda1s=NULL, lambda_length=NULL, lambda_ratio=Inf, mode=NULL, param1=NULL, param2=NULL, h=NULL, hp=NULL, verbose=TRUE, verbosetext="", tol=1e-6, maxit=1000, BIC_refit=TRUE, warmstart=TRUE, diagonal_multiplier=NULL, eBIC_gammas=c(0,0.5,1), return_raw=FALSE){
  ## BIC_refit: calculate BIC (with refit) or not
  ## return_raw: return the raw estimates or not
  ## If elts is given, centered, scale, mode, param1, param2, h, hp are all ignored
  ## If both h and hp given, mode, param1, param2 are ignored; h must be continuous and positive almost thinningwhere and h(0)=0; hp must be the almost thinningwhere derivative of h
  ## Must provide at least one of: 1. elts, 2. mode & param1, 3. h & hp
  n <- nrow(x); p <- ncol(x)
  if (is.null(diagonal_multiplier)) {
    diagonal_multiplier <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
  }
  if (!is.null(elts) && setting != elts$setting){
    stop(paste("The setting you chose was ", setting, ", but elts$setting was ", elts$setting, ".", sep=""))
  }
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
    if (!is.null(elts) && !is.null(elts$diagonal_multiplier) && abs(diagonal_multiplier-elts$diagonal_multiplier) > tol){
      warning("diagonal_multiplier and elts$diagonal_multiplier do not agree. Using elts$diagonal_multiplier instead.")
      diagonal_multiplier <- elts$diagonal_multiplier
    } else {
      if (diagonal_multiplier < 1+tol && p > n){ ## ill-behaved
        warning("p > n and diagonal_multiplier should be larger than 1.")
      } else if (diagonal_multiplier > 1-tol){ ## Allows numerical error; if diagonal_multiplier slightly below 1, set it to 1.
        diagonal_multiplier <- max(1, diagonal_multiplier)
      } else
        stop("diagonal_multiplier must be at least 1.")
    }
  }
  if (!symmetric %in% c("symmetric", "and", "or")){
    stop("Parameter symmetric must be one of \"symmetric\", \"and\", or \"or\".")
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
  if (setting != "gaussian" && is.null(elts) && (is.null(mode) || is.null(param1)) && (is.null(h) || is.null(hp)))
    stop("At least one of 1. elts, 2. mode & param1, or 3. h & hp has to be provided.")
  if (setting != "gaussian" && is.null(elts) && !is.null(mode)){
    if (!is.null(h)){
      if (is.null(hp))
        stop("hp must also be provided if h is given.")
      else{
        warning("mode and h/hp should not be provided at the same time. Using h/hp instead.\n")
        if (abs(h(0)) > tol)
          stop("h(0) must be equal to 0.")
      }
    } else{
      if (is.null(param1))
        stop("param1 (and param2 optionally) must be provided with mode.")
      h_hp <- get_h_hp(mode=mode, para=param1, para2=param2); h <- h_hp$h; hp <- h_hp$hp; remove(h_hp)
      if (is.null(h))
        stop("Mode not supported.")
      if (length(h(1))==0)
        stop("Error occurred in generating h and hp. Possibly due to invalid param1 and/or param2.")
      #if (abs(h(0)) > tol)
      #  stop(paste("h(0)=", h(0), ", larger than 0. Stopped.", sep=""))
    }
  }
  output("Calculating elements necessary for estimation.", verbose, verbosetext)
  if (is.null(elts))
    elts <- get_elts(h, hp, x, setting, centered=centered, profiled_if_noncenter = is.infinite(lambda_ratio), scale=scale, diagonal_multiplier=diagonal_multiplier)
  if (is.null(lambda1s)){
    output("Calculating lower bound for lambda.", verbose, verbosetext)
    lambda_lo <- test_lambda_bounds2(elts, symmetric, lambda_ratio, lower=TRUE, verbose=verbose, tol=tol, maxit=maxit)
    output("Calculating upper bound for lambda.", verbose, verbosetext)
    lambda_hi <- lambda_max(elts = elts, symmetric=symmetric, lambda_ratio = lambda_ratio)
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
  if (!elts$centered) {etas <- matrix(0, nrow=lambda_length, ncol=p)
  } else {etas <- NULL}
  convergeds <- numeric(lambda_length)
  iters <- numeric(lambda_length)
  output("Calculating estimates.", verbose, verbosetext)
  if (verbose){
    #  pb <- txtProgressBar(style=3)
    checkpoints <- ceiling(c(0.1, 0.2, 0.5, 1:10) * lambda_length / 10)
  }
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
    if (verbose && lambda_length > 10 && lambda_index %in% checkpoints)
      output(paste(floor(lambda_index/lambda_length*100), "% done", sep=""), verbose, verbosetext)
    #setTxtProgressBar(pb, lambda_index / lambda_length)
    if (length(res$edges) == elts$p*(elts$p-1) && lambda_index < lambda_length){ ## If all edges are selected for some lambda, end early
      BICs[(lambda_index+1):lambda_length, ] <- Inf
      convergeds[(lambda_index+1):lambda_length] <- 0#convergeds[lambda_index]
      iters[(lambda_index+1):lambda_length] <- 0
      for (li in (lambda_index+1):lambda_length){
        edgess[[li]] <- edgess[[li-1]]
        if (!is.null(etas))
          etas[li,] <- NA#etas[li-1,]
        if (return_raw)
          raw_estimates[[li]] <- matrix(NA,p,p)
      }
      break
    }
  }
  #if (verbose)
  #  close(pb)
  output("Done.", verbose, verbosetext)
  return (list("edgess"=edgess,
               "etas"=switch(changed_from_nc_to_c+1, etas, matrix(0, nrow=lambda_length, ncol=p)),
               "BICs"=BICs[,1:length(eBIC_gammas)],
               "BIC_refits"=switch(BIC_refit+1, NULL, BICs[,(length(eBIC_gammas)+1):(2*length(eBIC_gammas))]),
               "lambda1s"=lambda1s,
               "lambda2s"=switch(1+elts$centered+2*(!elts$centered && elts$profiled), lambda1s/lambda_ratio, NULL, rep(0, length(lambda1s))),
               "converged"=convergeds, "iters"=iters, "raw_estimates"=raw_estimates, "symmetric"=symmetric))
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' res_nc_np <- get_results(elts_NC_NP, symmetric="symmetric",
#'                lambda1=0.35, lambda2=2, previous_res=NULL, is_refit=FALSE)
#' refit(res_nc_np, elts_NC_NP)
#'
#' elts_NC_P <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                centered=FALSE, profiled=TRUE, diag=dm)
#' res_nc_p <- get_results(elts_NC_P, symmetric="symmetric",
#'               lambda1=0.35, lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#' refit(res_nc_p, elts_NC_P)
#'
#' elts_C <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'             centered=TRUE, diag=dm)
#' res_c <- get_results(elts_C, symmetric="or", lambda1=0.35,
#'            lambda2=NULL, previous_res=NULL, is_refit=FALSE)
#' refit(res_c, elts_C)
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
    return (test$crit)
  }
  else {
    res$eta[-res$eta_support] <- 0
    test <- get_results(elts, res$symmetric, lambda1=0, lambda2=0, tol=res$tol, maxit=res$maxit, previous_res=res, is_refit=TRUE)
    return (test$crit)
  }
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' res_nc_np <- get_results(elts_NC_NP, symmetric="symmetric",
#'                lambda1=0.35, lambda2=2, previous_res=NULL,
#'                is_refit=FALSE)
#' eBIC(res_nc_np, elts_NC_NP, BIC_refit=TRUE, gammas=c(0,0.5,1))
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
#' if (!requireNamespace("tmvtnorm", quietly = TRUE)){
#'   stop("Please install package \"tmvtnorm\" first.", call. = FALSE)
#' }
#' require(tmvtnorm)
#' n <- 50
#' p <- 30
#' h_hp <- get_h_hp("min_pow", 1, 3)
#' mu <- rep(0, p)
#' K <- diag(p)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#'
#' elts_NC_NP <- get_elts(h_hp$h, h_hp$hp, x, setting="trun_gaussian",
#'                 centered=FALSE, profiled=FALSE, diag=dm)
#' res_nc_np <- get_results(elts_NC_NP, symmetric="symmetric", lambda1=0.35,
#'                lambda2=2, previous_res=NULL, is_refit=FALSE)
#' get_crit_nopenalty(elts_NC_NP, previous_res=res_nc_np)
#' @export
get_crit_nopenalty <- function(elts, exclude=NULL, exclude_eta=NULL, previous_res=NULL){
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
    if (elts$setting != "gaussian"){
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
#' The Cramér-Rao lower bound (times \code{n}) for estimating the mean parameter from a univariate truncated normal sample with known variance parameter.
#'
#' The Cramér-Rao lower bound (times \code{n}) on the variance for estimating the mean parameter \code{mu} from a univariate truncated normal sample, assuming the true variance parameter \code{sigmasq} is known.
#'
#' @param mu The mean parameter.
#' @param sigmasq The variance parameter.
#' @return A number, the Cramér-Rao lower bound.
#' @details The Cramér-Rao lower bound in this case is defined as \eqn{\sigma^4/var(X-\mu)}.
#' @export
crbound_mu <- function(mu,sigmasq){
  sigma <- sqrt(sigmasq)
  #return (-1/(exp(-mu^2/sigma^2)*(2*(exp(mu^2/2/sigma^2)*mu*sqrt(2*pi)+sigma-2*exp(mu^2/sigma^2)*pi*sigma)+(-exp(mu^2/2/sigma^2)*mu*sqrt(2*pi)+4*exp(mu^2/sigma^2)*pi*sigma)*(2-2*stats::pnorm(mu/sigma))-exp(mu^2/sigma^2)*pi*sigma*(2-2*stats::pnorm(mu/sigma))^2)/(pi*sigma^3*(2*stats::pnorm(mu/sigma,0,1))^2)))
  EZc2 <- -exp(-mu^2/2/sigma^2)*mu*sigma/sqrt(2*pi)+sigma^2/2+sigma^2*(1-stats::pnorm(0,mu,sigma)*2)/2
  E1 <- 1-stats::pnorm(0,mu,sigma)
  EZc <- exp(-mu^2/2/sigma^2)*sigma/sqrt(2*pi)
  return (sigma^4 / (EZc2/E1 - (EZc/E1)^2))
}

#' The Cramér-Rao lower bound (times \code{n}) for estimating the variance parameter from a univariate truncated normal sample with known mean parameter.
#'
#' The Cramér-Rao lower bound (times \code{n}) on the variance for estimating the variance parameter \code{sigmasq} from a univariate truncated normal sample, assuming the true mean parameter \code{mu} is known.
#'
#' @name crbound_sigma
#' @param mu The mean parameter.
#' @param sigmasq The variance parameter.
#' @return A number, the Cramér-Rao lower bound .
#' @details The Cramér-Rao lower bound in this case is defined as \eqn{4\sigma^8/var((X-\mu)^2)}.
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
  h_hp <- get_h_hp(mode,param1,param2)
  hx <- h_hp$h(x)
  hpx <- h_hp$hp(x)
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
  sigma <- sqrt(sigmasq)
  inte <- function(f,from=0,to=Inf){stats::integrate(f,lower=from,upper=to,rel.tol=1e-10)$value}
  adaptinte <- function(f,from=0,to=Inf){cubature::adaptIntegrate(function(t){x<-t/(1-t);1/(1-t)^2*f(x)},lowerLimit=from/(from+1),upperLimit=ifelse(is.infinite(to),1,to/(to+1)),tol=1e-10)$integral}
  all_inte <- function(f,from=0,to=Inf){tryCatch(inte(f,from=from,to=to), error=function(e){adaptinte(f,from,to)})}
  h_hp <- get_h_hp(mode, param1, param2); h <- h_hp$h; hp <- h_hp$hp
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



