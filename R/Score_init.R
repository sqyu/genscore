require(Matrix) ## Required by bdiag
require(igraph)
require(zoo) ## rollmean


#' Random generator of matrices with given eigenvalues.
#'
#' Random generator of matrices with given eigenvalues.
#'
#' @param p A positive integer >= 2, the dimension.
#' @param ev A vector of length \code{p}, the eigenvalues of the output matrix.
#' @param seed A number, the seed for the generator. Ignored if \code{NULL}.
#' @details The function randomly fills a \code{p} by \code{p} matrix with independent \eqn{Normal(0,1)} entries, takes the \code{Q} matrix from its \code{QR} decomposition, and returns \eqn{Q' diag(ev) Q}.
#' @return A \code{p} by \code{p} matrix whose eigenvalues equal to \code{ev}.
#' @examples
#' p <- 30
#' eigen_values <- (0.1*p-1+1:p)/p
#' K <- ran_mat(p, ev=eigen_values, seed=1)
#' sort(eigen(K)$val)-eigen_values
#' @importFrom Rdpack reprompt
#' @export
ran_mat <- function (p, ev = stats::runif(p, 0, 10), seed=NULL) {
  if (length(ev) != p) stop("length of ev must equal to p.")
  p <- as.integer(p)
  if (p < 2) stop("p must be a positive integer >= 2.")
  if (!is.null(seed)) set.seed(seed)
  Z <- matrix(ncol=p, stats::rnorm(p^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  #R <- qr.R(decomp); d <- diag(R); ph <- d / abs(d); O <- Q %*% diag(ph); Z <- t(O) %*% diag(ev) %*% O
  Z <- t(Q) %*% diag(ev) %*% Q
  return(Z)
}

#' Random generator of inverse covariance matrices.
#'
#' Random generator of inverse covariance matrices.
#'
#' @param mode A string, see details.
#' @param p A positive integer >= 2, the dimension.
#' @param seed A number, the seed for the generator. Ignored if \code{NULL} or \code{mode == "band"} or \code{mode == "chain"}.
#' @param spars A number, see details. Ignored if \code{mode == "chain"}. Default to 1.
#' @param eig A positive number, the minimum eigenvalue of the returned matrix. Default to 0.1.
#' @param subgraphs A positive integer, the number of subgraphs for the \code{"sub"} mode. Note that \code{p} must be divisible by \code{subgraphs}.
#' @return A \code{p} by \code{p} inverse covariance matrix. See details.
#' @details
#' \describe{
#' The function generates an inverse covariance matrix according to the \code{mode} argument as follows. The diagonal entries of the matrix are set to the same value such that the minimum eigenvalue of the returned matrix is equal to \code{eig}.
#'     \item{"random"}{Takes the \code{Q} matrix from the \code{QR} decomposition of a \code{p} by \code{p} random matrix with independent \eqn{Normal(0,1)} entries, and calculates \eqn{Q' diag(ev) Q}. Randomly zeros out its upper triangular entries using independent uniform Bernoulli(\code{spars}) variables, and then symmetrizes the matrix using the upper triangular part.}
#'     \item{"sub"}{Constructs a block diagonal matrix with \code{subgraphs} disconnected subgraphs with equal number of nodes. In each subgraph, takes each entry independently from \eqn{Uniform(0.5,1)}, and randomly zeros out its upper triangular entries using independent uniform Bernoulli(\code{spars}) variables, and finally symmetrizes the matrix using the upper triangular part. The construction from Section 4.2 of \insertCite{lin16;textual}{genscore}.}
#'     \item{"er"}{Constructs an Erdős–Rényi game with probability \code{spars}, and sets the edges to independent \eqn{Uniform(0.5,1)} variables, and finally symmetrizes the matrix using the lower triangular entries.}
#'     \item{"band"}{Constructs a banded matrix so that the \code{(i,j)}-th matrix is nonzero if and only if \eqn{|i-j|<=spars}, and is equal to \eqn{1-|i-j|/(spars+1)} if \eqn{i!=j}.}
#'     \item{"chain"}{A chain graph, where the \code{(i,j)}-th matrix is nonzero if and only if \eqn{|i-j|<=1}, and is equal to 0.5 if \eqn{|i-j|==1}. A special case of the \code{"band"} construction with \code{spars} equal to 1.}
#' }
#' @examples
#' p <- 100
#' K1 <- cov_cons("random", p, seed = 1, spars = 0.05, eig = 0.1)
#' K2 <- cov_cons("sub", p, seed = 2, spars = 0.5, eig = 0.1, subgraphs=10)
#' K3 <- cov_cons("er", p, seed = 3, spars = 0.05, eig = 0.1)
#' K4 <- cov_cons("band", p, spars = 2, eig = 0.1)
#' K5 <- cov_cons("chain", p, eig = 0.1)
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{lin16}{genscore}
#' @export
cov_cons <- function(mode, p, seed=NULL, spars=1, eig=0.1, subgraphs=1){
  p <- as.integer(p); subgraphs <- as.integer(subgraphs)
  if (eig <= 0) stop("eig must be positive.")
  if (p < 2) stop("p must be a positive integer >= 2.")
  if (!is.null(seed)) set.seed(seed)
  if (mode == "random"){
    ## spars: p in binomial for the whole graph
    K <- ran_mat(p)
    K[upper.tri(K)] <- K[upper.tri(K)] * matrix(stats::rbinom(p^2,1,spars),p,p)[upper.tri(diag(p))]
    K[lower.tri(K)] <- t(K)[lower.tri(K)]
    K <- K + diag(p) * (eig - min(eigen(K)$values))
  } else if (mode == "sub"){ # From Section 4.2 of Lin et al (2016)
    ## spars: p in binomial for each subgraph
    if (!requireNamespace("Matrix", quietly = TRUE)){
      stop("Please install package \"Matrix\" first.", call. = FALSE)
    }
    if (subgraphs < 1 || p %% subgraphs) {stop("subgraphs must be a positive integer and p must be an exact multiple of subgraphs.")}
    p_sub <- p / subgraphs
    K <- as.matrix(Matrix::bdiag(lapply(1:subgraphs, function(x){mat <- matrix(stats::runif(p_sub^2, 0.5, 1) * stats::rbinom(p_sub^2, 1, spars), p_sub, p_sub);
      mat[upper.tri(mat, diag = TRUE)] <- 0; mat <- mat + t(mat);
      mat <- mat + diag(p_sub) * (eig - min(eigen(mat)$values))})))
  } else if (mode == "er"){
    ## spars: p in binomial for the whole graph
    if (!requireNamespace("igraph", quietly = TRUE)){
      stop("Please install package \"igraph\" first.", call. = FALSE)
    }
    K <- as.matrix(igraph::get.adjacency(igraph::erdos.renyi.game(p, spars))) # Not sure why t(K) below would cause an error otherwise
    K <- K * matrix(stats::runif(p^2, 0.5, 1), p, p); K[upper.tri(K)] <- 0; K <- K + t(K)
    K <- K + diag(p) * (eig - min(eigen(K)$values))
  } else if (mode == "band"){
    ## spars: bandwidth, i.e. Kij != 0 iff |i-j| <= spars
    ## non-random so seed not needed
    K <- matrix(0,p,p)
    entries <- seq(1,0,length.out=spars+2)[2:(spars+1)]
    for (i in 1:spars)
      K[(1:(p-i)-1)*p+(1+i):p] <- K[((1+i):p-1)*p+1:(p-i)] <- entries[i]
    K <- K + diag(p) * (eig - min(eigen(K)$values))
  } else if (mode == "chain"){
    ## equivalent to a 1-band graph
    return (cov_cons(mode="band", p=p, seed=seed, spars=1, eig=eig))
  } else{
    stop("Mode not supported.")
  }
  return (K)
}

#' Calculates the true and false positive rates given the estimated and true edges.
#'
#' Calculates the true and false positive rates given the estimated and true edges.
#'
#' @param edges A vector of indices corresponding to the estimated edges. Should not contain the diagonals.
#' @param true_edges A vector of indices corresponding to the true edges.
#' @param p A positive integer, the dimension.
#' @return A vector containing the true positive rate and the false positive rate.
#' @examples
#' n <- 40
#' p <- 50
#' mu <- rep(0, p)
#' tol <- 1e-8
#' K <- cov_cons(mode="sub", p=p, seed=1, spars=0.2, eig=0.1, subgraphs=10)
#' true_edges <- which(abs(K) > tol & diag(p) == 0)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' set.seed(1)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' est <- estimate(x, "trun_gaussian", elts=NULL, centered=TRUE,
#'          symmetric="symmetric", lambda_length=100, mode="min_pow",
#'          param1=1, param2=3, diagonal_multiplier=dm,)
#' # Apply tp_fp to each estimated edges set for each lambda
#' TP_FP <- t(sapply(est$edgess, function(edges){tp_fp(edges, true_edges, p)}))
#' par(mfrow=c(1,1), mar=c(5,5,5,5))
#' plot(c(), c(),  ylim=c(0,1), xlim=c(0,1), cex.lab=1, main = "ROC curve",
#'   xlab="False Positives", ylab="True Positives")
#' points(TP_FP[,2], TP_FP[,1], type="l")
#' points(c(0,1), c(0,1), type = "l", lty = 2)
#'
#' @export
tp_fp <- function(edges, true_edges, p){
  edges <- sort(edges); true_edges <- sort(true_edges)
  check_input <- function(e, name){
    if (length(e)){
      if (e[length(e)] > p^2) {stop(name, " should not contain indices larger than p^2.")}
      if (e[1] < 1) {stop(name, " should not contain numbers smaller than 1.")}
      if (e[length(e)] == p^2 || any(e%%p == ceiling(e/p))) {stop("The diagonals should not be included in ", name, ".")}
    }
  }
  check_input(edges, "edges"); check_input(true_edges, "true_edges")
  return (c(length(intersect(true_edges, edges)) / length(true_edges),
            length(setdiff(edges, true_edges)) / (p^2-p-length(true_edges))))
}

#' Calculates the AUC of an ROC curve.
#'
#' Calculates the area under an ROC curve (AUC).
#'
#' @param tpfp A matrix with two columns, the true positive and the false positive rates.
#' @return A number between 0 and 1, the area under the curve (AUC).
#' @examples
#' n <- 40
#' p <- 50
#' mu <- rep(0, p)
#' tol <- 1e-8
#' K <- cov_cons(mode="sub", p=p, seed=1, spars=0.2, eig=0.1, subgraphs=10)
#' true_edges <- which(abs(K) > tol & diag(p) == 0)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' set.seed(1)
#' x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'        lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'        burn.in.samples = 100, thinning = 10)
#' est <- estimate(x, "trun_gaussian", elts=NULL, centered=TRUE,
#'          symmetric="symmetric", lambda_length=100, mode="min_pow",
#'          param1=1, param2=3, diagonal_multiplier=dm)
#' # Apply tp_fp to each estimated edges set for each lambda
#' TP_FP <- t(sapply(est$edgess, function(edges){tp_fp(edges, true_edges, p)}))
#' par(mfrow=c(1,1), mar=c(5,5,5,5))
#' auc <- AUC(TP_FP)
#' plot(c(), c(),  ylim=c(0,1), xlim=c(0,1), cex.lab=1,
#'   main=paste("ROC curve, AUC",round(auc,4)), xlab="False Positives",
#'   ylab="True Positives")
#' points(TP_FP[,2], TP_FP[,1], type="l")
#' points(c(0,1), c(0,1), type = "l", lty = 2)
#' @export
AUC <- function(tpfp){
  if (!requireNamespace("zoo", quietly = TRUE)){
    stop("Please install package \"zoo\" first.", call. = FALSE)
  }
  if (min(tpfp) < 0 || max(tpfp) > 1) {stop("All values in tpfp must be between 0 and 1.")}
  if (ncol(tpfp) != 2) {stop("tpfp must be a matrix of two columns, namely the true and false positive rates.")}
  tpfp <- tpfp[order(tpfp[,2]), ]
  sum(diff(tpfp[,2]) * zoo::rollmean(tpfp[,1],2))
}

#' Finds the max index in a vector that does not exceed a target number.
#'
#' Finds the max index in a vector that does not exceed a target number.
#'
#' @param vals A vector of numbers.
#' @param target A number. Must not be smaller than \code{vals[start]}.
#' @param start A number, the starting index; default to 1. Must be such that \code{vals[start] <= target}.
#' @return The max index \code{i} such that \code{vals[i] <= target} and \code{i >= start}.
#' @examples
#' for (i in 1:100) {
#'    vals <- 1:i
#'    for (start in 1:i)
#'       for (target in seq(start, i+0.5, by=0.5))
#'          if (find_max_ind(vals, target, start) != floor(target))
#'             stop()
#' }
#' @export
find_max_ind <- function(vals, target, start=1){
  if (start < 1 || start > length(vals)) stop("start must be between 1 and length(vals).")
  if (vals[start] > target) stop("vals[start] must be <= target.")
  if (vals[length(vals)] <= target) return (length(vals))
  l <- start; r <- length(vals)
  while (r-l>1){
    m <- floor((r+l)/2)
    if (vals[m] <= target) l <- m  # l always satisfies vals[l] <= target
    else r <- m
  }
  return (ifelse(vals[r]<=target,r,l))
}


#' Takes the vertical average of ROC curves.
#'
#' Takes the vertical average of ROC curves using algorithm 3 from \insertCite{faw06;textual}{genscore}. The resulting ROC curve preserves the average AUC.
#'
#' @param rocs A list of ROC curves, each of which is a matrix with two columns corresponding to the true positive and false positive rates, respectively.
#' @param num_true_edges A positive integer, the number of true edges
#' @param p A positive integer, the dimension
#' @return The averaged ROC curve, a matrix with 2 columns and \code{(p^2-p-num_true_edges+1)} rows.
#' @references
#' \insertRef{faw06}{genscore}
#' @examples
#' n <- 40
#' p <- 50
#' mu <- rep(0, p)
#' tol <- 1e-8
#' K <- cov_cons(mode="sub", p=p, seed=1, spars=0.2, eig=0.1, subgraphs=10)
#' true_edges <- which(abs(K) > tol & diag(p) == 0)
#' dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
#' ROCs <- list()
#' par(mfrow=c(2,2), mar=c(5,5,5,5))
#' for (i in 1:3){
#'   set.seed(i)
#'   x <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = solve(K),
#'          lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs",
#'          burn.in.samples = 100, thinning = 10)
#'   est <- estimate(x, "trun_gaussian", elts=NULL, centered=TRUE,
#'            symmetric="symmetric", lambda_length=100, mode="min_pow",
#'            param1=1, param2=3, diag=dm)
#'   # Apply tp_fp to each estimated edges set for each lambda
#'   TP_FP <- t(sapply(est$edgess, function(edges){tp_fp(edges, true_edges, p)}))
#'   ROCs[[i]] <- TP_FP
#'   plot(c(), c(),  ylim=c(0,1), xlim=c(0,1), cex.lab=1,
#'     main=paste("ROC, trial ",i,", AUC ",round(AUC(TP_FP),4),sep=""),
#'     xlab="False Positives", ylab="True Positives")
#'   points(TP_FP[,2], TP_FP[,1], type="l")
#'   points(c(0,1), c(0,1), type = "l", lty = 2)
#' }
#' average_ROC <- avgrocs(ROCs, length(true_edges), p)
#' plot(c(), c(),  ylim=c(0,1), xlim=c(0,1), cex.lab=1,
#'   main=paste("Average ROC, AUC",round(AUC(average_ROC),4)),
#'   xlab="False Positives", ylab="True Positives")
#' points(average_ROC[,2], average_ROC[,1], type="l")
#' points(c(0,1), c(0,1), type = "l", lty = 2)
#' @export
avgrocs <- function(rocs, num_true_edges, p){
  ### rocs: list of roc curves; each curve has 2 columns
  if (length(num_true_edges) != 1) {stop("Attention: the second argument has been changed from true_edges to num_true_edges!")}
  ngrid <- (p^2 - p - num_true_edges)
  if (ngrid < 0) {stop("num_true_edges must be at most p^2-p.")}
  if (any(sapply(rocs, ncol) != 2)) {stop("Each element in rocs must be a matrix of two columns, namely the true and false positive rates.")}
  if (any(sapply(rocs, min) < 0) || any(sapply(rocs, max) > 1)) {stop("All values in rocs must be between 0 and 1.")}
  tpravg <- numeric(ngrid+1)
  pointers <- rep(1, length(rocs))
  for (i in 1:length(rocs)) ## Sort each ROC by increasing FPR
    rocs[[i]] <- rocs[[i]][order(rocs[[i]][,2]), ]
  interpolate <- function(rocp1, rocp2, X){
    slope <- (rocp2[1] - rocp1[1])/(rocp2[2] - rocp1[2])
    return (rocp1[1] + slope*(X-rocp1[2]))
  }
  tpr_for_fpr <- function(fpr, roc, curve_index){
    #pointers[curve_index] <<- find_max_ind(roc[,2], fpr+(1e-6)/ngrid, pointers[curve_index])
    #Search the largest index in roc[,2] that is <= fpr+tol. Binary search in the previous line might be even slower than linear search below.
    while (pointers[curve_index] < nrow(roc) && roc[pointers[curve_index]+1, 2] <= fpr+(1e-6)/ngrid)
      pointers[curve_index] <<- pointers[curve_index] + 1
    if (abs(roc[pointers[curve_index], 2]-fpr)*ngrid < 1e-6) # If exact match of fpr
      return (roc[pointers[curve_index], 1])
    else{ # Otherwise interpolate
      if (pointers[curve_index] < nrow(roc)) ## Needed since roc[pointers[curve_index]+1, ] might be out of bound
        return (interpolate(roc[pointers[curve_index], ], roc[pointers[curve_index]+1, ], fpr))
      else
        return (interpolate(roc[pointers[curve_index], ], c(1,1), fpr))
    }
  }
  for (i in 1:(ngrid+1))
    tpravg[i] <- mean(sapply(1:length(rocs), function(ci){tpr_for_fpr((i-1)/ngrid, rocs[[ci]], ci)}))
  return (cbind(tpravg, (0:ngrid)/ngrid))
}

#' Computes the sum of absolute differences in the finite non-NA/NULL elements between two vectors.
#'
#' Computes the sum of absolute differences in the finite non-\code{NA}/\code{NULL} elements between two vectors.
#'
#' @param l1 A vector.
#' @param l2 A vector.
#' @param relative A boolean, default to \code{FALSE}. If \code{TRUE}, returns the relative difference (sum of absolute differences divided by the elementwise minimum between \code{l1} and \code{l2}).
#' @return The sum of (relative) absolute differences in \code{l1} and \code{l2}, or a positive integer if two vectors differ in length or hold \code{NA}, \code{NULL} or \code{Inf} values in different places.
diff_vecs <- function(l1, l2, relative=FALSE){
  if (length(l1) != length(l2)) {return (abs(length(l1)-length(l2)))}
  tmp1 <- which(is.null(l1) | is.na(l1) | is.infinite(l1))
  tmp2 <- which(is.null(l2) | is.na(l2) | is.infinite(l2))
  if (length(tmp1) != length(tmp2) || (length(tmp1) && length(tmp2) && any(tmp1 != tmp2))) {return(abs(length(tmp1)-length(tmp2))+abs(sum(tmp1)-sum(tmp2)))}
  if (relative)
    return (sum(abs((l1[-tmp1]-l2[-tmp2])/pmin(abs(l1[-tmp1]), abs(l2[-tmp2])))))
  else
    return (sum(abs(l1[-tmp1]-l2[-tmp2])))
}

#' Computes the sum of absolute differences between two lists.
#'
#' Computes the sum of absolute differences between two lists using diff_vecs().
#'
#' @param l1 A list.
#' @param l2 A list.
#' @param name A string, default to \code{NULL}. If not \code{NULL}, computes the differences in the \code{l1[[name]]} and \code{l2[[name]]}.
#' @return
#' Returns the sum of absolute differences between \code{l1} and \code{l2} if \code{name} is \code{NULL}, or that between \code{l1[[name]]} and \code{l2[[name]]} otherwise. If \code{name} is not \code{NULL} and if \code{name} is in exactly one of \code{l1} and \code{l2}, returns \code{Inf}; if \code{name} is in neither, returns \code{NA}. Exception: Returns a positive integer if the two elements compared hold \code{NA}, \code{NULL} or \code{Inf} values in different places.
diff_lists <- function(l1, l2, name=NULL){
  if (length(l1) != length(l2)) {return (abs(length(l1)-length(l2)))}
  if (any(sapply(l1, length) != sapply(l2, length))) {return (sum(abs(sapply(l1, length) - sapply(l2, length))))}
  if (is.null(name))
    return (sum(abs(sapply(1:length(l1), function(i){diff_vecs(l1[[i]], l2[[i]])}))))
  else{
    if ((name %in% l1) && (name %in% l2))
      return (sum(abs(sapply(1:length(l1), function(i){diff_vecs(l1[[i]][[name]], l2[[i]][[name]])}))))
    if (xor(name %in% l1, name %in% l2))
      return (Inf)
    return (NA)
  }
}

#' Compares two lists returned from get_results().
#'
#' Compares two lists returned from \code{get_results()}.
#'
#' @param res A res list returned from \code{get_results()}.
#' @param res2 A res list returned from \code{get_results()}.
#' @return A list of numbers all of which should be close to 0 if \code{res} and \code{res2} are expected to be the same.
#' @export
compare_two_sub_results <- function(res, res2){
  d <- c("p"=res$p-res2$p,
         "K"=sum(abs(res$K-res2$K)), "lambda1"=sum(abs(res$lambda1-res2$lambda1)),
         "tol"=abs(res$tol-res2$tol), "maxit"=abs(res$maxit-res2$maxit),
         "converged"=abs(res$converged-res2$converged),
         "crit"=abs(res$crit-res2$crit),
         "diagonals_with_multiplier"=sum(abs(res$diagonals_with_multiplier-res2$diagonals_with_multiplier)),
         "edges"=diff_vecs(res$edges, res2$edges))
  for (name in c("eta"))
    if (!is.null(res[[name]]) || !is.null(res2[[name]])){
      if (is.null(res[[name]]) || is.null(res2[[name]]))
        d[name] <- Inf
      else
        d[name] <- sum(abs(res[[name]]-res2[[name]]))
    }
  if (!is.null(res[["eta_support"]]) || !is.null(res2[["eta_support"]])){
    if (is.null(res[["eta_support"]]) || is.null(res2[["eta_support"]]))
      d["eta_support"] <- Inf
    else
      d["eta_support"] <- diff_vecs(res$eta_support, res2$eta_support)
  }
  return (d)
}

#' Compares two lists returned from estimate().
#'
#' Compares two lists returned from \code{estimate}().
#'
#' @param res A res list returned from \code{estimate()}.
#' @param res2 A res list returned from \code{estimate()}.
#' @return A list of numbers all of which should be close to 0 if \code{res} and \code{res2} are expected to be the same.
#' @export
compare_two_results <- function(res, res2){
  d <- c("edgess"=diff_lists(res$edgess, res2$edgess),
         "BICs"=sum(abs(res$BICs-res2$BICs)),
         "K"=diff_lists(res$raw_estimates, res2$raw_estimates),
         "lambda1s"=diff_vecs(res$lambda1s, res2$lambda1s),
         "lambda2s"=diff_vecs(res$lambda2s, res2$lambda2s)
         )
  if (!is.null(res[["etas"]]) || !is.null(res2[["etas"]])){
    if (is.null(res[["etas"]]) || is.null(res2[["etas"]]))
      d["etas"] <- Inf
    else
      d["etas"] <- diff_vecs(res$etas, res2$etas)
  }
  for (name in c("BICs", "BIC_refits"))
    d[name] <- diff_vecs(res[[name]], res2[[name]], relative=TRUE)
  return (d)
}

