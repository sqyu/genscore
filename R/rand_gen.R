
#' Generates translated and truncated exponential variables.
#' 
#' Generates translated and truncated exponential variables.
#' 
#' @param n An integer, the number of samples to return.
#' @param lo A double, the lower limit of the distribution, cannot be \code{-Inf}.
#' @param hi A double, the upper limit of the distribution.
#' @details Returns \code{n} random variables from the translated truncated exponential distribution with density \eqn{\exp(-(x-lo))/(1-\exp(lo-hi))}{exp(-(x-lo))/(1-exp(lo-hi))} on \code{[lo,hi]}.
#' @return \code{n} random variables from the translated truncated exponential distribution.
#' @examples
#' hist(rexp_truncated(1e4, 0, Inf), breaks=200)
#' hist(rexp_truncated(1e4, 10, 12), breaks=200)
#' hist(rexp_truncated(1e4, -2, 2), breaks=200)
#' hist(rexp_truncated(1e4, -10, 0), breaks=200)
#' hist(rexp_truncated(1e4, -100, Inf), breaks=200)
#' hist(rexp_truncated(1e4, -100, -95), breaks=200)
#' @export
rexp_truncated <- function(n, lo, hi) {
  if (lo > hi) stop("lo must be smaller than hi.")
  if (is.infinite(lo)) stop("lo cannot be infinite.")
  lo - log(stats::runif(n, exp(lo-hi), 1))
}

#' Generates centered laplace variables with scale 1.
#' 
#' Generates centered laplace variables with scale 1.
#' 
#' @param n An integer, the number of samples to return.
#' @param lo A double, the lower limit of the distribution.
#' @param hi A double, the upper limit of the distribution.
#' @details Returns \code{n} random variables from the truncated laplace distribution with density proportional to \eqn{\exp(-|x|)}{exp(-|x|)} on \code{[lo,hi]}.
#' @return \code{n} random variables from the truncated laplace distribution.
#' @examples
#' hist(rlaplace_truncated_centered(1e4, -Inf, Inf), breaks=200)
#' hist(rlaplace_truncated_centered(1e4, 0, Inf), breaks=200)
#' hist(rlaplace_truncated_centered(1e4, 10, 12), breaks=200)
#' hist(rlaplace_truncated_centered(1e4, -2, 2), breaks=200)
#' hist(rlaplace_truncated_centered(1e4, -10, 0), breaks=200)
#' hist(rlaplace_truncated_centered(1e4, -100, Inf), breaks=200)
#' hist(rlaplace_truncated_centered(1e4, -100, -95), breaks=200)
#' @export
rlaplace_truncated_centered <- function(n, lo, hi) {
  if (n < 1) return (c())
  if (lo > hi) stop("lo must be smaller than hi.")
  if (lo >= 0) {rexp_truncated(n, lo, hi)
  } else if (hi <= 0) {-rexp_truncated(n, -hi, -lo)
  } else {
    left_prob <- (1 - exp(lo)) / (2 - exp(lo) - exp(-hi))
    replicate(n, if (stats::runif(1,0,1) < left_prob) {
      -rexp_truncated(1, 0, -lo)
    } else {rexp_truncated(1, 0, hi)})
  }
}

#' Finds the index of the bin a number belongs to using binary search.
#' 
#' Finds the index of the bin a number belongs to using binary search.
#' 
#' @param arr A vector of size at least 2.
#' @param l An integer between 1 and \code{length(arr)}. Must be smaller than \code{1}.
#' @param r An integer between 1 and \code{length(arr)}. Must be larger than \code{l}.
#' @param x A number. Must be within the range of [\code{arr[l]}, \code{arr[r]}].
#' @details Finds the smallest index \code{i} such that \code{arr[i] <= x <= arr[i+1]}.
#' @return The index \code{i} such that \code{arr[i] <= x <= arr[i+1]}.
#' @examples
#' binarySearch_bin(1:10, 1, 10, seq(1, 10, by=0.5))
#' binarySearch_bin(1:10, 5, 8, seq(5, 8, by=0.5))
#' @export
binarySearch_bin <- Vectorize(function(arr, l, r, x){
  if (r > l + 1) {
    mid <- floor((l + r) / 2)
    if (arr[mid] >= x)
      return (binarySearch_bin(arr, l, mid, x))
    return (binarySearch_bin(arr, mid, r, x))
  }
  return (l)
}, vectorize.args=c("x"))

#' Finds the index of the bin a number belongs to using naive search.
#' 
#' Finds the index of the bin a number belongs to using naive search.
#' 
#' @param arr A vector of size at least 2.
#' @param x A number. Must be within the range of [\code{arr[1]}, \code{arr[length(arr)]}].
#' @details Finds the smallest index \code{i} such that \code{arr[i] <= x <= arr[i+1]}.
#' @return The index \code{i} such that \code{arr[i] <= x <= arr[i+1]}.
#' @examples
#' naiveSearch_bin(1:10, seq(1, 10, by=0.5))
#' @export
naiveSearch_bin <- Vectorize(function(arr, x){
  for (i in seq_len(length(arr)-1))
    if (arr[i + 1] >= x)
      return (i)
  stop("Invalid x.")
}, vectorize.args=c("x"))

#' Finds the index of the bin a number belongs to.
#' 
#' Finds the index of the bin a number belongs to.
#' 
#' @param arr A vector of size at least 2.
#' @param x A number. Must be within the range of [\code{arr[1]}, \code{arr[length(arr)]}].
#' @details Finds the smallest index \code{i} such that \code{arr[i] <= x <= arr[i+1]}. Calls \code{binarySearch_bin()} if \code{length(arr) > 8} and calls \code{naiveSearch_bin()} otherwise.
#' @return The index \code{i} such that \code{arr[i] <= x <= arr[i+1]}.
#' @examples
#' search_bin(1:10, seq(1, 10, by=0.5))
#' @export
search_bin <- function(arr, x){
  if (length(arr) < 2 || any(x < arr[1]) || any(x > arr[length(arr)]))
    stop("arr must have at least 2 elements and x must be in its range.")
  if (length(arr) > 8) return (binarySearch_bin(arr, 1, length(arr), x))
  else return (naiveSearch_bin(arr, x))
}

#' Generates laplace variables truncated to a finite union of intervals.
#' 
#' Generates laplace variables truncated to a finite union of intervals.
#' 
#' @param n An integer, the number of samples to return.
#' @param lefts A vector of numbers, must have the same length as \code{rights}. A non-empty vector of numbers (may contain \code{-Inf}), the left endpoints of a domain defined as a union of intervals. It is required that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}.
#' @param rights A vector of numbers, must have the same length as \code{lefts}. A non-empty vector of numbers (may contain \code{Inf}), the right endpoints of a domain defined as a union of intervals. It is required that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}.
#' @param m A number, the location parameter of the laplace distribution.
#' @param s A number, the scale/dispersion parameter of the laplace distribution.
#' @details Returns \code{n} random variables from the truncated laplace distribution with density proportional to \eqn{\exp(-|x-m|/s)}{exp(-|x-m|/s)} truncated to the domain defined by the union of [\code{lefts[i]}, \code{rights[i]}].
#' @return \code{n} random variables from the truncated laplace distribution.
#' @examples
#' hist(rlaplace_truncated(1e4, -Inf, Inf), breaks=200)
#' hist(rlaplace_truncated(1e4, c(0, 5), c(2, Inf), m=2, s=3), breaks=200)
#' hist(rlaplace_truncated(1e4, c(-Inf, 0, 3), c(-3, 1, 12), m=8, s=4), breaks=200)
#' hist(rlaplace_truncated(1e4, c(-5, 0), c(-2, 2), s=0.8), breaks=200)
#' hist(rlaplace_truncated(1e4, c(-10, 1), c(-7, 10), m=-4), breaks=200)
#' hist(rlaplace_truncated(1e4, c(-Inf, 100), c(-100, Inf), m=100), breaks=200)
#' hist(rlaplace_truncated(1e4, c(-Inf, 100), c(-100, Inf), m=-100), breaks=200)
#' hist(rlaplace_truncated(1e4, c(-100, -90), c(-95, -85), s=2), breaks=200)
#' @export
rlaplace_truncated <- function(n, lefts, rights, m=0, s=1) {
  lefts <- (lefts - m) / s; rights <- (rights - m) / s
  if (length(lefts) == 1)
    return (rlaplace_truncated_centered(n, lefts, rights) * s + m)
  norm_consts <- rep(length(lefts))
  for (i in 1:length(lefts)) {
    if (lefts[i] < 0) {
      if (rights[i] > 0) # left < 0 < right: exp(-|x|) on [left, right]
        norm_consts[i] <- 2 - exp(lefts[i]) - exp(-rights[i])
      else # left < right < 0: exp(x) on [left, right]
        norm_consts[i] <- exp(rights[i]) - exp(lefts[i])
    } else # 0 < left < right: exp(-x) on [left, right]
      norm_consts[i] <- exp(-lefts[i]) - exp(-rights[i])
  }
  norm_consts <- norm_consts / sum(norm_consts)
  bins <- data.frame(table(search_bin(c(0, cumsum(norm_consts)), stats::runif(n, 0, 1)))) # bin index -> number of samples in the bin
  bins[,1] <- as.integer(as.character(bins[,1]))
  return (Reduce("c", sapply(1:nrow(bins), 
          function(i){rlaplace_truncated_centered(bins[i,2], lefts[bins[i,1]], rights[bins[i,1]]) * s + m})))
}

#' Generates random numbers from a finite union of intervals.
#' 
#' Generates random numbers from a finite union of intervals.
#' 
#' @param n An integer, the number of samples to return.
#' @param lefts A vector of numbers, must have the same length as \code{rights}. A non-empty vector of numbers (may contain \code{-Inf}), the left endpoints of a domain defined as a union of intervals. It is required that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}.
#' @param rights A vector of numbers, must have the same length as \code{lefts}. A non-empty vector of numbers (may contain \code{Inf}), the right endpoints of a domain defined as a union of intervals. It is required that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}.
#' @details For each sample, a random bin \code{i} is uniformly chosen from 1 through \code{length(lefts)}; if the \code{lefts[i]} and \code{rights[i]} define a finite interval, a random uniform variable is drawn from the interval; if the interval is infinite, a truncated laplace variable with location 0 and scale 1 is drawn. Used for randomly generating initial points for generators of truncated multivariate distributions.
#' @return \code{n} random numbers from the union of intervals.
#' @examples
#' hist(random_init_uniform(1e4, -Inf, Inf), breaks=200)
#' hist(random_init_uniform(1e4, c(0, 5), c(2, Inf)), breaks=200)
#' hist(random_init_uniform(1e4, c(-Inf, 0, 3), c(-3, 1, 12)), breaks=200)
#' hist(random_init_uniform(1e4, c(-5, 0), c(-2, 2)), breaks=200)
#' hist(random_init_uniform(1e4, c(-10, 1), c(-7, 10)), breaks=200)
#' hist(random_init_uniform(1e4, c(-Inf, 100), c(-100, Inf)), breaks=200)
#' hist(random_init_uniform(1e4, c(-100, -90), c(-95, -85)), breaks=200)
#' @export
random_init_uniform <- function(n, lefts, rights){
  num_int <- length(lefts)
  xinit <- replicate(n, {
    bi <- sample(num_int, 1); # Randomly pick one bin
    if ((bi == 1 && is.infinite(lefts[bi])) || 
        bi == num_int && is.infinite(rights[bi])) # If bin originally infinite
      rlaplace_truncated_centered(1, lefts[bi], rights[bi]) # Laplace
    else # If bin originally finite
      stats::runif(1, lefts[bi], rights[bi]) # Uniform
  })
  return (xinit)
}

#' Generates a random point in the (p-1)-simplex.
#' 
#' Generates a random point in the \code{(p-1)}-simplex.
#' 
#' @param p An integer, the dimension.
#' @param rdist A function that generates a random number when called using \code{rdist(1)}. Must return a non-zero number with a large enough probability.
#' @param tol A small positive number. Samples are regenerated until each of the \code{p} components are at least of size \code{tol}. 
#' @param maxtimes An integer, maximum number of attempts.
#' @details \code{p} numbers are generated from \code{rdist} and their absolute values are normalized to sum to 1. This will be repeated up to \code{maxtimes} times until all \code{p} components are larger than or equal to \code{tol}.
#' @return A random point (\code{p}-vector) in the \code{(p-1)}-simplex, i.e. \code{sum(x) == 1 && x > 0}.
#' @examples
#' random_init_simplex(100, stats::rnorm, 1e-10, 100)
#' @export
random_init_simplex <- function(p, rdist=stats::rnorm, tol=1e-10, maxtimes=100) {
  # rdist: A random generator such that rdist(1) returns a non-zero number with a large probability.
  count <- 0
  x <- abs(replicate(p, rdist(1)))
  x <- x / sum(x)
  while (any(x <= tol)) {
    if (count > maxtimes)
      stop("Failed to generate x >= tol after ", maxtimes, " times. Try using another rdist() that generates a non-zero number with a larger probability, or set tol to a smaller number.")
    x <- abs(replicate(p, rdist(1)))
    x <- x / sum(x)
    count <- count + 1
  }
  return (x)
}

#' Randomly generate an initial point in the domain defined by a single polynomial with no negative coefficient.
#' 
#' Randomly generate an initial point in the domain defined by a single polynomial with no negative coefficient.
#' 
#' @param domain A list returned from \code{make_domain()} that represents the domain. Currently only supports \code{domain$type == "polynomial" && length(domain$ineqs) == 1}. If \code{domain$ineqs[[1]]$uniform == FALSE}, \code{domain$ineqs[[1]]$coeffs} must not contain negative numbers.
#' @details If inequality is uniform, find the uniform bound for each component and generate each coordinate using \code{random_init_uniform()}.
#' Otherwise, first randomly generate centered laplace variables for components with coefficient 0 (free variables).
#' Then assign a \code{quota} of \code{eq$const / length(nonzero_coefficient)} to each coordinate (so that each \cr
#'    \code{frac_pow(x[i], eq$power_numers[i], eq$power_denoms[i], eq$abs) * eq$coeffs[i]} is compared to \code{quota}).
#' Deal with components with \code{exp()} term first, and generate each coordinate while fulfilling \code{quota} if possible; if not, randomly generate from \cr
#'    \code{[-0.01,0.01]/abs(eq$power_numers[i])}.
#' Then recalculate the new \code{quota} which subtracts the exp() terms from \code{eq$const}, and this time divided by the number of remaining components.
#' If \code{quota} becomes negative and \code{eq$larger == FALSE}, each component, after \code{frac_pow()} is assumed to give a negative number.
#' This is not possible if the term has the form x^{even_number/even_number}, or if the term is not log() in the case where \code{eq$nonnegative == TRUE || eq$abs == TRUE}.
#' Change quota to a positive smaller in absolute value for these bad terms and generate.
#' Finally, recalculate quota as before and generate the rest of the "good" components.
#' 
#' In some extreme domains the function may fail to generate a point within the domain.
#' Also, it is not guaranteed that the function returns a point in an area with a high probability density. 
#' @return A \code{p} vector inside the domain defined by \code{domain}.
#' @examples
#' p <- 30
#' poly_d <- function(ex, abs, nng){
#'    return (make_domain("polynomial", p=p, 
#'                        ineqs=list(list(expression=ex, abs=abs, nonnegative=nng))))
#' }
#' 
#' random_init_polynomial(poly_d(paste("sum(exp(x))<=", p*1.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("sum(exp(x))<=", p*1.01), abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d(paste("sum(exp(x))>", p*1.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("sum(exp(x))>", p*1.01), abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d(paste("sum(log(x))<=", 0.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("sum(log(x))>", 0.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("sum(x^2)<=", 0.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("sum(x^2)>", 0.01), abs=TRUE, nng=TRUE))
#' 
#' random_init_polynomial(poly_d(paste("exp(x)<=", 1.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("exp(x)<=", 1.01), abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d(paste("exp(x)>", 1.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("exp(x)>", 1.01), abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d(paste("log(x)<=", 0.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("log(x)>", 0.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("x^2<=", 0.01), abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(paste("x^2>", 0.01), abs=TRUE, nng=TRUE))
#' 
#' random_init_polynomial(poly_d("x1^2+x2^2+log(x3)<-2", abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d("x1^2+x2^2+log(x3)>-2", abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+x2^2+x3^(1/3)<-2", abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+x2^2+x3^(1/3)>-2", abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+1.2*exp(2*x2)+2.3*exp(-3*x3)<-2", abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+1.2*exp(2*x2)+2.3*exp(-3*x3)<2", abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+1.2*exp(2*x2)+2.3*exp(-3*x3)>-2", abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+2.3*log(x4)+1.3*exp(2*x2)+0.7*exp(-3*x3)<-2", 
#'                        abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d("x1^(3/5)+2.3*log(x4)+1.3*exp(2*x2)+0.7*exp(-3*x3)>-2", 
#'                        abs=FALSE, nng=FALSE))
#' random_init_polynomial(poly_d(
#'    "x1^(3/5)+0.9*x2^(2/3)+2.7*x3^(-3/2)+1.1*x4^(-5)+1.1*exp(2x5)+1.3*exp(-3x6)+0.7*log(x7)<-2", 
#'                        abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d(
#'    "x1^(3/5)+0.9*x2^(2/3)+2.7*x3^(-3/2)+1.1*x4^(-5)+1.1*exp(2x5)+1.3*exp(-3x6)+0.7*log(x7)<-2", 
#'                        abs=FALSE, nng=TRUE))
#' random_init_polynomial(poly_d(
#'    "x1^(3/5)+0.9*x2^(2/3)+2.7*x3^(-3/2)+1.1*x4^(-5)+1.1*exp(2x5)+1.3*exp(-3x6)+0.7*log(x7)>-2", 
#'                        abs=TRUE, nng=FALSE))
#' random_init_polynomial(poly_d(
#'    "x1^(3/5)+0.9*x2^(2/3)+2.7*x3^(-3/2)+1.1*x4^(-5)+1.1*exp(2x5)+1.3*exp(-3x6)+0.7*log(x7)>2", 
#'                        abs=TRUE, nng=TRUE))
#' random_init_polynomial(poly_d(
#'    "x1^(3/5)+0.9*x2^(2/3)+2.7*x3^(-3/2)+1.1*x4^(-5)+1.1*exp(2x5)+1.3*exp(-3x6)+0.7*log(x7)>2", 
#'                        abs=FALSE, nng=FALSE))
#' @useDynLib genscore, .registration = TRUE
#' @export
random_init_polynomial <- function(domain) {
  if (!"checked" %in% names(domain) || domain$type != "polynomial")
    stop("domain must be an object returned by make_domain() with type == 'polynomial'.")

  if (length(domain$ineqs) > 1)
    stop("Random initialization for domains defined by more than 1 inequality not ",
         "supported. Please provide with your own initialization.")
  eq <- domain$ineqs[[1]]
  tol <- 1e-10
  
  get_bound_and_runif <- function(n, eq, a, b, c){
    bounds <- .C("poly_domain_1d_for_R", a=as.integer(a),
                 b=as.integer(b), c=as.double(c),
                 larger=as.integer(eq$larger), abs=as.integer(eq$abs),
                 nonnegative=as.integer(eq$nonnegative),
                 num_intervals=as.integer(0),
                 lefts=as.double(c(0,0)), rights=as.double(c(0,0)),
                 print=as.integer(FALSE))
    if (bounds$num_intervals == 0)
      stop("Domain is empty. Please check your boundary conditions. If you believe this is an error, please provide your own xinit.")
    return (random_init_uniform(n, bounds$lefts[1:bounds$num_intervals],
                                bounds$rights[1:bounds$num_intervals]))
  }
  
  if (eq$uniform) {
    xinit <- get_bound_and_runif(domain$p, eq, eq$power_numers, eq$power_denoms, eq$const)
    if (!in_bound(xinit, domain))
      stop("Failed to initialize an x inside the boundary. Please provide your own initialization. Generated: ", paste(xinit, collapse=", "), ".")
    return (xinit)
  }
  if (any(eq$coeffs < 0))
    stop("Random initialization for polynomial domains with negative coefficients not ",
         "supported. Please provide with your own initialization.")
  xinit <- numeric(domain$p)
  bad_idx <- which(abs(eq$coeffs) < tol)
  xinit[bad_idx] <- rlaplace_truncated_centered(length(bad_idx), 
                                                 ifelse(eq$nonnegative, 0, -Inf), Inf)
  good_idx <- setdiff(1:domain$p, bad_idx)
  if (length(good_idx) == 0)
    return (xinit)
  quota <- eq$const / length(good_idx)
  for (i in 1:length(eq$power_numers)) { # Make power numers and denoms coprime
    tmp <- makecoprime(eq$power_numers[i], eq$power_denoms[i])
    eq$power_numers[i] <- tmp[1]; eq$power_denoms[i] <- tmp[2]
  }
  if (length(eq$power_denoms) == 1) { # If uniform power and feasible
    # Assumes coeffs[i] > 0
    xinit[good_idx] <- sapply(good_idx, function(i){
      get_bound_and_runif(1, eq, eq$power_numers, eq$power_denoms, quota / eq$coeffs[i])})
    if (!in_bound(xinit, domain))
      stop("Failed to initialize an x inside the boundary. Please provide your own initialization. Generated: ", paste(xinit, collapse=", "), ".")
    return (xinit)
  }
  bad_idx <- intersect(good_idx, which(eq$power_numers != 0 & eq$power_denoms == 0)) # indices with exp()
  if (length(bad_idx) > 0) { ### Deal with indices with exp() first
    if (quota <= 0) { # If exp(...)*coeffs compared to quota <= 0
      if (eq$larger == FALSE)  # Infeasible: exp() always > 0, so cannot satisfy exp() < quota; randomly generate from [-0.01,0.01]/power_numer
        xinit[bad_idx] <- sapply(abs(eq$power_numers[bad_idx]), 
                                 function(co){stats::runif(1, ifelse(eq$nonnegative, 0, -0.01/co), 0.01/co)})
      else # Feasible: exp() always > 0, so can always satisfy exp() > quota; randomly generate laplace
        xinit[bad_idx] <- rlaplace_truncated_centered(length(bad_idx),
                                                      ifelse(eq$nonnegative, 0, -Inf), Inf) / abs(eq$power_numers[bad_idx])
    } else { # If exp(...)*coeffs compared to quota > 0
      if ((!eq$nonnegative) && (!eq$abs)) # Feeasible: If not nonneg and not abs, exp() can be in (0, +Inf) and so quota > 0 can be satisfied
        xinit[bad_idx] <- sapply(bad_idx, function(i){
          get_bound_and_runif(1, eq, eq$power_numers[i], 0, quota / eq$coeffs[i])})
      else { # If exp(power_numers*|x|) * coeffs compared to quota > 0
        ind_bounds <- quota / eq$coeffs # exp(power_numers[i]*|x|) compared to ind_bounds[i]
        if (eq$larger) # exp() > ind_bounds feasible for power_numers > 0 or ind_bounds < 1
          bad_idx_feas <- intersect(bad_idx, which(eq$power_numers > 0 | ind_bounds < 1))
        else # exp() < ind_bounds feasible for power_numers < 0 or ind_bounds > 1
          bad_idx_feas <- intersect(bad_idx, which(eq$power_numers < 0 | ind_bounds > 1))
        bad_idx_infeas <- setdiff(bad_idx, bad_idx_feas)
        if (length(bad_idx_feas) > 0)
          xinit[bad_idx_feas] <- sapply(bad_idx_feas, function(i){ # Generate within quota
          get_bound_and_runif(1, eq, eq$power_numers[i], 0, ind_bounds[i])})
        if (length(bad_idx_infeas) > 0)
          xinit[bad_idx_infeas] <- sapply(abs(eq$power_numers[bad_idx_infeas]), # Infeasible; generate from [-0.01,0.01]/power_numer
                                        function(co){stats::runif(1, ifelse(eq$nonnegative, 0, -0.01/co), 0.01/co)})
      }
    }
    good_idx <- setdiff(good_idx, bad_idx)
    quota <- (eq$const - sum(eq$coeffs[-good_idx] *
                               frac_pow(xinit[-good_idx], eq$power_numers[-good_idx],
                                        eq$power_denoms[-good_idx], eq$abs))) / length(good_idx)
  }
  if (length(good_idx) > 0 && quota <= 0 && eq$larger == FALSE) { ### If quota <= 0 and !eq$larger, deal with terms that cannot produce a negative number first.
    if (eq$nonnegative || eq$abs) { # If nonnegative and abs, only terms log(|x|) can be negative, since |x|^... or exp(|x|) cannot be <= eq$const < 0.
      bad_idx <- intersect(good_idx, which((eq$power_numers != 0 | eq$power_denoms != 0)))
    } else 
      # If eq$nonnegative and eq$abs both FALSE, indices that can produce negative numbers after frac_pow are those with odd power numers and denoms, or log()
      bad_idx <- intersect(good_idx, 
                           which((eq$power_numers %% 2 == 0 | eq$power_denoms %% 2 == 0) & 
                             (eq$power_numers != 0 | eq$power_denoms != 0)))

   if (length(bad_idx) > 0) {
      good_idx <- setdiff(good_idx, bad_idx)
      if (length(good_idx) == 0)
        stop("Domain is empty. Please check your boundary conditions. If you believe this is an error, please provide your own xinit.")
      # Set quota for bad indices to a smaller value
      quota_bad <- abs(quota) * max(1/2, length(good_idx) / length(bad_idx))
      xinit[bad_idx] <- sapply(bad_idx, function(i){
        get_bound_and_runif(1, eq, eq$power_numers[i], eq$power_denoms[i], quota_bad / eq$coeffs[i])})
      quota <- (eq$const - sum(eq$coeffs[-good_idx] * 
                                 frac_pow(xinit[-good_idx], eq$power_numers[-good_idx],
                                       eq$power_denoms[-good_idx], eq$abs))) / length(good_idx)
    }
  }
  if (length(good_idx) > 0)
    xinit[good_idx] <- sapply(good_idx, function(i){
    get_bound_and_runif(1, eq, eq$power_numers[i], eq$power_denoms[i], quota / eq$coeffs[i])})
  if (!in_bound(xinit, domain))
    stop("Failed to initialize an x inside the boundary. Please provide your own initialization. Generated: ", paste(xinit, collapse=", "), ".")
  return (xinit)
}
