## 1. Deal with potential error in rab_arms: Retry
## 4. Allow fractional powers for negative numbers in density calculations/estimation
## 5. scaling for x is not appropriate for boundaries that are not scale-invariant

#dyn.load("src/genscore.so")

#' Finds the greatest (positive) common divisor of two integers.
#'
#' Finds the greatest (positive) common divisor of two integers; if one of them is 0, returns the absolute value of the other number.
#'
#' @param a An integer.
#' @param b An integer.
#' @return The greatest (positive) common divisor of two integers; if one of them is 0, returns the absolute value of the other number.
#' @examples
#' gcd(1, 2)
#' gcd(1, -2)
#' gcd(12, -18)
#' gcd(-12, 18)
#' gcd(15, 0)
#' gcd(0, -15)
#' gcd(0, 0)
#' @export
gcd <- function(a, b) {
  a <- as.integer(a); b <- as.integer(b)
  if (b == 0) return (abs(a))
  return (gcd(b, a %% b))
}

#' Makes two integers coprime.
#'
#' Divides both integers by their greatest common divisor, switching their signs if the second integer is negative. If either integer is 0, returns without modification.
#'
#' @param a An integer.
#' @param b An integer.
#' @return The greatest (positive) common divisor of two integers; if one of them is 0, returns the absolute value of the other number.
#' @examples
#' makecoprime(1, 2)
#' makecoprime(1, -2)
#' makecoprime(12, -18)
#' makecoprime(-12, 18)
#' makecoprime(15, 0)
#' makecoprime(0, -15)
#' makecoprime(0, 0)
#' @export
makecoprime <- function(a, b) {
  if (a == 0 || b == 0) return (c(a, b))
  tmp <- gcd(a,b); a <- a/tmp; b <- b/tmp
  if (b < 0) {return (c(-a, -b))}
  return (c(a, b))
}

#' Evaluate x^(a/b) and |x|^(a/b) with integer a and b with extension to conventional operations.
#'
#' Evaluate x^(a/b) and |x|^(a/b) with integer a and b with extension to conventional operations (listed under details) that would otherwise result in \code{NaN}.
#'
#' @param x A number or a vector of numbers.
#' @param a An integer.
#' @param b An integer.
#' @param abs TRUE or FALSE.
#' @return A vector of numbers of the same size as \code{x}. See details.
#' @details
#' Replace \code{x} by \code{abs(x)} below if \code{abs == TRUE}.
#' If \code{a == 0 && b == 0}, returns \code{log(x)}.
#' If \code{a != 0 && b == 0}, returns \code{exp(a*x)}.
#' Otherwise, for \code{b != 0}, evaluates \code{x^(a/b)} with the following extensions.
#' \code{0^0} evaluates to \code{1}.
#' If \code{x < 0}, returns \code{(-1)^a * |x|^(a/b)} if \code{b} is odd, or \code{NaN} otherwise.
#' If \code{x == 0 && a < 0}, returns \code{NaN}.
#' @examples
#' frac_pow(-5:5, 3, 2, TRUE) - abs(-5:5)^(3/2)
#' frac_pow(-5:5, 5, 3, FALSE) - sign(-5:5)^5*abs(-5:5)^(5/3)
#' frac_pow(-5:5, 2, 3, FALSE) - ((-5:5)^2)^(1/3)
#' frac_pow(c(-5:-1,1:5), 0, 0, TRUE) - log(abs(c(-5:-1,1:5)))
#' frac_pow(-5:5, 0, 1, FALSE) - 1
#' frac_pow(-5:5, 3, 0, FALSE) - exp(3*-5:5)
#' @export
frac_pow <- Vectorize(function(x, a, b, abs){
  if (b == 0) {
    if (a == 0) {
      if (abs) return (log(abs(x)))
      if (x > 0) return (log(x))
      stop("x^(0/0) treated as log(x) which is undefined with x < 0. Got ", x, ".")
    } else {
      if (abs) return (exp(a*abs(x)))
      else return (exp(a*x))
    }
  }
  if (a == 0) return (1)
  tmp <- makecoprime(a, b); a <- tmp[1]; b <- tmp[2]
  if (x < 0) {
    if (abs) return ((-x)^(a/b))
    if (b %% 2) {return ((-1)^a*(-x)^(a/b))}
    stop("x^(a/b) is undefined with x < 0 and b even. Got ", x, ".")
  }
  if (x == 0 && a / b < 0) return (NaN)
  return (x^(a/b))
})

#' Checks if two equally sized numeric vectors satisfy the requirements for being left and right endpoints of a domain defined as a union of intervals.
#'
#' Checks if two equally sized numeric vectors satisfy the requirements for being left and right endpoints of a domain defined as a union of intervals.
#'
#' @param lefts A non-empty vector of numbers (may contain \code{-Inf}), the left endpoints of a domain defined as a union of intervals.
#' @param rights A non-empty vector of numbers (may contain \code{Inf}), the right endpoints of a domain defined as a union of intervals. Must have the same size as \code{lefts}.
#' @return \code{NULL}. Program stops if \code{lefts} and \code{rights} do not define valid left and right endpoints.
#' @details
#' Both \code{lefts} and \code{rights} must be non-empty and should have the same length.
#' Suppose \code{lefts} and \code{rights} both have length l, [lefts[1], rights[1]], ..., [lefts[l], rights[l]] must be an increasing and non-overlapping set of valid intervals, meaning lefts[i] <= rights[i] <= lefts[j] for any i < j (singletons and overlapping at the boundary points are allowed).
#' \code{Inf} is not allowed in \code{lefts} and \code{-Inf} is not allowed in \code{rights}.
#' @examples
#' ## [-4,-3], [-2,-1], [0,1], [2,3], [4,5]
#' check_endpoints(lefts=c(-4,-2,0,2,4), rights=c(-3,-1,1,3,5))
#' #check_endpoints(lefts=c(), rights=c()) # Cannot be empty
#' #check_endpoints(lefts=c(-4,-2,0,2,4), rights=c(-3,-1,1,3)) # Unequal length
#' #check_endpoints(lefts=c(Inf), rights=c(Inf)) # No Inf in lefts, otherwise invalid interval
#' #check_endpoints(lefts=c(-Inf), rights=c(-Inf)) # No -Inf in rights, otherwise invalid interval
#' #check_endpoints(lefts=c(0, 1), rights=c(2, 3)) # [0,2] and [1,3] overlap, not allowed
#' #check_endpoints(lefts=c(2, 0), rights=c(3, 1)) # [2,3], [0,1] not increasing, not allowed
#' ## Singletons and overlapping at the boundary points allowed
#' check_endpoints(lefts=c(0, 1, 2), rights=c(0, 2, 3))
#' @export
check_endpoints <- function(lefts, rights){
  if (length(lefts) == 0 || length(rights) == 0)
    stop("lefts and rights must not be empty.")
  if (length(lefts) != length(rights))
    stop("lefts and rights must have the same length.")
  if (any(is.infinite(lefts) & lefts > 0))
    stop("No positive infinity allowed in lefts.")
  if (any(is.infinite(rights) & rights < 0))
    stop("No negative infinity allowed in rights.")
  if (sum(is.infinite(lefts)) > 1 || sum(is.infinite(rights)) > 1)
    stop("Only one negative infinity is allowed in lefts and one positive infinity in rights to avoid overlapping intervals.")
  pts <- c(rbind(lefts, rights))
  if (any(pts[2:length(pts)]-pts[1:(length(pts)-1)] < 0))
    stop("Points in lefts and rights must be increasing and intervals must be non-overlapping.")
}

#' Maximum between finite_infinity and 10 times the max abs value of finite elements in \code{lefts} and \code{rights}.
#'
#' Maximum between \code{finite_infinity} and 10 times the max abs value of finite elements in \code{lefts} and \code{rights}.
#'
#' @param lefts A non-empty vector of numbers (may contain \code{-Inf}), the left endpoints of a domain defined as a union of intervals. Must pass \code{check_endpoints(lefts, rights)}.
#' @param rights A non-empty vector of numbers (may contain \code{Inf}), the right endpoints of a domain defined as a union of intervals. Must pass \code{check_endpoints(lefts, rights)}.
#' @param finite_infinity A finite positive number. \code{Inf} will be truncated to \code{finite_infinity} if applicable. See details.
#' @return A double, larger than or equal to \code{finite_infinity}.
#' @details
#' Since we assume that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}, the function takes the maximum between \code{finite_infinity} and 10 times the absolute values of \code{lefts[1]}, \code{lefts[length(lefts)]}, \code{rights[1]}, and \code{rights[length(rights)]}, if they are finite.
#' @examples
#' # Does not change sincee 1000 > 12 * 10
#' update_finite_infinity_for_uniform(c(-10,-5,0,5,9), c(-8,-3,2,7,12), 1000)
#' # Changed to 12 * 10
#' update_finite_infinity_for_uniform(c(-10,-5,0,5,9), c(-8,-3,2,7,12), 10)
#' # Changed to 12 * 10
#' update_finite_infinity_for_uniform(c(-Inf,-5,0,5,9), c(-8,-3,2,7,12), 10)
#' # Changed to 9 * 10
#' update_finite_infinity_for_uniform(c(-Inf,-5,0,5,9), c(-8,-3,2,7,Inf), 10)
#' @export
update_finite_infinity_for_uniform <- function(lefts, rights, finite_infinity){
  if (is.infinite(finite_infinity))
    stop("finite_infinity must be finite.")
  if (length(lefts) > 1)
    finite_infinity <- max(finite_infinity, 10 * abs(rights[1]), 10 * abs(lefts[length(lefts)]))
  if (is.finite(lefts[1]))
    finite_infinity <- max(finite_infinity, 10 * abs(lefts[1]))
  if (is.finite(rights[length(rights)]))
    finite_infinity <- max(finite_infinity, 10 * abs(rights[length(rights)]))
  return (finite_infinity)
}

#' Creates a list of elements that defines the domain for a multivariate distribution.
#'
#' Creates a list of elements that define the domain for a multivariate distribution.
#'
#' @param type A string, the domain type. Currently support \code{"R"}, \code{"R+"}, \code{"uniform"}, \code{"polynomial"}, \code{"simplex"}. See details.
#' @param p An integer, the dimension of the domain.
#' @param lefts Optional, required if \code{type == "uniform"} and must have the same length as \code{rights}. A non-empty vector of numbers (may contain \code{-Inf}), the left endpoints of a domain defined as a union of intervals. It is required that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}.
#' @param rights Optional, required if \code{type == "uniform"} and must have the same length as \code{lefts}. A non-empty vector of numbers (may contain \code{Inf}), the right endpoints of a domain defined as a union of intervals. It is required that \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}.
#' @param ineqs Optional, required if \code{type == "polynomial"}. A list of lists, each sublist representing an inequality that defines the domain. Each sublist must contain fields \code{abs} (logical) and \code{nonnegative} (logical), and in addition either a single \code{expression} (string), or all of the following: \code{uniform} (logical), \code{larger} (logical), \code{power_numers} (1 or \code{p} integers), \code{power_denoms} (1 or \code{p} integers), \code{const} (a number), \code{coeffs} (1 or \code{p} numbers).
#' @param rule Optional, required if \code{type == "polynomial" && length(ineqs) > 1}. A string containing inequality numbers, spaces, parentheses, '&' and '|' only. Used to indicate the logic operation on how to combine the domains defined by each inequality, i.e. "(1 & 2 && 3) || 4 | 5". Chained operations not separated by parentheses are only allowed for the same type of operation ('&'/'|'), i.e. "1 & 2 | 3" is not allowed; it should be either "(1 & 2) | 3" or "1 & (2 | 3)".
#' @return A list containing the elements that define the domain.
#' For all types of domains, the following are returned.
#'   \item{type}{A string, same as the input.}
#'   \item{p}{An integer, same as the input.}
#'   \item{p_deemed}{An integer, equal to \code{p-1} if \code{type == "simplex"} or \code{p} otherwise.}
#'   \item{checked}{A logical, \code{TRUE}. Used in other functions to test whether a list is returned by this function.}
#' In addition,
#' \itemize{
#'    \item{For \code{type == "simplex"}, returns in addition
#'      \describe{
#'         \item{\code{simplex_tol}}{Tolerance used for simplex domains. Defaults to \code{1e-10}.}
#'      }}
#'    \item{For \code{type == "uniform"}, returns in addition
#'      \describe{
#'        \item{\code{lefts}}{A non-empty vector of numbers, same as the input.}
#'        \item{\code{rights}}{A non-empty vector of numbers, same as the input.}
#'        \item{\code{left_inf}}{A logical, indicates whether \code{lefts[1]} is \code{-Inf}.}
#'        \item{\code{right_inf}}{A logical, indicates whether \code{rights[length(rights)]} is \code{Inf}.}
#'      }}
#'    \item{For \code{type == "polynomial"}, returns in addition
#'      \describe{
#'        \item{\code{rule}}{A string, same as the input if provided and valid; if not provided and \code{length(ineqs) == 1}, set to \code{"1"} by default.}
#'       \item{\code{postfix_rule}}{A string, \code{rule} in postfix notation (reverse Polish notation) containing numbers, \code{" "}, \code{"&"} and \code{"|"} only.}
#'       \item{\code{ineqs}}{A list of lists, each sublist representing one inequality containing the following fields:
#'         \describe{
#'         \item{\code{uniform}}{A logical, indicates whether the inequality should be uniformly applied to all components (e.g. \code{"x>1"} -> \code{"x1>1 && ... && xp>1"}).}
#'         \item{\code{larger}}{A logical, indicates whether the polynomial should be larger or smaller than the constant (e.g. \code{TRUE} for \code{x1 + ... + xp > C}, and \code{FALSE} for \code{x1 + ... + xp < C}).}
#'         \item{\code{const}}{A number, the constant the polynomial should be compared to  (e.g. \code{2.3} for \code{x1 + ... + xp > 2.3}).}
#'         \item{\code{abs}}{A logical, indicates whether one should evaluate the polynomial in \code{abs(x)} instead of \code{x}.}
#'         \item{\code{nonnegative}}{A logical, indicates whether the domain of this inequality should be restricted to the non-negative orthant.}
#'         \item{\code{power_numers}}{A single integer or a vector of \code{p} integers. \code{x[i]} will be raised to the power of \code{power_numers[i] / power_denoms[i]} (or without subscript if a singer integer). Note that \code{x^(0/0)} is interpreted as \code{log(x)}, and \code{x^(n/0)} as \code{exp(n*x)} for \code{n} non-zero. For a negative \code{x}, \code{x^(a/b)} is defined as \code{(-1)^a*|x|^(a/b)} if \code{b} is odd, or \code{NaN} otherwise.}
#'         \item{\code{power_denoms}}{A single integer or a vector of \code{p} integers.}
#'         \item{\code{coeffs}}{\code{NULL} if \code{uniform == TRUE}. A vector of \code{p} doubles, where \code{coeffs[i]} is the coefficient on \code{x[i]} in the inequality}
#'      }}}}}
#' @details The following types of domains are supported:
#' \describe{
#'   \item{\code{"R"}}{The entire \code{p}-dimensional real space. Equivalent to \code{"uniform"} type with \code{lefts=-Inf} and \code{rights=Inf}.}
#'   \item{\code{"R+"}}{The non-negative orthant of the \code{p}-dimensional real space. Equivalent to \code{"uniform"} type with \code{lefts=0} and \code{rights=Inf}.}
#'   \item{\code{"uniform"}}{A union of finitely many disjoint intervals as a uniform domain for all components. The left endpoints should be specified through \code{lefts} and the right endpoints through \code{rights}. The intervals must be disjoint and strictly increasing, i.e. \code{lefts[i] <= rights[i] <= lefts[j]} for any \code{i < j}. E.g. \code{lefts=c(0, 10)} and \code{rights=c(5, Inf)} represents the domain ([0,5]v[10,+Inf])^p.}
#'   \item{\code{"simplex"}}{The standard \code{p-1}-simplex with all components positive and sum to 1, i.e. \code{sum(x) == 1} and \code{x > 0}.}
#'   \item{\code{"polynomial"}}{A finite intersection/union of domains defined by comparing a constant to a polynomial with at most one term in each component and no interaction terms (e.g. \code{x1^3+x1^2>1} or \code{x1*x2>1} not supported). The following is supported: \code{{x1^2 + 2*x2^(3/2) > 1} && ({3.14*x1 - 0.7*x3^3 < 1} || {-exp(3*x2) + 3.7*log(x3) + 2.4*x4^(-3/2)})}.}
#'      To specify a polynomial-type domain, one should define the \code{ineqs}, and in case of more than one inequality, the logical \code{rule} to combine the domains defined by each inequality.
#'      \item{\code{rule}}{A logical rule in infix notation, e.g. \code{"(1 && 2 && 3) || (4 && 5) || 6"}, where the numbers represent the inequality numbers starting from 1. \code{"&&"} and \code{"&"} are not differentiated, and similarly for \code{"||"} and \code{"|"}. Chained operations are only allowed for the same operation (\code{"&"} or \code{"|"}), so instead of \code{"1 && 2 || 3"} one should write either \code{"(1 && 2) || 3"} or \code{"1 && (2 || 3)"} to avoid ambiguity.}
#'      \item{\code{ineqs}}{A list of lists, each sublist represents one inequality, and must contain the following fields:
#'         \describe{
#'           \item{\code{abs}}{A logical, indicates whether one should evaluate the polynomial in \code{abs(x)} instead of \code{x} (e.g. \code{"sum(x) > 1"} with \code{abs == TRUE} is interpreted as \code{sum(abs(x)) > 1}).}
#'           \item{\code{nonnegative}}{A logical, indicates whether the domain of this inequality should be restricted to the non-negative orthant.}
#'         }
#'
#'         In addition, one must in addition specify either a single string \code{"expression"} (highly recommended, detailed below), or all of the following fields (discouraged usage):
#'         \describe{
#'           \item{\code{uniform}}{A logical, indicates whether the inequality should be uniformly applied to all components (e.g. \code{"x>1"} -> \code{"x1>1 && ... && xp>1"}).}
#'           \item{\code{larger}}{A logical, indicates whether the polynomial should be larger or smaller than the constant (e.g. \code{TRUE} for \code{x1 + ... + xp > C}, and \code{FALSE} for \code{x1 + ... + xp < C}).}
#'           \item{\code{const}}{A number, the constant the polynomial should be compared to  (e.g. \code{2.3} for \code{x1 + ... + xp > 2.3}).}
#'           \item{\code{power_numers}}{A single integer or a vector of \code{p} integers. \code{x[i]} will be raised to the power of \code{power_numers[i] / power_denoms[i]} (or without subscript if a singer integer). Note that \code{x^(0/0)} is interpreted as \code{log(x)}, and \code{x^(n/0)} as \code{exp(n*x)} for \code{n} non-zero. For a negative \code{x}, \code{x^(a/b)} is defined as \code{(-1)^a*|x|^(a/b)} if \code{b} is odd, or \code{NaN} otherwise. (Example: \code{c(2,3,5,0,-2)} for \code{x1^2+2*x2^(3/2)+3*x3^(5/3)+4*log(x4)+5*exp(-2*x)>1}).}
#'           \item{\code{power_denoms}}{A single integer or a vector of \code{p} integers.  (Example: \code{c(1,2,3,0,0)} for \code{x1^2+2*x2^(3/2)+3*x3^(5/3)+4*log(x4)+5*exp(-2*x)>1}).}
#'           \item{\code{coeffs}}{Required if \code{uniform == FALSE}. A vector of \code{p} doubles, where \code{coeffs[i]} is the coefficient on \code{x[i]} in the inequality.}
#'        }
#'        The user is recommended to use a single \code{expression} for ease and to avoid potential errors. The user may safely skip the explanations and directly look at the examples to get a better understanding.\cr
#'
#'        The expression should have the form \code{"POLYNOMIAL SIGN CONST"}, where \code{"SIGN"} is one of \code{"<"}, \code{"<="}, \code{">"}, \code{">="}, and \code{"CONST"} is a single number (scientific notation allowed).\cr
#'
#'        \code{"POLYNOMIAL"} must be (1) a single term (see below) in \code{"x"} with no coefficient (e.g. \code{"x^(2/3)"}, \code{"exp(3x)"}), or (2) such a term surrounded by \code{"sum()"} (e.g. \code{"sum(x^(2/3))"}, \code{"sum(exp(3x))"}), or (3) a sum of such terms in \code{"x1"} through \code{"xp"} (one term max for each component) with or without coefficients, separated by the plus or the minus sign (e.g. \cr\code{"2.3x1^(2/3)-3.4exp(x2)+x3^(-3/5)"}).\cr
#'
#'        For (1) and (2), the term must be in one of the following forms: \code{"x"}, \code{"x^2"}, \code{"x^(-2)"}, \code{"x^(2/3)"}, \code{"x^(-2/3)"}, \code{"log(x)"}, \code{"exp(x)"}, \code{"exp(2x)"}, \code{"exp(2*x)"}, \code{"exp(-3x)"}, where the \code{2} and \code{3} can be changed to any other non-zero integers.\cr
#'        For (3), each term should be as above but in \code{"x1"}, ..., \code{"xp"} instead of \code{"x"}, following an optional double number and optionally a \code{"*"} sign.\cr
#'
#'        Examples: For \code{p=10}, \cr
#'           (1) \code{"x^2 > 2"} defines the domain \code{abs(x1) > sqrt(2) && ... && abs(x10) > sqrt(2)}.\cr
#'           (2) \code{"sum(x^2) > 2"} defines the domain \code{x1^2 + ... + x10^2 > 2}.\cr
#'           (3) \code{"2.3x3^(2/3)-3.4x4+x5^(-3/5)+3.7*x6^(-4)-1.9*log(x7)+1.3e5*exp(-3x8)}\cr
#'           \code{-2*exp(x9)+0.5exp(2*x10) <= 2"} defines a domain using a polynomial in \code{x3} through \code{x10}, and \code{x1} and \code{x2} are thus allowed to vary freely.\cr
#'
#'        Note that \code{">"} and \code{">="} are not differentiated, and so are \code{"<"} and \code{"<="}.
#'    }}
#' @examples
#' p <- 30
#' # The 30-dimensional real space R^30
#' domain <- make_domain("R", p=p)
#'
#' # The non-negative orthant of the 30-dimensional real space, R+^30
#' domain <- make_domain("R+", p=p)
#'
#' # x such that sum(x^2) > 10 && sum(x^(1/3)) > 10 with x allowed to be negative
#' domain <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list("expression"="sum(x^2)>10", abs=FALSE, nonnegative=FALSE),
#'                       list("expression"="sum(x^(1/3))>10", abs=FALSE, nonnegative=FALSE)))
#' # Alternatively
#' domain2 <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list(uniform=FALSE, power_numers=2, power_denoms=1, const=10, coeffs=1,
#'                  larger=1, abs=FALSE, nonnegative=FALSE),
#'                  list(uniform=FALSE, power_numers=1, power_denoms=3, const=10, coeffs=1,
#'                  larger=1, abs=FALSE, nonnegative=FALSE)))
#'
#'
#' # ([0, 1] v [2,3]) ^ p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#'
#' # x such that {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
#' # Alternatively
#' domain2 <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list(uniform=FALSE, power_numers=1, power_denoms=1, const=1,
#'                  coeffs=c(1,rep(0,p-1)), larger=1, abs=FALSE, nonnegative=TRUE),
#'                  list(uniform=FALSE, power_numers=1, power_denoms=1, const=1,
#'                  coeffs=c(0,1,rep(0,p-2)), larger=0, abs=FALSE, nonnegative=TRUE),
#'                  list(uniform=TRUE, power_numers=1, power_denoms=0, const=1.3,
#'                  larger=1, abs=FALSE, nonnegative=FALSE)))
#'
#'
#' # x in R_+^p such that {sum(log(x))<2 || (x1^(2/3)-1.3x2^(-3)<1 && exp(x1)+2.3*x2>2)}
#' domain <- make_domain("polynomial", p=p, rule="1 || (2 && 3)",
#'        ineqs=list(list("expression"="sum(log(x))<2", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x1^(2/3)-1.3x2^(-3)<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x1)+2.3*x2^2>2", abs=FALSE, nonnegative=TRUE)))
#' # Alternatively
#' domain2 <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list(uniform=FALSE, power_numers=0, power_denoms=0, const=2,
#'                  coeffs=1, larger=0, abs=FALSE, nonnegative=TRUE),
#'                  list(uniform=FALSE, power_numers=c(2,-3,rep(1,p-2)), power_denoms=c(3,rep(1,p-1)),
#'                  const=1, coeffs=c(1.0,-1.3,rep(0,p-2)), larger=0, abs=FALSE, nonnegative=TRUE),
#'                  list(uniform=FALSE, power_numers=c(1,2,rep(1,p-2)), power_denoms=c(0,rep(1,p-1)),
#'                  const=2, coeffs=c(1,2.3,rep(0,p-2)), larger=1, abs=FALSE, nonnegative=TRUE)))
#'
#'
#' # x in R_+^p such that {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' # Alternatively
#' domain2 <- make_domain("polynomial", p=p,
#'        ineqs=list(list(uniform=FALSE, power_numers=1, power_denoms=1, const=1,
#'                  coeffs=1:p, larger=0, abs=FALSE, nonnegative=TRUE)))
#'
#'
#' # The (p-1)-simplex
#' domain <- make_domain("simplex", p=p)
#'
#' # The l-1 ball {sum(|x|) < 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"="sum(x)<1", abs=TRUE, nonnegative=FALSE)))
#' @export
make_domain <- function(type, p, lefts=NULL, rights=NULL, ineqs=NULL, rule=NULL) {
  if (type == "R" || type == "R+") {
    return (list("type"=type, "p"=p, "p_deemed"=p, "checked"=TRUE))
  } else if (type == "simplex") {
    return (list("type"=type, "p"=p, "p_deemed"=p-1, "checked"=TRUE, "simplex_tol"=1e-10))
  } else if (type == "uniform") {
    if (is.null(lefts) || is.null(rights))
      stop("You must specify lefts and rights for uniform-type domain.")
    check_endpoints(lefts, rights)
    if (length(lefts) == 1 && rights == Inf) {
      if (lefts == -Inf) {
        warning("Domain type automatically changed to R.")
        return (make_domain("R", p))
      } else if (lefts == 0) {
        warning("Domain type automatically changed to R+.")
        return (make_domain("R+", p))
      }
    }
    return (list("type"=type, "p"=p, "p_deemed"=p, "lefts"=lefts, "rights"=rights,
                 "left_inf"=is.infinite(lefts[1]), "right_inf"=is.infinite(rights[length(rights)]),
                 "checked"=TRUE))
  } else if (type == "polynomial") {
    if (is.null(ineqs))
      stop("You must specify \"ineqs\" (inequalities) for polynomial-type domain.")
    if (!is.list(ineqs) || !all(sapply(ineqs, is.list)) || length(ineqs) == 0)
      stop("\"ineqs\" must be a list of at least one list(s) (inequalities).")
    domain <- list("type"=type, "p"=p, "p_deemed"=p, "checked"=TRUE)

    if (is.null(rule)) {
      if (length(ineqs) > 1)
        stop("You must specify the rule to combine the regions defined by polynomials.")
      else
        rule <- "1"
    } else if ((!is.character(rule)) || length(rule) != 1)
      stop("rule must be a single string.")
    domain$rule <- rule
    domain$ineqs <- ineqs
    domain$postfix_rule <- get_postfix_rule(rule, length(ineqs))

    # Equation-specific
    for (ei in 1:length(ineqs)){
      ineq <- ineqs[[ei]]
      for (field in c("abs", "nonnegative")) {
        if (!field %in% names(ineq))
          stop("Inequality ", ei, " must have element '", field, "'.")
        if (!ineq[[field]] %in% c(0,1))
          stop("In inequality ", ei, ", ", field, " must be 0 (FALSE) or 1 (TRUE).")
      }
      if ("expression" %in% names(ineq)) {
        # If ineq provided as an expression, parse it
        domain$ineqs[[ei]] <- c(parse_ineq(ineq$expression, p),
                                    list("abs"=ineq$abs, "nonnegative"=ineq$nonnegative))
      } else {
        # If ineq elements are directly specifieed
        for (field in c("uniform", "larger", "power_numers", "power_denoms", "const"))
          if (!field %in% names(ineq))
            stop("If ineq expression not directly specified by ineq$expression, inequality ", ei, " must have element '", field, "'.")

        if (!ineq$larger %in% c(0,1))
          stop("In ineq ", ei, " 'larger' must be 0 (FALSE) or 1 (TRUE).")
        if (!ineq$uniform %in% c(0,1))
          stop("In ineq ", ei, " 'uniform' must be 0 (FALSE) or 1 (TRUE).")

        if (ineq$uniform) {
          # For uniform inequalities, power_numers and power_denoms must be scalars
          if (length(ineq$power_numers) != 1 || length(ineq$power_denoms) != 1)
            stop("Since inequality ", ei, " is uniform, its power_numers and power_denoms must both be a single number.")
        } else {
          # For non-uniform ineqs, power_numers and power_denoms must have lengths 1 or p
          if (length(ineq$power_numers) == 1 && length(ineq$power_denoms) == p) {
            # If one has length 1 and the other has length p, duplicate so that both have length p
            ineq$power_numers <- rep(ineq$power_numers, p)
            domain$ineqs[[ei]]$power_numers <- ineq$power_numers
          }
          else if (length(ineq$power_numers) == p && length(ineq$power_denoms) == 1) {
            # If one has length 1 and the other has length p, duplicate so that both have length p
            ineq$power_denoms <- rep(ineq$power_denoms, p)
            domain$ineqs[[ei]]$power_denoms <- ineq$power_denoms
          }
          else if ((length(ineq$power_numers) != 1 && length(ineq$power_numers) != p) ||
                   (length(ineq$power_denoms) != 1 && length(ineq$power_denoms) != p))
            stop("In inequality ", ei, " power_numers and power_denoms must have length 1 or p.")

          # For non-uniform ineqs, coefficients for each component must be specified.
          if (!"coeffs" %in% names(ineq))
            stop("If ineq expression not directly specified by ineq$expression and is also not uniform, inequality ", ei, " must have element 'coeffs'.")
          if (length(ineq$coeffs) == 1)
            domain$ineqs[[ei]]$coeffs <- rep(ineq$coeffs, p)
          else if (length(ineq$coeffs) != p)
            stop("In inequality ", ei, " coeffs must have length 1 or p.")
          if (any(ineqs$coeffs == 0))
            stop("Coefficient 0 not allowed.")
        }

        if (any(c(ineq$power_numers, ineq$power_denoms) %% 1 != 0))
          stop("In inequality ", ei, " power_numers and domain$power_denoms must all be integers.")
        if (any(ineq$power_denoms == 0)) {
          if (any(ineq$power_numers == 0 & ineq$power_denoms == 0))
            warning("In inequality ", ei, ": As a reminder, x^(0/0) will be treated as log(x).")
          if (any(ineq$power_numers != 0 & ineq$power_denoms == 0))
            warning("In ineq ", ei, ": As a reminder, x^(1/0) will be treated as exp(x).")
        }
        if (any(ineq$power_numers == 0 & ineq$power_denoms != 0))
          warning("In ineq ", ei, ": x^(0/n) with n != 0 will always return 1. If this is not intended, simply set the corresponding coefficients to 0 to exclude variables from the inequality.")
        if (length(ineq$const) != 1)
          stop("In ineq ", ei, " const must be a scalar.")
        if (ineq$const <= 0 && ineq$larger == FALSE) {
          if (ineq$nonnegative || ineq$abs ||
              all(ineq$power_numers %% 2 == 0 |
                  ineq$power_denoms %% 2 == 0))
            stop("Domain for inequality ", ei, " is empty. Please check your boundary conditions.")
        }
      }
    }
    return (domain)
  } else
    stop("Only R, R+, uniform, polynomial and simplex domain types are supported.")
}

#' Returns the character at a position of a string.
#'
#' Returns the character at a position of a string.
#'
#' @param string A string.
#' @param position A positive number.
#' @return A character
#' @details
#' Calls \code{substr(string, position, position)}.
#' @examples
#' s_at("123", 1)
#' s_at("123", 2)
#' s_at("123", 3)
#' s_at("123", 4)
#' s_at("123", 0)
#' @export
s_at <- function(string, position) {substr(string, position, position)}

#' Replaces consecutive "&"s and "|"s in a string to a single & and |.
#'
#' Replaces consecutive \code{"&"}s and \code{"|"}s in a string to a single \code{"&"} and \code{"|"}.
#'
#' @param rule A string containing positive integers, parentheses, and \code{"&"} and \code{"|"} only.
#' @details Applied to \code{domain$rule} if \code{domain$type == "polynomial"}.
#' @return A string with extra \code{"&"}s and \code{"|"}s removed.
#' @examples
#' beautify_rule("(1 & 2 && 3 &&& 4) | 5 || 6 ||| 7")
#' @export
beautify_rule <- function(rule) {
  # (1 || 2) && 3 -> (1 | 2) & 3
  base::gsub("&+", "&", base::gsub("[\\|]+", "\\|", rule))
}

#' Changes a logical expression in infix notation to postfix notation using the shunting-yard algorithm.
#'
#' Changes a logical expression in infix notation to postfix notation using the shunting-yard algorithm.
#'
#' @param rule A string containing positive integers, parentheses, and \code{"&"} and \code{"|"} only. \code{"&&"} and \code{"&"} are not differentiated, and similarly for \code{"||"} and \code{"|"}. Chained operations are only allowed for the same operation (\code{"&"} or \code{"|"}), so instead of \code{"1 && 2 || 3"} one should write either \code{"(1 && 2) || 3"} or \code{"1 && (2 || 3)"} to avoid ambiguity.
#' @param num_eqs An integer, must be larger than or equal to the largest integer appearing in \code{rule}.
#' @details Applied to \code{domain$rule} if \code{domain$type == "polynomial"}, and internally calls \code{beautify_rule()}.
#' @return \code{rule} in postfix notation.
#' @examples
#' get_postfix_rule("1 & 2 && 3", 3)
#' get_postfix_rule("1 & (2 || 3)", 3)
#' get_postfix_rule("(1 & 2) || 3 | (4 & (5 || 6) && 7) | 8 | (9 && (10 || 11 | 12) & 13)", 13)
#' #get_postfix_rule("1 && 2 & 3 && 4", 3) # Error, ineq number 4 appearing in \code{rule}.
#' # Error, ambigious rule. Change to either \code{"1 & (2 | 3)"} or \code{"(1 & 2) | 3"}.
#' #get_postfix_rule("1 & 2 | 3", 3)
#' @useDynLib genscore, .registration = TRUE
#' @export
get_postfix_rule <- function(rule, num_eqs) {
  tmp <- .C("shunting_yard", num_eqs=as.integer(num_eqs),
            infix_pt=beautify_rule(rule), postfix_pt="", errno=as.integer(0))
  if (tmp$errno)
    stop("Error occurred in C -> shunting_yard() called in get_postfix_rule().")
  return (tmp$postfix_pt)
}

#' Returns whether a vector or each row of a matrix falls inside a domain.
#'
#' Returns whether a vector or each row of a matrix falls inside a domain.
#'
#' @param x A vector of length or a matrix of number of columns equal to \code{domain$p} if \code{domain$type != "simplex"}, or either \code{domain$p} or \code{domain$p-1} otherwise.
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @details Returns whether a vector or each row of a matrix falls inside a domain.
#' If \code{domain$type == "simplex"}, if the length/number of columns is \code{domain$p}, returns \code{all(x > 0) && abs(sum(x) - 1) < domain$simplex_tol}; if the dimension is \code{domain$p-1}, returns \code{all(x > 0) && sum(x) < 1}.
#' @return A logical vector of length equal to the number of rows in \code{x} (\code{1} if \code{x} is a vector).
#' @examples
#' p <- 30
#' n <- 10
#'
#' # The 30-dimensional real space R^30, assuming probability of
#' domain <- make_domain("R", p=p)
#' in_bound(1:p, domain)
#' in_bound(matrix(1:(p*n), ncol=p), domain)
#'
#' # The non-negative orthant of the 30-dimensional real space, R+^30
#' domain <- make_domain("R+", p=p)
#' in_bound(matrix(1:(p*n), ncol=p), domain)
#' in_bound(matrix(1:(p*n) * (2*rbinom(p*n, 1, 0.98)-1), ncol=p), domain)
#'
#' # x such that sum(x^2) > 10 && sum(x^(1/3)) > 10 with x allowed to be negative
#' domain <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list("expression"="sum(x^2)>10", abs=FALSE, nonnegative=FALSE),
#'                       list("expression"="sum(x^(1/3))>10", abs=FALSE, nonnegative=FALSE)))
#' in_bound(rep((5/p)^3, p), domain)
#' in_bound(rep((10/p)^3, p), domain)
#' in_bound(rep((15/p)^3, p), domain)
#' in_bound(rep((5/p)^(1/2), p), domain)
#' in_bound(rep((10/p)^(1/2), p), domain)
#' in_bound(rep((15/p)^(1/2), p), domain)
#'
#' # ([0, 1] v [2,3]) ^ p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#' in_bound(c(0.5, 2.5)[rbinom(p, 1, 0.5)+1], domain)
#' in_bound(c(rep(0.5, p/2), rep(2.5, p/2)), domain)
#' in_bound(c(rep(0.5, p/2), rep(2.5, p/2-1), 4), domain)
#'
#' # x such that {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
#' in_bound(c(1.5, (log(1.3)+1)/2, rep(log(1.3)*2, p-2)), domain)
#' in_bound(c(0.5, (log(1.3)+1)/2, rep(log(1.3)*2, p-2)), domain)
#' in_bound(c(1.5, log(1.3)/2, rep(log(1.3)*2, p-2)), domain)
#' in_bound(c(1.5, (log(1.3)+1)/2, rep(log(1.3)/2, p-2)), domain)
#'
#' # x in R_+^p such that {sum(log(x))<2 || (x1^(2/3)-1.3x2^(-3)<1 && exp(x1)+2.3*x2>2)}
#' domain <- make_domain("polynomial", p=p, rule="1 || (2 && 3)",
#'        ineqs=list(list("expression"="sum(log(x))<2", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x1^(2/3)-1.3x2^(-3)<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x1)+2.3*x2^2>2", abs=FALSE, nonnegative=TRUE)))
#' in_bound(rep(exp(1/p), p), domain)
#' in_bound(c(1, 1, rep(1e5, p-2)), domain)
#'
#' # x in R_+^p such that {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' in_bound(0.5/p/1:p, domain)
#' in_bound(2/p/1:p, domain)
#' in_bound(rep(1/p, p), domain)
#' in_bound(rep(1/p^2, p), domain)
#'
#' # The (p-1)-simplex
#' domain <- make_domain("simplex", p=p)
#' x <- abs(matrix(rnorm(p*n), ncol=p))
#' x <- x / rowSums(x)
#' in_bound(x, domain) # TRUE
#' in_bound(x[,1:(p-1)], domain) # TRUE
#' x2 <- x
#' x2[,1] <- -x2[,1]
#' in_bound(x2, domain) # FALSE since the first component is now negative
#' in_bound(x2[,1:(p-1)], domain) # FALSE since the first component is now negative
#' x3 <- x
#' x3[,1] <- x3[,1] + domain$simplex_tol * 10
#' in_bound(x3, domain) # FALSE since the rows do not sum to 1
#' in_bound(x3[,1:(p-1)], domain) # TRUE since the first (p-1) elts in each row still sum to < 1
#' x3[,1] <- x3[,1] + x3[,p]
#' in_bound(x3[,1:(p-1)], domain) # FALSE since the first (p-1) elts in each row now sum to > 1
#'
#' # The l-1 ball {sum(|x|) < 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"="sum(x)<1", abs=TRUE, nonnegative=FALSE)))
#' in_bound(rep(0.5/p, p)*(2*rbinom(p, 1, 0.5)-1), domain)
#' in_bound(rep(1.5/p, p)*(2*rbinom(p, 1, 0.5)-1), domain)
#' @export
in_bound <- function(x, domain) { # x: vector or 2-dimensional matrix (n x p)
  if (!"checked" %in% names(domain))
    stop("domain must be an object returned by make_domain().")
  if (length(dim(x)) == 2)
    return (apply(x, 1, function(x){in_bound(x, domain)}))
  if (length(x) != domain$p && domain$type != "simplex")
    stop("x should have length ", domain$p, ".")
  if (domain$type == "R")
    return (TRUE)
  if (domain$type == "R+")
    return (all(x > 0))
  if (domain$type == "uniform")
    return (all(sapply(x, function(xx){(xx >= domain$lefts[1]) && (xx <= domain$rights[length(domain$rights)]) && (xx <= domain$rights[max(which(xx >= domain$lefts))])})))
  if (domain$type == "simplex") {
    if (length(x) == domain$p_deemed)
      return (all(x > 0) && sum(x) < 1)
    if (length(x) == domain$p)
      return (all(x > 0) && abs(sum(x) - 1) < domain$simplex_tol)
    stop("For simplex, x should have length ", domain$p, " or ", domain$p_deemed, " (excluding the last coordinate).")
  }
  if (domain$type == "polynomial") {
    in_bound_eq <- integer(length(domain$ineqs)) # A list of bool for each ineq
    for (eq_i in 1:length(domain$ineqs)) {
      eq <- domain$ineqs[[eq_i]]
      if (eq$nonnegative && any(x < 0))
        in_bound_eq[eq_i] <- FALSE
      else
        in_bound_eq[eq_i] <- tryCatch(
          if (eq$uniform) {
            s <- frac_pow(x, eq$power_numers, eq$power_denoms, eq$abs)
            (eq$larger && all(s >= eq$const)) || (!eq$larger && all(s <= eq$const))
          } else {
            s <- sum(eq$coeffs * frac_pow(x, eq$power_numers, eq$power_denoms, eq$abs))
            (eq$larger && s >= eq$const) || (!eq$larger && s <= eq$const)
          }, error = function(e) FALSE)
    }
    if (length(domain$ineqs) == 1) # If only one ineq
      return (in_bound_eq)
    # Replace ineq numbers in domain$postfix_rule by their bools
    pos <- 1
    new_string <- ""
    stack <- c()
    while (pos <= nchar(domain$postfix_rule)) {
      ch <- s_at(domain$postfix_rule, pos)
      if (grepl("[[:digit:]]",  ch)) {
        eq_num <- 0
        # If this position is a number, finish reading the ineq number
        while (pos <= nchar(domain$postfix_rule) && grepl("[[:digit:]]", ch)) {
          eq_num <- eq_num * 10 + as.integer(ch)
          pos <- pos + 1; ch <- s_at(domain$postfix_rule, pos)
        }
        # If ineq number out of bound
        if (eq_num <= 0 || eq_num > length(domain$ineqs))
          stop("In in_bound(): Equation ", eq_num, " out of bound. There are only ", length(domain$ineqs), " ineqs in total. Got postfix_rule: ", domain$postfix_rule, ".")
        stack <- c(stack, in_bound_eq[eq_num])
      } else {
        if (ch == "&" || ch == "|") {
          if (length(stack) < 2)
            stop("In in_bound(): Stack size smaller than two encountered. Got postfix_rule: ", domain$postfix_rule, ", and stack: ", stack, ".")
          if (ch == "&")
            stack <- c(stack[seq_len(length(stack)-2)], stack[length(stack)-1] && stack[length(stack)])
          else
            stack <- c(stack[seq_len(length(stack)-2)], stack[length(stack)-1] || stack[length(stack)])
        } else if (!grepl("[ \t\n\r\v\f]", ch))
          stop("In in_bound(): Invalid character '", ch, "' encountered when evaluating ineqs in postfix_rule: ", domain$postfix_rule, ".")
        pos <- pos + 1
      }
    }
    if (length(stack) != 1)
      stop("In in_bound(): Final stack should have exactly 1 element. Got: ", paste(stack, collapse=", "), ". Newstring: ", new_string, ".")
    # And evaluate by parsing the string; new_string should only contain 0, 1, &, | and spaces
    return (stack)
  }
}

#' Returns a list to be passed to C that represents the domain.
#'
#' Returns a list to be passed to C that represents the domain.
#'
#' @param domain A list returned from \code{make_domain()} that represents the domain.
#' @details Construct a list to be read by C code that represents the domain.
#' @return A list of the following elements.
#'   \item{\code{num_char_params}}{An integer, length of \code{char_params}.}
#'   \item{\code{char_params}}{A vector of string (\code{char *} or \code{char **}) parameters.}
#'   \item{\code{num_int_params}}{An integer, length of \code{int_params}.}
#'   \item{\code{int_params}}{A vector of integer (\code{int}) parameters.}
#'   \item{\code{num_double_params}}{An integer, length of \code{double_params}.}
#'   \item{\code{double_params}}{A vector of double (\code{double}) parameters.}
#' @examples
#' p <- 30
#' # The 30-dimensional real space R^30
#' domain <- make_domain("R", p=p)
#' domain_for_C(domain)
#'
#' # The non-negative orthant of the 30-dimensional real space, R+^30
#' domain <- make_domain("R+", p=p)
#' domain_for_C(domain)
#'
#' # x such that sum(x^2) > 10 && sum(x^(1/3)) > 10 with x allowed to be negative
#' domain <- make_domain("polynomial", p=p, rule="1 && 2",
#'        ineqs=list(list("expression"="sum(x^2)>10", abs=FALSE, nonnegative=FALSE),
#'                       list("expression"="sum(x^(1/3))>10", abs=FALSE, nonnegative=FALSE)))
#' domain_for_C(domain)
#'
#' # ([0, 1] v [2,3]) ^ p
#' domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))
#' domain_for_C(domain)
#'
#' # x such that {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
#' domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
#'        ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
#' domain_for_C(domain)
#'
#' # x in R_+^p such that {sum(log(x))<2 || (x1^(2/3)-1.3x2^(-3)<1 && exp(x1)+2.3*x2>2)}
#' domain <- make_domain("polynomial", p=p, rule="1 || (2 && 3)",
#'        ineqs=list(list("expression"="sum(log(x))<2", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="x1^(2/3)-1.3x2^(-3)<1", abs=FALSE, nonnegative=TRUE),
#'                       list("expression"="exp(x1)+2.3*x2^2>2", abs=FALSE, nonnegative=TRUE)))
#' domain_for_C(domain)
#'
#' # x in R_+^p such that {x in R_+^p: sum_j j * xj <= 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"=paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
#'                      abs=FALSE, nonnegative=TRUE)))
#' domain_for_C(domain)
#'
#' # The (p-1)-simplex
#' domain <- make_domain("simplex", p=p)
#' domain_for_C(domain)
#'
#' # The l-1 ball {sum(|x|) < 1}
#' domain <- make_domain("polynomial", p=p,
#'        ineqs=list(list("expression"="sum(x)<1", abs=TRUE, nonnegative=FALSE)))
#' domain_for_C(domain)
#' @export
domain_for_C <- function(domain){
  if (!"checked" %in% names(domain))
    stop("domain must be an object returned by make_domain().")
  if (domain$type == "R" || domain$type == "R+") {
    list(num_char_params=as.integer(1), char_params=toString(c(domain$type)),
         num_int_params=as.integer(0), int_params=as.integer(c()),
         num_double_params=as.integer(0), double_params=as.double(c()))
  } else if (domain$type == "uniform") {
    if (domain$left_inf) domain$lefts[1] <- 0
    if (domain$right_inf) domain$rights[length(domain$rights)] <- 0
    list(num_char_params=as.integer(1),
         char_params=toString(c("uniform")),
         num_int_params=as.integer(2),
         int_params=as.integer(c(domain$left_inf, domain$right_inf)),
         num_double_params=as.integer(2*length(domain$lefts)),
         double_params=as.double(c(domain$lefts, domain$rights)))
  } else if (domain$type == "polynomial") {
    list(num_char_params=as.integer(2),
         char_params=c("polynomial", domain$postfix_rule),
         num_int_params=as.integer(1 + 5 * length(domain$ineqs) +
                                     2 * sum(sapply(domain$ineqs,
                                                    function(eq){
                                                      length(eq$power_numers)}))),
         int_params=as.integer(
           c(length(domain$ineqs),
             Reduce("c", sapply(domain$ineqs, function(eq){
               c(eq$uniform, eq$larger, eq$abs, length(eq$power_numers)==1,
                 eq$nonnegative, eq$power_numers, eq$power_denoms)})))),
         num_double_params=as.integer(length(domain$ineqs) +
                                        sum(sapply(domain$ineqs, function(eq){1-eq$uniform}) * domain$p)),
         double_params=as.double(Reduce("c", sapply(domain$ineqs, function(eq){c(eq$coeffs, eq$const)}))))
  } else if (domain$type == "simplex") {
    list(num_char_params=as.integer(1), char_params=c("simplex"),
         num_int_params=0, int_params=integer(0), num_double_params=0, double_params=numeric(0))
  }
}

#' Parses an ineq expression into a list of elements that represents the ineq.
#'
#' Parses an ineq expression into a list of elements that represents the ineq.
#'
#' @param s A string, an ineq expression. Please refer \code{make_domain()}.
#' @param p An integer, the dimension.
#' @details Please refer \code{make_domain()} for the syntax of the expression.
#' @return A list containing the following elements:
#'  \item{uniform}{A logical, indicates whether the ineq is a uniform expression that applies to each component independently (e.g. \code{x^2>1}, \code{exp(3*|x|)<3.4}).}
#'  \item{const}{A number, the constant side of the ineq that the variable side should compare to (e.g. \code{1.3} in \code{x1^2+2*x2^3>1.3}).}
#'  \item{larger}{A logical, indicates whether the variable side of the expression should be larger or smaller than \code{const}.}
#'  \item{power_numers}{A single number or a vector of length \code{p}. The numerators of the powers in the ineq for each component (e.g. \code{c(2,3,5,0,-2)} for \cr
#'       \code{x1^2+2*x2^(3/2)+3*x3^(5/3)+4*log(x4)+5*exp(-2*x)>1}).}
#'  \item{power_denoms}{A single number or a vector of length \code{p}. The denominators of the powers in the ineq for each component (e.g. \code{c(1,2,3,0,0)} for \cr
#'       \code{x1^2+2*x2^(3/2)+3*x3^(5/3)+4*log(x4)+5*exp(-2*x)>1}).}
#'  \item{coeffs}{A vector of length \code{p} that represents the coefficients in the ineq associated with each component. Returned only if \code{uniform == FALSE}.}
#' @examples
#' p <- 30
#' parse_ineq("sum(x^2)>10", p)
#' parse_ineq("sum(x^(1/3))>10", p)
#' parse_ineq("x1>1", p)
#' parse_ineq("x2<1", p)
#' parse_ineq("exp(x)>1.3", p)
#' parse_ineq("sum(log(x)) < 2", p)
#' parse_ineq("x1^(2/3)-1.3x2^(-3)<1", p)
#' parse_ineq("exp(x1)+2.3*x2^2 > 2", p)
#' parse_ineq(paste(paste(sapply(1:p,
#'                            function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"), p)
#'
#' parse_ineq("0.5*x1^(-2/3)-x3^3 + 2log(x2)- 1.3e4exp(-25*x6)+x8-.3x5^(-3/-4) >= 2", 8)
#' parse_ineq("0.5*x1^(-2/3)-x2^(4/-6)+2e3x3^(-6/9) < 3.5e5", 3)
#' parse_ineq("x^(-2/3)<=3e3", 5)
#' parse_ineq("sum(x^(-2/3))<=3e3", 5)
#' parse_ineq("x<=3e3", 5)
#' parse_ineq("sum(x)<=3e3", 5)
#' parse_ineq("exp(-23x)<=3e3", 5)
#' parse_ineq("sum(exp(-23x))<=3e3", 5)
#' @export
parse_ineq <- function(s, p) {
  s <- base::gsub("[[:space:]]", "", s) # Remove space
  s <- base::gsub(">=", ">", base::gsub("<=", "<", s)) # Remove = from >= / <=
  const_loc <- stringr::str_locate(s, "[<>][-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$") # e.g. location of >2.34 or <-.3e6
  if (is.na(const_loc[1]))
    stop("The inequality must end with '>' or '<' followed by a constant, e.g. '... > 2.3' or '... <= -1.5e3'.")
  larger <- (s_at(s, const_loc[1]) == ">") # TRUE if > and FALSE if <
  const <- as.double(substr(s, const_loc[1]+1, const_loc[2])) # constant after > or <
  s <- substr(s, 1, const_loc[1]-1) # Remove the >/< const part

  if (grepl("(x\\D+)|(x$)", base::gsub("exp", "", s))) { # x directly followed by a non-digit (except in 'exp') must be a uniform expression like x>=2, x^2>=2, x^(-2/3)>=2, sum(exp(3x))<=2
    if (grepl("x\\d+", s))
      stop("The inequality must either be a single uniform term in x (e.g. x^2>1), or contain terms in x with an index (e.g. x1^2+x2^3>1); mixed expressions are not allowed. Got '", s, "'.")
    is_sum <- grepl("^sum\\(.+\\)$", s) # If has the form sum(a_uniform_term)
    if (is_sum)
      s <- substr(s, 5, nchar(s)-1) # Remove sum()

    tmp <- read_uniform_term(s)
    if (is.null(tmp))
      stop("Error occurred while parsing inequality as a single uniform term in x (e.g. x^2>1). Got '", s, "'.")
    if (is_sum)
      return (list("uniform"=FALSE, "larger"=larger, "power_numers"=tmp$power_numers, "power_denoms"=tmp$power_denoms,
                 "coeffs"=rep(1, p), "const"=const))
    else
      return (list("uniform"=TRUE, "larger"=larger, "power_numers"=tmp$power_numers, "power_denoms"=tmp$power_denoms,
                     "const"=const))
  }

  term_results <- list() # List of elements in each term
  appeared <- rep(FALSE, p) # Indicates whether each component has already appeared
  while (s != "") {
    res <- read_one_term(s)
    if (s == res$s)
      stop("Error occurred in parsing the inequality starting with: ", s)
    if (res$idx > p)
      stop("Component index (subscript to x) cannot exceed the domain dimension ", p, ".")
    if (appeared[res$idx])
      stop("Multiple terms encountered for x", res$idx, ". Please combine them.")
    s <- res$s
    res$s <- NULL
    term_results[[length(term_results) + 1]] <- res
    appeared[res$idx] <- TRUE
  }
  power_numers <- power_denoms <- rep(1, p)
  coeffs <- rep(0, p)
  for (res in term_results) {
    if (res$coef == 0) stop("Coefficient 0 not allowed.")
    power_numers[[res$idx]] <- res$power_numer
    power_denoms[[res$idx]] <- res$power_denom
    coeffs[[res$idx]] <- res$coef
  }
  if (length(unique(power_numers)) == 1 && length(unique(power_denoms)) == 1) {
    power_numers <- power_numers[1]
    power_denoms <- power_denoms[1]
  }
  return (list("uniform"=FALSE, "larger"=larger, "power_numers"=power_numers, "power_denoms"=power_denoms,
          "coeffs"=coeffs, "const"=const))
}

#' Attempts to parse a single term in x into power_numer and power_denom.
#'
#' Attempts to parse a single term in x into power_numer and power_denom.
#'
#' @param s A string, the variable side of an inequality expression. Please refer to \code{make_domain()}.
#' @details Returns \code{NULL} if \code{s} is not a single uniform term in \code{x} (e.g. \code{x^2} is uniform, while \code{x1^2+x2^2} is not uniform).
#' @return
#'    \item{power_numers}{The uniform numerator in the power (e.g. \code{-2} for \code{x^(-2/3)}).}
#'    \item{power_denoms}{The uniform denominator in the power (e.g. \code{3} for \code{x^(-2/3)}).}
#' @examples
#' p <- 30
#' read_uniform_term("x^2")
#' read_uniform_term("x^(1/3)")
#' read_uniform_term("exp(x)")
#' read_uniform_term("log(x)")
#' read_uniform_term("x^(-2/3)")
#' read_uniform_term("x")
#' read_uniform_term("exp(-23x)")
#' @export
read_uniform_term <- function(s) {
  # Matches uniform expressions for all components: x, log(x), exp(x), exp(23x), x^2, x^(-5/-3), no coef allowed
  # If not a uniform expression, returns NULL
  if (s == "x")
    return (list("power_numers"=1, "power_denoms"=1))
  else if (s == "log(x)")
    return (list("power_numers"=0, "power_denoms"=0))
  else if (grepl("^exp\\(", s)) { # If starts with exp(
    tmp <- read_exponential(s, has_index=FALSE)
    if (is.null(tmp) || tmp$s != "") return (NULL)
    else # If has exactly the form exp(x), exp(-x), exp(2x), exp(-123*x) with no leftover terms
      return (list("power_numers"=tmp$power_numer, "power_denoms"=0L))
  }
  else if (substr(s, 1, 2) == "x^") {
    tmp <- read_exponent(substr(s, 2, nchar(s)))
    if (is.null(tmp) || tmp$s != "") return (NULL)
    else # If s has exactly the form "x^2", "x^(2/3)", "x^(-2/3)", no leftover terms
      return (list("power_numers"=tmp$power_numer, "power_denoms"=tmp$power_denom))
  }
  return (NULL)
}

#' Parses the exponent part into power_numer and power_denom.
#'
#' Parses the exponent part into power_numer and power_denom.
#'
#' @param s A string. Must be of the form "" (empty string), "^2", "^(-5/3)" followed by other terms (starting with "+" or "-").
#' @details Parses the exponential part of the first term into power_numer and power_denom and returns the rest of the terms. Please refer to the examples. \code{s} must be of the form "", "^2", "^(-5/3)" followed by other terms, e.g. "+x2^2", "^2+x2^2", "^(-5/3)+x2^2". Assuming these come from "x1+x2^2", "x1^2+x2^2" and "x1^(-5/3)+x2^2", respectively, these will parsed into \code{power_numer=1, power_denom=1}, \code{power_numer=2, power_denom=1}, and \code{power_numer=-5, power_denom=3}, respectively.
#' @return A list containing the following elements:
#'    \item{power_numer}{An integer, the numerator of the power.}
#'    \item{power_denom}{An integer, the denominator of the power.}
#'    \item{s}{A string, the rest of the unparsed string.}
#' If parsing is unsuccessful, \code{NULL} is returned.
#' @examples
#' read_exponent("")
#' read_exponent("^(-2*4)") # Unsuccessful parsing, returns \code{NULL}.
#' read_exponent("+x2^(2/3)+x3^(-3/4)") # Comes from "x1+x2^(2/3)+x3^(-3/4)"
#' read_exponent("^2+x2^(2/3)+x3^(-3/4)") # Comes from "x1^2+x2^(2/3)+x3^(-3/4)"
#' read_exponent("^(2/3)+x2^(2/3)+x3^(-3/4)") # Comes from "x1^(2/3)+x2^(2/3)+x3^(-3/4)"
#' read_exponent("^(-2)+x2^(2/3)+x3^(-3/4)") # Comes from "x1^(-2)+x2^(2/3)+x3^(-3/4)"
#' read_exponent("^(-2/3)+x2^(2/3)+x3^(-3/4)") # Comes from "x1^(-2/3)+x2^(2/3)+x3^(-3/4)"
#' @export
read_exponent <- function(s) { # "" (exponent=1), "^2", "^(-5/3)" followed by other terms
  if (s_at(s, 1) != "^") { # No exponent
    power_numer <- power_denom <- 1
  } else if (s_at(s, 2) == "(") { # xj^(...)
    exponent_loc <- stringr::str_locate(s, "\\^\\([-]?[0-9]+(/[-]?[0-9]+)?\\)")  # ^(2), ^(-2), ^(2/3), ^(-2/3), ^(2/-3), ^(-2/-3) all allowed, although the last two are not recommended
    exp_string <- substr(s, 3, exponent_loc[2] - 1) # 2, -2, 2/3, -2/3, 2/-3, -2/-3
    power <- as.integer(strsplit(exp_string, "/")[[1]])
    if (length(power) == 1) { # 2, -2
      power_numer <- power
      power_denom <- 1
    } else if (length(power) == 2) { # 2/3, -2/3, 2/-3, -2/-3
      tmp <- makecoprime(power[1], power[2])
      power_numer <- tmp[1]
      power_denom <- tmp[2]
    } else stop(paste("Error occurred when dealing exponent in: ", s))
    s <- substr(s, exponent_loc[2] + 1, nchar(s))
  } else { # xj^INTEGER
    exponent_loc <- stringr::str_locate(s, "\\^[-]?[0-9]+") # ^2, ^-2
    power_numer <- as.integer(substr(s, 2, exponent_loc[2]))
    power_denom <- 1
    s <- substr(s, exponent_loc[2] + 1, nchar(s))
  }
  if (is.na(s))
    return (NULL)
  return (list("power_numer"=power_numer, "power_denom"=power_denom, "s"=s))
}

#' Parses the integer coefficient in an exponential term.
#'
#' Parses the integer coefficient in an exponential term.
#'
#' @param s A string that starts with one of the following forms: \code{exp(x)}, \code{exp(-x)}, \code{exp(2x)}, \code{exp(-2x)}, \code{exp(12*x)}, \code{exp(-123*x)}, followed by other terms. If \code{has_index == TRUE}, the first term should be rewritten in \code{x} with an index (e.g. \code{exp(x1)}, \code{exp(-2*x2)}).
#' @param has_index A logical, indicates whether the term is written in a component (e.g. \code{x1}, \code{x2}) as opposed to a uniform term (i.e. \code{x}).
#' @details Parses the coefficient in the first exponential term and returns the rest of the terms.
#' @return A list containing the following elements:
#'    \item{power_numer}{An integer, the integer coefficient inside the first exponential term.}
#'    \item{idx}{An integer, the index of the term matched (e.g. \code{3} for \code{exp(2*x3)}). \code{NULL} if \code{has_index == FALSE}.}
#'    \item{s}{A string, the rest of the unparsed string.}
#' If parsing is unsuccessful, \code{NULL} is returned.
#' @examples
#' # Unsuccessful parsing, not starting with exponential, returns \code{NULL}.
#' read_exponential("x", FALSE)
#' # Unsuccessful parsing, not starting with exponential, returns \code{NULL}.
#' read_exponential("x1^2+exp(2x2)", TRUE)
#' read_exponential("exp(x)", FALSE)
#' read_exponential("exp(x1)", TRUE)
#' read_exponential("exp(-x)", FALSE)
#' read_exponential("exp(-x1)+x2^2", TRUE)
#' read_exponential("exp(2x)", FALSE)
#' read_exponential("exp(2x1)+x2^(-2/3)", TRUE)
#' read_exponential("exp(-2x)", FALSE)
#' read_exponential("exp(-2x1)+exp(x3)", TRUE)
#' read_exponential("exp(12x)", FALSE)
#' read_exponential("exp(12x2)+x3^(-3)+x4^2", TRUE)
#' read_exponential("exp(-12x)", FALSE)
#' read_exponential("exp(-12x3)+x1^(2/5)+log(x2)", TRUE)
#' read_exponential("exp(123*x)", FALSE)
#' read_exponential("exp(123*x1)+x2^4", TRUE)
#' read_exponential("exp(-123*x)", FALSE)
#' read_exponential("exp(-123*x4)+exp(2*x3)", TRUE)
#' @export
read_exponential <- function(s, has_index) {
  # If has_index == FALSE, matches exp(x), exp(-x), exp(2x), exp(-2x), exp(12x), exp(-12x), exp(123*x), exp(-123*x)
  # If has_index == TRUE, matches above with x followed by an integer component index (e.g. x1, x2)
  if (has_index) {
    exp_loc <- stringr::str_locate(s, "^exp\\((([-]?[0-9]+\\*)|([-]?[0-9]*))?x\\d+\\)")
  } else {
    exp_loc <- stringr::str_locate(s, "^exp\\((([-]?[0-9]+\\*)|([-]?[0-9]*))?x\\)")
  }
  if (is.na(exp_loc[1]))
    return (NULL)
  inside_parentheses <- substr(s, 5, exp_loc[2]-1)
  power_numer <- base::gsub("\\*", "", stringr::str_extract(inside_parentheses, ".*(?=x)")) # Coefficient on x, e.g. "", "-", "2", "-2*"
  if (power_numer == "")
    power_numer <- 1L
  else if (power_numer == "-")
    power_numer <- -1L
  else
    power_numer <- as.integer(power_numer)
  if (has_index)
    idx <- as.integer(stringr::str_extract(inside_parentheses, "(?<=x).*"))
  else
    idx <- NULL
  s <- substr(s, exp_loc[2]+1, nchar(s))
  return (list("power_numer"=power_numer, "idx"=idx, "s"=s))
}

#' Parses the first term of a non-uniform expression.
#'
#' Parses the first term of a non-uniform expression.
#'
#' @param s A string, the variable side of a non-uniform inequality expression (i.e. terms must be rewritten in e.g. \code{x1}, \code{x2} as opposed to \code{x}).
#' @details Parses the first term in a non-uniform expression and returns the rest of the terms.
#' @return A list containing the following elements:
#'    \item{idx}{An integer, the index of the first term (e.g. \code{3} for \code{1.3*x3^(-2/5))}).}
#'    \item{power_numer}{An integer, the power_numer of the first term.}
#'    \item{power_denom}{An integer, the power_denom of the first term.}
#'    \item{coef}{A number, the coefficient on the first term (e.g. \code{1.3} for \code{1.3*x3^(-2/5)}).}
#'    \item{s}{A string, the rest of the unparsed string.}
#' @examples
#' read_one_term("0.5*x1+x2^2")
#' read_one_term("2e3x1^(2/3)-1.3x2^(-3)")
#' read_one_term("2exp(3x1)+2.3*x2^2")
#' read_one_term(paste(sapply(1:10, function(j){paste(j, "x", j, "^", (11-j), sep="")}), collapse="+"))
#' read_one_term("0.5*x1^(-2/3)-x3^3 + 2log(x2)- 1.3e4exp(-25*x6)+x8-.3x5^(-3/-4)")
#' read_one_term("-1e-4x1^(-2/3)-x2^(4/-6)+2e3x3^(-6/9) < 3.5e5")
#' @export
read_one_term <- function(s) { # Assumes no whitespace, read one term from an inequality
  coef_loc <- stringr::str_locate(s, "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?") # Match an integer/double coefficient, potentially with scientific notation
  if (is.na(coef_loc[1])) { # No explicit coef, 1 or -1 depending on sign
    if (s_at(s, 1) == "-") {
      coef <- -1
      s <- substr(s, 2, nchar(s)) # Remove the - sign
    } else if (s_at(s, 1) == "+") {
      coef <- 1
      s <- substr(s, 2, nchar(s))
    } else {
      coef <- 1
    }
  } else {
    coef <- as.double(substr(s, coef_loc[1], coef_loc[2]))
    s <- substr(s, coef_loc[2]+1, nchar(s))
  }
  if (s_at(s, 1) == "*") # Ignore the multiple sign
    s <- substr(s, 2, nchar(s))
  if (grepl("^log\\(x\\d+\\)", s)) {
    power_numer <- power_denom <- 0L
    logx_loc <- stringr::str_locate(s, "log\\(x\\d+\\)")
    idx <- as.integer(substr(s, 6, logx_loc[2]-1)) # log(xIDX)
    s <- substr(s, logx_loc[2]+1, nchar(s))
  } else if (grepl("^exp\\(", s)) {
    tmp <- read_exponential(s, has_index=TRUE)
    if (is.null(tmp))
      stop("Terms with exp must have the form exp(x1), exp(-x2), exp(2x3), exp(-12x4), exp(123*x5), exp(-1234*x5), or x with no index if a uniform expression. Got '", s, "'.")
    power_numer <- tmp$power_numer
    power_denom <- 0L
    idx <- tmp$idx
    s <- tmp$s
  } else { # Regular power
    if (s_at(s, 1) != "x" || (!grepl("\\d", s_at(s, 2))))
      stop("Each term on the left hand side of the inequality string must be '[optional - sign][optional double coefficient with/without *]' followed by 'log(xj)' OR 'exp([optional coef*]xj)' OR 'xj[^optional exponent]', where j is the component index (starting with 1) and the optional exponent following the ^ sign must be a single integer, or (integer/integer), or (-integer/integer).",
           " Or a uniform expression with no index and no coefficient can be used, e.g. the entire expression can be 'x^2 >= 5', 'x^(-3/5) <= 1.3', 'log(x) > 3.9', 'exp(-3x) > 4.9', exp(x) < 2.6', 'x < 0', but not '3.5*x > 1'. Got '", s, "'.")
    xname_loc <- stringr::str_locate(s, "x\\d+")
    idx <- as.integer(substr(s, xname_loc[1]+1, xname_loc[2]))
    s <- substr(s, xname_loc[2]+1, nchar(s))
    if (is.null(s))
      stop("Error occurred in parsing the inequality.")
    tmp <- read_exponent(s)
    power_numer <- tmp$power_numer
    power_denom <- tmp$power_denom
    s <- tmp$s
  }
  return (list("idx"=idx, "power_numer"=power_numer, "power_denom"=power_denom,
               "coef"=coef, "s"=s))
}

#' Finds the intersection between two unions of intervals.
#'
#' Finds the intersection between two unions of intervals.
#'
#' @param A A list of vectors of size 2, each representing an interval. It is required that \code{A[[i]][1] <= A[[i]][2] <= A[[j]][1]} for any \code{i < j}.
#' @param B A list of vectors of size 2, each representing an interval. It is required that \code{A[[i]][1] <= A[[i]][2] <= A[[j]][1]} for any \code{i < j}.
#' @details Finds the intersection between the union of all intervals in \code{A} and the union of all intervals in \code{B}.
#' @return A list of vectors of size 2, whose union represents the intersection between \code{A} and \code{B}.
#' @examples
#' interval_intersection(list(c(1.2,1.5), c(2.3,2.7)),
#'        list(c(0.6,1.4), c(2.5,3.6), c(6.3,6.9)))
#' interval_intersection(list(c(-0.3,0.55), c(2.35,2.8)),
#'        list(c(0.54,0.62), c(2.5,2.9)))
#' interval_intersection(list(c(0,1)), list(c(1,2)))
#' interval_intersection(list(c(0,1+1e-8)), list(c(1,2)))
#' interval_intersection(list(c(0,1), c(2,3)),
#'        list(c(1,2)))
#' interval_intersection(list(c(0,1+1e-8), c(2-1e-8,3)),
#'        list(c(1,2)))
#' interval_intersection(list(c(0,1)), list())
#' @export
interval_intersection <- function(A, B) {
  ai <- bi <- 1; res <- list()
  while (ai <= length(A) && bi <= length(B)) {
    if (A[[ai]][2] > B[[bi]][1] && B[[bi]][2] > A[[ai]][1])
      res <- c(res, list(c(max(A[[ai]][1], B[[bi]][1]), min(A[[ai]][2], B[[bi]][2]))))
    if (A[[ai]][2] > B[[bi]][2]) bi <- bi + 1
    else if (A[[ai]][2] < B[[bi]][2]) ai <- ai + 1
    else {ai <- ai + 1; bi <- bi + 1}
  }
  return (res)
}

#' Finds the union betweeen two unions of intervals.
#'
#' Finds the union betweeen two unions of intervals.
#'
#' @param A A list of vectors of size 2, each representing an interval. It is required that \code{A[[i]][1] <= A[[i]][2] <= A[[j]][1]} for any \code{i < j}.
#' @param B A list of vectors of size 2, each representing an interval. It is required that \code{A[[i]][1] <= A[[i]][2] <= A[[j]][1]} for any \code{i < j}.
#' @details Finds the union between the union of all intervals in \code{A} and the union of all intervals in \code{B}.
#' @return A list of vectors of size 2, whose union represents the union between \code{A} and \code{B}.
#' @examples
#' interval_union(list(c(1.2,1.5), c(2.3,2.7)),
#'        list(c(0.6,1.4), c(2.5,3.6), c(6.3,6.9)))
#' interval_union(list(c(-0.3,0.55), c(2.35,2.8)),
#'        list(c(0.54,0.62), c(2.5,2.9)))
#' interval_union(list(c(0,1)), list(c(1,2)))
#' interval_union(list(c(0,1-1e-8)), list(c(1,2)))
#' interval_union(list(c(0,1), c(2,3)),
#'        list(c(1,2)))
#' interval_union(list(c(0,1-1e-8), c(2+1e-8,3)),
#'        list(c(1,2)))
#' interval_union(list(c(0,1)), list())
#' @export
interval_union <- function(A, B) {
  if (length(A) == 0) return (B)
  if (length(B) == 0) return (A)
  lefts <- sort(c(sapply(A, function(x){x[1]}), sapply(B, function(x){x[1]})))
  rights <- sort(c(sapply(A, function(x){x[2]}), sapply(B, function(x){x[2]})))
  res <- list(c(lefts[1], rights[1]))
  for (i in 2:length(lefts)) { # length(lefts) >= 2
    if (lefts[i] <= res[[length(res)]][2]) res[[length(res)]][2] <- rights[i]
    else res <- c(res, list(c(lefts[i], rights[i])))
  }
  return (res)
}
