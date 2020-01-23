library(genscore)
if (!requireNamespace("tmvtnorm", quietly = TRUE))
  stop("Please install package \"tmvtnorm\".")
library(tmvtnorm)

readline("For a full tutorial on the package, please refer to the vignette instead.\n(Press enter to continue)")

## Domain definition
p <- 30
# The 30-dimensional real space R^30
domain <- make_domain("R", p=p)

# The non-negative orthant of the 30-dimensional real space, R+^30
domain <- make_domain("R+", p=p)

# x such that sum(x^2) > 10 && sum(x^(1/3)) > 10 with x allowed to be negative
domain <- make_domain("polynomial", p=p, rule="1 && 2",
                      ineqs=list(list("expression"="sum(x^2)>10", abs=FALSE, nonnegative=FALSE),
                                 list("expression"="sum(x^(1/3))>10", abs=FALSE, nonnegative=FALSE)))
# Alternatively, equivalent domain
domain2 <- make_domain("polynomial", p=p, rule="1 && 2",
                       ineqs=list(list(uniform=FALSE, power_numers=2, power_denoms=1, const=10, 
                                       coeffs=1, larger=1, abs=FALSE, nonnegative=FALSE),
                                  list(uniform=FALSE, power_numers=1, power_denoms=3, const=10, 
                                       coeffs=1, larger=1, abs=FALSE, nonnegative=FALSE)))

# ([0, 1] v [2,3]) ^ p
domain <- make_domain("uniform", p=p, lefts=c(0,2), rights=c(1,3))

# x such that {x1 > 1 && log(1.3) < x2 < 1 && x3 > log(1.3) && ... && xp > log(1.3)}
domain <- make_domain("polynomial", p=p, rule="1 && 2 && 3",
                      ineqs=list(list("expression"="x1>1", abs=FALSE, nonnegative=TRUE),
                                 list("expression"="x2<1", abs=FALSE, nonnegative=TRUE),
                                 list("expression"="exp(x)>1.3", abs=FALSE, nonnegative=FALSE)))
# Alternatively, equivalent domain
domain2 <- make_domain("polynomial", p=p, rule="1 && 2",
                       ineqs=list(list(uniform=FALSE, power_numers=1, power_denoms=1, const=1, 
                                       coeffs=c(1,rep(0,p-1)), larger=1, abs=FALSE, nonnegative=TRUE),
                                  list(uniform=FALSE, power_numers=1, power_denoms=1, const=1, 
                                       coeffs=c(0,1,rep(0,p-2)), larger=0, abs=FALSE, nonnegative=TRUE),
                                  list(uniform=TRUE, power_numers=1, power_denoms=0, const=1.3,
                                       larger=1, abs=FALSE, nonnegative=FALSE)))

# x in R_+^p such that {sum(log(x))<2 || (x1^(2/3)-1.3x2^(-3)<1 && exp(x1)+2.3*x2>2)}
domain <- make_domain("polynomial", p=p, rule="1 || (2 && 3)",
                      ineqs=list(list("expression"="sum(log(x))<2", abs=FALSE, nonnegative=TRUE),
                                 list("expression"="x1^(2/3)-1.3x2^(-3)<1", abs=FALSE, nonnegative=TRUE),
                                 list("expression"="exp(x1)+2.3*x2^2>2", abs=FALSE, nonnegative=TRUE)))
# Alternatively, equivalent domain
domain2 <- make_domain("polynomial", p=p, rule="1 && 2",
                       ineqs=list(list(uniform=FALSE, power_numers=0, power_denoms=0, const=2, 
                                       coeffs=1, larger=0, abs=FALSE, nonnegative=TRUE),
                                  list(uniform=FALSE, power_numers=c(2,-3,rep(1,p-2)), power_denoms=c(3,rep(1,p-1)), 
                                       const=1, coeffs=c(1.0,-1.3,rep(0,p-2)), larger=0, abs=FALSE, nonnegative=TRUE),
                                  list(uniform=FALSE, power_numers=c(1,2,rep(1,p-2)), power_denoms=c(0,rep(1,p-1)), 
                                       const=2, coeffs=c(1,2.3,rep(0,p-2)), larger=1, abs=FALSE, nonnegative=TRUE)))

# x in R_+^p such that {x in R_+^p: sum_j j * xj <= 1}
domain <- make_domain("polynomial", p=p, 
                      ineqs=list(list("expression"=paste(paste(sapply(1:p,
                        function(j){paste(j, "x", j, sep="")}), collapse="+"), "<1"),
                        abs=FALSE, nonnegative=TRUE)))
# Alternatively, equivalent domain
domain2 <- make_domain("polynomial", p=p,
                       ineqs=list(list(uniform=FALSE, power_numers=1, power_denoms=1, const=1, 
                                       coeffs=1:p, larger=0, abs=FALSE, nonnegative=TRUE)))

# The (p-1)-dimensional simplex
domain <- make_domain("simplex", p=p)

# The l-1 ball {sum(|x|) < 1}
domain <- make_domain("polynomial", p=p, 
                      ineqs=list(list("expression"="sum(x)<1", abs=TRUE, nonnegative=FALSE)))

readline("Press enter to continue to the next example.")




## Gaussian Graphical Models on the non-negative orthant R_+^p
n <- 200
p <- 50
tol <- 1e-8
domain <- make_domain("R+", p=p)
# Mean parameter
mu <- rnorm(p) * 0.5
# Inverse covariance matrix
K <- cov_cons(mode="sub", p=p, seed=1, spars=0.2, eig=0.1, subgraphs=5)
# True edge set
true_edges <- which(abs(K) > tol & diag(p) == 0)
# Generate a sample from truncated GGMs
set.seed(1)
x <- rtmvnorm(n, mean = mu, sigma = solve(K), lower = rep(0, p), upper = rep(Inf, p), algorithm = "gibbs", burn.in.samples = 100, thinning = 10)

### Estimate the graph using estimate()
# max number of iterations
maxit <- 1000
# number of lambdas
nlambda <- 100
# diagonal multiplier
dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n)))) # diagonal multiplier
# h function
h_mode <- "min_pow"
h_param1 <- 1
h_param2 <- 3
# estimate assuming mean unknown, using the symmetric profiled estimator
est <- estimate(x, setting="gaussian", domain=domain, centered=FALSE, symmetric="symmetric", 
                scale="norm", lambda_length=nlambda, lambda_ratio=Inf, mode=h_mode, 
                param1=h_param1, param2=h_param2, verbose=TRUE, tol=tol, maxit=maxit, 
                BIC_refit=TRUE, diagonal_multiplier=dm,  return_raw=TRUE)

# The estimated graph for the largest automatically chosen lambda
print(est$edgess[[1]]) # empty graph
# The estimated graph for the smallest automatically chosen lambda
print(length(est$edgess[[nlambda]]) == p*(p-1)) # complete graph
# Estimated graph for the 50th lambda
print(length(est$edgess[[nlambda/2]]))
readline("Press enter to continue to the next example.")

### BIC scores
print(est$BICs[1:10,])
# Refitted BIC scores
print(est$BIC_refits[1:10, ])
readline("Press enter to continue to the next example.")

### lambda values
print(est$lambda1s) # lambda_Ks used
print(est$lambda2s) # lambda_etas, all 0 since with lambda_ratio=Inf we are in the profiled setting
readline("Press enter to continue to the next example.")

### eta estimates
print(est$etas[1:5,1:5]) # First 5 components for the first 5 lambdas
print(est$raw_estimates[[nlambda/2]][1:5,1:5]) # First few entries of the estimated K for the 50th lambda
readline("Press enter to continue to the next example.")

### ROC curve
ROC <- t(sapply(est$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
# The area under the curve
AUC(ROC)
# Plot refitted BIC and ROC
par(mfrow=c(1,2), mar=c(4,4,4,4))
colors_ <- rainbow(nlambda)
plot_BIC_ind <- !is.infinite(est$BIC_refits[,3])
plot(est$lambda1s[plot_BIC_ind], est$BIC_refits[plot_BIC_ind,3], col=colors_[plot_BIC_ind], pch=16, cex=1.5, xlab="lambda", ylab="eBIC", main="Refitted BIC (when exists)")
plot(c(),c(), ylim = c(0,1), xlim = c(0,1), cex.lab=1, xlab = "False Positives", ylab = "True Positives", main = "ROC curve")
points(ROC[,2], ROC[,1], type="l")
points(ROC[,2], ROC[,1], pch=16, cex=1.5, col=colors_)
points(c(0,1), c(0,1), type="l", lty=2)
readline("Press enter to continue to the next example.")

### Estimation using estimate() with get_elts()
# Get the h and h' functions
h_hp <- get_h_hp(mode=h_mode, para=h_param1, para2=h_param2) # h(x)=min(x,3), h'(x)=(x<=3)
# Sufficient statistics, named "elts"
elts_P <- get_elts(h_hp, x, setting="gaussian", domain=domain, centered=FALSE, 
                   profiled_if_noncenter=Inf, scale="norm", diagonal_multiplier=dm, 
                   use_C=TRUE, tol=tol) # P stands for Profiled
### Estimate using estimate() with elts_P
est2 <- estimate(x, setting="gaussian", elts=elts_P, domain=domain, centered=FALSE, 
                 symmetric="symmetric", scale="norm", lambda_length=nlambda, 
                 lambda_ratio=Inf, h_hp=h_hp, verbose=TRUE, tol=tol, 
                 maxit=maxit, BIC_refit=TRUE, diagonal_multiplier=dm,  return_raw=TRUE)
# Compare to the earlier results; should be all close to 0
print(compare_two_results(est, est2))
readline("Press enter to continue to the next example.")


### Estimation for one lambda value
i <- 10
# Re-estimate at the 10th lambda
lambda1 <- est$lambda1s[i]
sub_res <- get_results(elts_P, symmetric="symmetric", lambda1=lambda1, lambda2=0, 
                       tol=tol, maxit=maxit, previous_res=NULL, is_refit=FALSE)
# Check if the result matches the previous one
print(sum(abs(est$raw_estimates[[i]] - sub_res$K)))
print(sum(abs(est$etas[i,] - sub_res$eta)))
print(all(sub_res$edges == est$edgess[[i]]))
readline("Press enter to continue to the next example.")
## Equivalently, use estimate() with one lambda value
est_sub <- estimate(x, setting="gaussian", centered=FALSE, domain=domain, 
                    symmetric="symmetric", scale="norm", lambda1s = lambda1, 
                    lambda_ratio=Inf, mode=h_mode, param1=h_param1, param2=h_param2, 
                    verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=TRUE, diagonal_multiplier=dm, 
                    return_raw=TRUE)
# Check if the result matches the previous one
print(sum(abs(est_sub$raw_estimates[[1]] - sub_res$K)))
print(sum(abs(est_sub$etas[1,] - sub_res$eta)))
print(all(est_sub$edgess[[1]] == est$edgess[[1]]))
readline("Press enter to continue to the next example.")


## The same function can also be used for fitting an unpenalized model restricted
## to a specific edge set when e.g. calculating the refitted eBIC by
## restricting to the estimated edge set. One can thus manually calculate
## the refitted eBIC in this manner, but we encourage using the function
## wrapper `estimate()` directly, which produces this information automatically.
est_refit <- get_results(elts_P, symmetric="symmetric", lambda1=lambda1, lambda2=0, 
                         tol=tol, maxit=maxit, previous_res=sub_res, is_refit=TRUE)
# Manually calculate the BIC score (eBIC with gamma = 0)
print(2*n*est_refit$crit + (length(est_refit$edges)/2+length(est_refit$eta_support))*log(n))
# Should match the number returned previously by estimate()
print(est$BIC_refit[i,1])
readline("Press enter to continue to the next example.")

### Aggregating multiple ROC curves
## The avgrocs() function aggregates multiple ROC curves using vertical
## averaging by adopting algorithm 3 from Fawcett (2006).
ROCs <- NULL
par(mfrow=c(2,2), mar=c(3,3,3,3))
for (i in 1:3){
  set.seed(i)
  x <- rtmvnorm(n, mean = mu, sigma = solve(K), lower = rep(0, p), 
                upper = rep(Inf, p), algorithm = "gibbs", burn.in.samples = 100, 
                thinning = 10)
  est <- estimate(x, setting="gaussian", centered=FALSE, domain=domain,
                  symmetric="symmetric", scale="norm", lambda_length=nlambda, 
                  lambda_ratio=Inf, mode=h_mode, param1=h_param1, param2=h_param2, 
                  verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=FALSE, 
                  diagonal_multiplier=dm, return_raw=FALSE)
  ROC <- t(sapply(est$edgess, function(edges){tp_fp(edges, true_edges, p)}))
  ROCs[[i]] <- ROC
  plot(c(), c(),  ylim=c(0,1), xlim=c(0,1), cex.lab=1, main=paste("ROC, trial ",i,", AUC ",round(AUC(ROC),4),sep=""), xlab="False Positives", ylab="True Positives")
  points(ROC[,2], ROC[,1], type="l")
  points(c(0,1), c(0,1), type = "l", lty = 2)
}
average_ROC <- avgrocs(ROCs, length(true_edges), p)
plot(c(), c(),  ylim=c(0,1), xlim=c(0,1), cex.lab=1, main=paste("Average ROC, AUC",round(AUC(average_ROC),4)), xlab="False Positives", ylab="True Positives")
points(average_ROC[,2], average_ROC[,1], type="l")
points(c(0,1), c(0,1), type = "l", lty = 2)
readline("Press enter to continue to the next example.")



## Exponential Square-Root Graphical Models on the non-negative orthant R_+^p
# Consider the exponential square-root graphical models from Inouye et al (2016),
# also defined as a special case of the a-b model in (1) in Yu et al (2019)
# with a=b=0.5.
x <- gen(n, a_numer=1, a_denom=2, b_numer=1, b_denom=2, abs=FALSE, eta=mu, K=K, 
         domain=domain, finite_infinity=100, seed=2, burn_in=1000, thinning=1000, verbose=TRUE, 
         remove_outofbound=TRUE) 

## Estimation
# As suggested in Yu et al (2019), we use h(x)=min(x,3)^1.5.
# We also try with a finite lambda_ratio as it is known to improve
# the result compared to the profiled estimator, although the best ratio
# is subject to careful tuning.
h_mode <- "min_pow"
h_param1 <- 1.5
h_param2 <- 3
est_exp <- estimate(x, setting="exp", centered=FALSE, domain=domain, symmetric="symmetric", 
                    scale="norm", lambda_length=nlambda, lambda_ratio=2, mode=h_mode, 
                    param1=h_param1, param2=h_param2, verbose=TRUE, tol=tol, maxit=maxit,
                    BIC_refit=TRUE, diagonal_multiplier=dm,  return_raw=FALSE)
ROC_exp <- t(sapply(est_exp$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
AUC(ROC_exp) # The area under the curve
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot_BIC_ind <- !is.infinite(est_exp$BIC_refits[,3])
plot(est_exp$lambda1s[plot_BIC_ind], est_exp$BIC_refits[plot_BIC_ind,3], col=colors_[plot_BIC_ind], pch=16, cex=1.5, xlab="lambda", ylab="eBIC", main="Refitted BIC (when exists)")
plot(c(),c(), ylim = c(0,1), xlim = c(0,1), cex.lab=1, xlab = "False Positives", ylab = "True Positives", main = "ROC curve")
points(ROC_exp[,2], ROC_exp[,1], type="l")
points(ROC_exp[,2], ROC_exp[,1], pch=16, cex=1.5, col=colors_)
points(c(0,1), c(0,1), type="l", lty=2)
readline("Press enter to continue to the next example.")


## Gamma Graphical Models on the non-negative orthant R_+^p
# Gamma models (a=0.5, b=0), again with h(x)=\min(x,3)^1.5 and
# `lambda_ratio=2`. Recall that the gamma graphical models require
# all entries in the linear parameter to be strictly larger than -1.
mu[mu <= -1] <- abs(mu[mu <= -1])
x <- gen(n, a_numer=1, a_denom=2, b_numer=0, b_denom=0, abs=FALSE, eta=mu, K=K, 
         domain=domain, finite_infinity=100, seed=3, burn_in=1000, thinning=1000, verbose=TRUE, 
         remove_outofbound=TRUE) 
h_mode <- "min_pow"
h_param1 <- 1.5
h_param2 <- 3
est_gamma <- estimate(x, setting="gamma", centered=FALSE, domain=domain, 
                      symmetric="symmetric", scale="norm", lambda_length=nlambda, 
                      lambda_ratio=2, mode=h_mode, param1=h_param1, param2=h_param2, 
                      verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=TRUE,
                      diagonal_multiplier=dm,  return_raw=FALSE)
ROC_gamma <- t(sapply(est_gamma$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
AUC(ROC_gamma) # The area under the curve
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot_BIC_ind <- !is.infinite(est_gamma$BIC_refits[,3])
plot(est_gamma$lambda1s[plot_BIC_ind], est_gamma$BIC_refits[plot_BIC_ind,3], col=colors_[plot_BIC_ind], pch=16, cex=1.5, xlab="lambda", ylab="eBIC", main="Refitted BIC (when exists)")
plot(c(),c(), ylim = c(0,1), xlim = c(0,1), cex.lab=1, xlab = "False Positives", ylab = "True Positives", main = "ROC curve")
points(ROC_gamma[,2], ROC_gamma[,1], type="l")
points(ROC_gamma[,2], ROC_gamma[,1], pch=16, cex=1.5, col=colors_)
points(c(0,1), c(0,1), type="l", lty=2)
readline("Press enter to continue to the next example.")


## General a-b Graphical Models on the non-negative orthant R_+^p
# Now choose a=1.5, b=0.5, and estimate using h(x)=\min(x,3)^0.5
# and `lambda_ratio=2`.
x <- gen(n, a_numer=3, a_denom=2, b_numer=3, b_denom=2, abs=FALSE, eta=mu, K=K, 
         domain=domain, finite_infinity=100, seed=4, burn_in=1000, thinning=1000, verbose=TRUE, 
         remove_outofbound=TRUE) 
h_mode <- "min_pow"
h_param1 <- 0.5
h_param2 <- 3
est_ab <- estimate(x, setting="ab_1.5_0.5", centered=FALSE, domain=domain,
                   symmetric="symmetric", scale="norm", lambda_length=nlambda, 
                   lambda_ratio=2, mode=h_mode, param1=h_param1, param2=h_param2, 
                   verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=TRUE, 
                   diagonal_multiplier=dm,  return_raw=FALSE)
ROC_ab <- t(sapply(est_ab$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
AUC(ROC_ab) # The area under the curve
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot_BIC_ind <- !is.infinite(est_ab$BIC_refits[,3])
plot(est_ab$lambda1s[plot_BIC_ind], est_ab$BIC_refits[plot_BIC_ind,3], col=colors_[plot_BIC_ind], pch=16, cex=1.5, xlab="lambda", ylab="eBIC", main="Refitted BIC (when exists)")
plot(c(),c(), ylim = c(0,1), xlim = c(0,1), cex.lab=1, xlab = "False Positives", ylab = "True Positives", main = "ROC curve")
points(ROC_ab[,2], ROC_ab[,1], type="l")
points(ROC_ab[,2], ROC_ab[,1], pch=16, cex=1.5, col=colors_)
points(c(0,1), c(0,1), type="l", lty=2)
readline("Press enter to continue to the next example.")


## (Untruncated) Gaussian Graphical Models on the entire R^p
# Now estimating the inverse covariance matrix from GGMs on the entire
# R^p. However, for Gaussian graphical models we use the penalized
# version of the original score matching from Hyv\"arinen (2005);
# see Lin et al (2016). Thus, the h and hp functions are ignored
# (where we fix h(x)=1). Here we estimate with `lambda_ratio=2`.
set.seed(5)
x <- rmvnorm(n, mean = mu, sigma = solve(K))
domain <- make_domain(type="R", p=p)
est_gauss <- estimate(x, setting="gaussian", centered=FALSE, domain=domain,
                      symmetric="symmetric", scale="norm", lambda_length=nlambda, 
                      lambda_ratio=2, verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=TRUE, 
                      diagonal_multiplier=dm,  return_raw=TRUE)
# Alternatively, get the summary statistics first (without specifying h and hp)
elts_gauss <- get_elts(NULL, x, setting="gaussian", centered=FALSE, domain=domain,
                       profiled_if_noncenter=FALSE, scale="norm", diagonal_multiplier=dm, 
                       use_C=TRUE, tol=tol)
# Then pass it to estimate()
est_gauss2 <- estimate(x, setting="gaussian", elts=elts_gauss, centered=FALSE, 
                       domain=domain, symmetric="symmetric", scale="norm", 
                       lambda_length=nlambda, lambda_ratio=2, verbose=TRUE, tol=tol, 
                       maxit=maxit, BIC_refit=TRUE, diagonal_multiplier=dm,  
                       return_raw=TRUE)
# Compare with the previous results; should be all close to 0
print(compare_two_results(est_gauss, est_gauss2))

# ROC curves
ROC_gauss <- t(sapply(est_gauss$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
AUC(ROC_gauss) # The area under the curve
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot_BIC_ind <- !is.infinite(est_gauss$BIC_refits[,3])
plot(est_gauss$lambda1s[plot_BIC_ind], est_gauss$BIC_refits[plot_BIC_ind,3], col=colors_[plot_BIC_ind], pch=16, cex=1.5, xlab="lambda", ylab="eBIC", main="Refitted BIC (when exists)")
plot(c(),c(), ylim = c(0,1), xlim = c(0,1), cex.lab=1, xlab = "False Positives", ylab = "True Positives", main = "ROC curve")
points(ROC_gauss[,2], ROC_gauss[,1], type="l")
points(ROC_gauss[,2], ROC_gauss[,1], pch=16, cex=1.5, col=colors_)
points(c(0,1), c(0,1), type="l", lty=2)
readline("Press enter to continue to the next example.")


## Aitchison $A^d$ Models on the Simplex
# Consider the $A^d$ models from Aitchison (1985), which has density proportional to
# exp(-1/2*t(log(x)) %*% K %*% log(x) - t(eta) %*% log(x)) on the (p-1)-dimensional simplex.
# In the original model, one assumes that all(rowSums(K) == 0) and all(colSums(K) == 0), 
# and all(eta > -1). This is called the "log_log_sum0" model in the package.

eta <- rep(0.5, p)
K <- -cov_cons("band", p, seed=1, spars=3)
diag(K) <- diag(K) - rowSums(K) # So that K has row and column sums 0
eigen(K)$val # Verify that K has one 0 and (p-1) positive eigenvalues 
true_edges <- which(abs(K) > tol & diag(p) == 0)
domain <- make_domain(type="simplex", p=p) # Simplex domain
x <- gen(n, a_numer=0, a_denom=0, b_numer=0, b_denom=0, abs=FALSE, eta=eta, K=K, 
         domain=domain, finite_infinity=100, seed=6, burn_in=1000, thinning=1000, 
         verbose=TRUE, remove_outofbound=TRUE)
h_mode <- "pow" # Simplex domains are bounded by nature, so no truncation needed
h_param1 <- 2 
# scale = "" since simplex not invariant to scaling
est_log_log_sum0 <- estimate(x, setting="log_log_sum0", domain=domain, centered=FALSE, 
                             symmetric="symmetric", scale="", lambda_length=nlambda, 
                             lambda_ratio=Inf, mode=h_mode, param1=h_param1, param2=NULL, 
                             verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=TRUE, 
                             diagonal_multiplier=NULL,  return_raw=FALSE)

ROC_log_log_sum0 <- t(sapply(est_log_log_sum0$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
AUC(ROC_log_log_sum0) # The area under the curve
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot_BIC_ind <- !is.infinite(est_log_log_sum0$BIC_refits[,3])
plot(est_log_log_sum0$lambda1s[plot_BIC_ind], est_log_log_sum0$BIC_refits[plot_BIC_ind,3], col=colors_[plot_BIC_ind], pch=16, cex=1.5, xlab="lambda", ylab="eBIC", main="Refitted BIC (when exists)")
plot(c(),c(), ylim = c(0,1), xlim = c(0,1), cex.lab=1, xlab = "False Positives", ylab = "True Positives", main = "ROC curve")
points(ROC_log_log_sum0[,2], ROC_log_log_sum0[,1], type="l")
points(ROC_log_log_sum0[,2], ROC_log_log_sum0[,1], pch=16, cex=1.5, col=colors_)
points(c(0,1), c(0,1), type="l", lty=2)
readline("Press enter to continue to the next example.")

# If one does not assume the row and column sums of K are all zero, the setting is 
# "log_log" instead of "log_log_sum0".
est_log_log <- estimate(x, setting="log_log", domain=domain, centered=FALSE, 
                        symmetric="symmetric", scale="", lambda_length=40, 
                        lambda_ratio=Inf, mode=h_mode, param1=h_param1, param2=NULL, 
                        verbose=TRUE, tol=tol, maxit=maxit, BIC_refit=TRUE, 
                        diagonal_multiplier=NULL,  return_raw=FALSE)
ROC_log_log <- t(sapply(est_log_log$edgess, tp_fp, true_edges=true_edges, p=p)) # ROC curve
AUC(ROC_log_log)
readline("Press enter to continue to the next example.")


# Univariate Truncated Normal Distributions
# Consider univariate truncated normal distributions with probability
# density function proportional to exp(-(x-mu)^2/(2*sigma^2)) on (0,Inf)$.
# First fix an $h$ function, and estimate the parameters from a random sample.
mode <- "min_log_pow"
param1 <- 1
param2 <- 2
n <- 1000
mu <- 0
sigma <- 1
set.seed(6)
x <- tmvtnorm::rtmvnorm(n, mean=mu, sigma=sigma^2, lower=c(0), upper=c(Inf))
## Can provide one of the parameters if the true value is known, or estimate both.
# Assuming sigma is known, estimated muhat should be close to the truth 0
print(mu_sigmasqhat(x, mode, param1, param2, sigma=sigma^2))
# Assuming mu is known, estimated sigmasqhat should be close to the truth 1
print(mu_sigmasqhat(x, mode, param1, param2, mu=mu))
# Assuming both are unknown, estimates should be close to the truth c(0,1)
print(mu_sigmasqhat(x, mode, param1, param2))
readline("Press enter to continue to the next example.")


## Variance Estimation and Confidence Intervals
# varhat() calculates the asymptotic variance of the estimator of
# one parameter assuming knowing the other.
# The asymptotic variance for mu assuming sigmasq known
print(varhat(mu, sigma^2, mode, param1, param2, est_mu=TRUE))
# The asymptotic variance for sigmasq assuming mu known
print(varhat(mu, sigma^2, mode, param1, param2, est_mu=FALSE))
# Plug in the estimated values and construct a confidence interval using
# the asymptotic variance evaluated at the estimates.
alpha <- 0.05
muhat <- mu_sigmasqhat(x, mode, param1, param2, sigma=sigma^2)[1]
# Assuming sigmasq is known, plug in muhat
sdmu <- sqrt(varhat(muhat, sigma^2, mode, param1, param2, est_mu=TRUE)/n)
# The 95% confidence interval for mu
print(muhat + c(-1,1)*qnorm(1-alpha/2)*sdmu)

# Do the same for sigma
sigmasqhat <- mu_sigmasqhat(x, mode, param1, param2, mu=mu)[2]
# Assuming mu is known, plug in sigmasqhat
sdsigmasq <- sqrt(varhat(mu, sigmasqhat, mode, param1, param2, est_mu=FALSE)/n)
# The 95% confidence interval for sigmasq
print(sigmasqhat + c(-1,1)*qnorm(1-alpha/2)*sdsigmasq)
readline("Press enter to continue to the next example.")

# Generate 10000 samples, each with sample size 1000,
# and test the coverage probability of the CIs above.
samples <- 10000
# 1000*10000 random samples from TN(mu, sigma)
x <- matrix(tmvtnorm::rtmvnorm(n*samples, mean=mu, sigma=sigma^2, lower=c(0), upper=c(Inf)), nrow=n, ncol=samples)

# One muhat estimate for each of the 10000 samples
muhats <- apply(x, 2, mu_sigmasqhat, mode=mode, param1=param1, param2=param2, sigmasq=sigma^2)[1,]
# One standard error for each of the 10000 samples, plugging in the estimates for mu and true sigmasq
ses_mu <- sapply(muhats, function(muhat){sqrt(varhat(muhat, sigma^2, mode, param1, param2, est_mu=TRUE)/n)})
# Coverage probability of the confidence intervals for mu, should be close to 0.95
print(mean(mu >= muhats - qnorm(1-alpha/2)*ses_mu & mu <= muhats + qnorm(1-alpha/2)*ses_mu))

# One sigmasqhat estimate for each of the 10000 samples
sigmasqhats <- apply(x, 2, mu_sigmasqhat, mode=mode, param1=param1, param2=param2, mu=mu)[2,]
# One standard error for each of the 10000 samples, plugging in the estimates for sigmasq and true mu
ses_sigmasq <- sapply(sigmasqhats, function(sigmasqhat){sqrt(varhat(mu, sigmasqhat^2, mode, param1, param2, est_mu=FALSE)/n)})
# Coverage probability of the confidence intervals for sigmasq, should be close to 0.95
print(mean(sigma^2 >= sigmasqhats - qnorm(1-alpha/2)*ses_sigmasq & sigma^2 <= sigmasqhats + qnorm(1-alpha/2)*ses_sigmasq))
readline("Press enter to continue to the next example.")

# Compare the asymptotic variance of the estimator to the empirical
# variance, as well as the Cram\'er-Rao lower bound.
trials <- 100
# Reshape x into a 1000*100*100 array
x <- array(x, dim=c(n,trials,samples/trials))
# The asymptotic variance for mu assuming sigmasq known
print(varhat(mu, sigma^2, mode, param1, param2, est_mu=TRUE))
# Empirical variance by averaging 100 variance estimates over 100 samples, each with sample size 1000
print(n*mean(apply(apply(x, c(2,3), mu_sigmasqhat, mode=mode, param1=param1, param2=param2, sigmasq=sigma^2)[1,,], 1, var)))
# The Cram\'er-Rao lower bound on the variance of unbiased estimators of mu
print(crbound_mu(mu, sigma))

# The asymptotic variance for sigmasq assuming mu known
print(varhat(mu, sigma^2, mode, param1, param2, est_mu=FALSE))
# Empirical variance by averaging 100 variance estimates over 100 samples, each with sample size 1000
print(n*mean(apply(apply(x, c(2,3), mu_sigmasqhat, mode=mode, param1=param1, param2=param2, mu=mu)[2,,], 1, var)))
# The Cram\'er-Rao lower bound on the variance of unbiased estimators of sigmasq
print(crbound_sigma(mu, sigma))
readline("Press enter to continue to the next example.")

## Plots from Yu et al (2019)
## The following plot reproduces Figure 1 from Yu et al (2019),
# which compares the asymptotic variances and efficiencies
# (relative to the Cram\'er-Rao lower bound) of our estimators with
# differnt h functions when estimating the mu parameter,
# assuming sigmasq=1 is known.
modes <- c("min_log_pow", "min_log_pow", "log_pow", "min_pow", "min_pow", "min_pow", "pow", "pow")
param1s <- c(1,1,1,1,1,2,1,2)
param2s <- c(1,2,-1,1,2,1,-1,-1)
h_names <- c("min(log(1+x),1)", "min(log(1+x),2)", "log(1+x)", "min(x,1)", "min(x,2)","min(x^2,1)","x","x^2")

mus <- seq(-4,4,length.out=100)
sigma <- 1
curves <- matrix(NA, nrow=length(h_names)+1, ncol=length(mus))
for (mu_i in 1:length(mus)){
  for (hi in 1:length(h_names))
    curves[hi,mu_i] <- varhat(mus[mu_i],sigma^2,modes[hi],param1s[hi],param2s[hi],est_mu=TRUE)
  curves[length(h_names)+1,mu_i] <- crbound_mu(mus[mu_i],sigma^2)
}

par(mfrow=c(1,2), mar=c(5,5,1,1), xpd=TRUE, bty="n")
plot(c(),c(), xlim=range(mus), ylim=log(range(curves)), xlab=expression(mu[0]), ylab="log(Asymptotic var)", cex.lab=1, cex.axis=1)
order <- c(1,2,3,4,7,6,5,8)
colors <- c("red","darkorange1","gold3","green","forestgreen",
            "hotpink","blue",gray.colors(1,start=0.7,end=0.7))[order]
for (hi in 1:length(h_names)){
  points(mus, log(curves[hi, ]), col = colors[hi], type="l", lwd=3)
}
points(mus, log(curves[length(h_names)+1,]), lty=2, type="l", lwd=3)
legend("topright", x.intersp = 0.2, y.intersp=0.7, inset = c(0.45,0.0), lty = c(rep(1,length(h_names)),2), ncol = 1, lwd=c(rep(3,length(h_names)),2), seg.len=1, bty = "n", text.width = .1, col=c(colors[match(1:length(h_names), order)],"black"), legend=c(h_names[match(1:length(h_names), order)],"C-R lower bound"), cex=1.2)

plot(c(),c(),xlim=range(mus),ylim=range(t(curves[length(h_names)+1,]/t(curves))), xlab=expression(mu[0]), ylab="Efficiency", cex.lab=1, cex.axis=1)
for (hi in 1:length(h_names)){
  points(mus, curves[length(h_names)+1,]/curves[hi,], col=colors[hi], type="l", lwd=3)
}
readline("Press enter to continue to the next example.")


## The following plot reproduces Figure 2 from Yu et al (2019),
# which compares the asymptotic variances and efficiencies
# (relative to the Cram\'er-Rao lower bound) of our estimators with
# differnt h functions when estimating the sigmasq parameter,
# assuming mu=0.5 is known.
sigmas <- exp(seq(log(0.1),log(10),length.out=100))
mu <- 0.5
curves2 <- matrix(NA, nrow=length(h_names)+1, ncol=length(sigmas))
for (sigma_i in 1:length(sigmas)){
  for (hi in 1:length(h_names))
    curves2[hi,sigma_i] <- varhat(mu,sigmas[sigma_i]^2,modes[hi],param1s[hi],param2s[hi],est_mu=FALSE)
  curves2[length(h_names)+1,sigma_i] <- crbound_sigma(mu,sigmas[sigma_i]^2)
}

par(mfrow=c(1,2), mar=c(5,5,1,1), xpd=TRUE, bty="n")
plot(c(),c(), xlim=range(sigmas), ylim=log(range(curves2)), xlab=expression(sigma[0]), ylab="log(Asymptotic var)", cex.lab=1, cex.axis=1)
order <- c(4,6,5,1,2,3,7,8)
colors <- c("red","darkorange1","gold3","green","forestgreen",
            "hotpink","blue",gray.colors(1,start=0.7,end=0.7))[order]
for (hi in 1:length(h_names)){
  points(sigmas, log(curves2[hi, ]), col = colors[hi], type="l", lwd=3)
}
points(sigmas, log(curves2[length(h_names)+1,]), type="l", lty=2, lwd=3)
legend("topright", x.intersp = 0.2, y.intersp=0.7, inset = c(0.5,0.3), lty = c(rep(1,length(h_names)),2), ncol = 1, lwd=c(rep(3,length(h_names)),2), seg.len=1, bty = "n", text.width = .1, col=c(colors[match(1:length(h_names), order)],"black"), legend=c(h_names[match(1:length(h_names), order)],"C-R lower bound"), cex=1)

plot(c(),c(),xlim=range(sigmas),ylim=range(t(curves2[length(h_names)+1,]/t(curves2))), xlab=expression(sigma[0]), ylab="Efficiency", cex.lab=1, cex.axis=1)
for (hi in 1:length(h_names)){
  points(sigmas, curves2[length(h_names)+1,]/curves2[hi,], col=colors[hi], type="l", lwd=3)
}
print("Demo completed.")
