/* Compile with R CMD SHLIB arms.c domain.c genscore.c sampling.c set_ops.c simplex_genscore.c tests.c utils.c -o genscore.so
 Update in 1002 from 0928: added sampler for the general a/b setting; in gibbs_one_round with gamma FALSE, A = Theta[j] now becomes Theta[j] * 2, reflecting the fact that the linear part in the formulation of the exp density changed from Theta%*%sqrt(x) to (Theta%*%((x^1/2-1)/(1/2))), as in the general a/b setting
 Thus, 1002.rsqrt_gibbs_one_round(...,Theta,...) is the same as 0928.rsqrt_gibbs_one_round(...,Theta*2,...)
 Update in genscore on 20190415: renamed rsqrt_gibbs_one_round to rexp_gamma_reject and renamed rsqrt_gibbs_one_round_ab to rab_arms
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <sys/param.h>
#include <math.h>
#include <stdio.h>
#include <Rmath.h>
#include <assert.h>
#include <R_ext/BLAS.h>
#include "utils.h"
#include "genscore.h"

double EPS = 1-1E-14;

// elts for centered Gaussian
void elts_gauss_c(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier){
	//// d is here although center, because we can reuse its calculation in forming g1
	int n = *nIn, p = *pIn, j,k,l, jn, lp, ln, jpp;
	// g1 <- crossprod(hpdx, x)/n; d <- colMeans(hdx); diag(g1) <- diag(g1) + d
	for (j=0; j<p; j++){ // Column
		d[j] = sum(n, hdx+j*n) / n;
		for (k=0; k<p; k++){ // Column
			g1[j*p+k] = in_order_dot_prod(n, hpdx+(k*n), x+j*n) / n;
		}
		g1[j*p+j] += d[j];
	}
	// Gamma is the big matrix flattened (or its upper-left flattened for non-profiled non-centered)
	// Will be modified later for profiled non-centered
	for (j=0; j<p; j++){
		jpp = j*p*p; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=l; k<p; k++){
				Gamma[jpp+lp+k] = Gamma[jpp+k*p+l] = in_order_tri_dot_prod(n, x+k*n, x+ln, hdx+jn) / n;
			}
			diagonals_with_multiplier[j*p+l] = *diagonal_multiplier * Gamma[jpp+lp+l]; //Gamma[jpp+lp+l] *= *diagonal_multiplier;
		}
	}
}

// elts for non-centered and non-profiled Gaussian; calls elts_gauss_c
void elts_gauss_np(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	elts_gauss_c(nIn, pIn, hdx, hpdx, x, g1, d, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	// Gamma12flat <- -crossprod(x,hdx)/n
	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = -in_order_dot_prod(n, x+k*n, hdx+j*n) / n;
		}
	}
	// g2 = -colMeans(hpdx)
	for (j=0; j<p; j++){ // Column
		g2[j] = -sum(n, hpdx+j*n) / n;
	}
}

// Takes the full Gamma and gs and convert them into Gamma11 - Gamma12 %*% Gamma22^(-1) %*% Gamma12 and g1 - Gamma12 %*% Gamma22^(-1) %*% g2
// eta can be solved by Gamma22^(-1) %*% (g2 - t(Gamma12) %*% vec(Khat)); namely newg2 - newGamma12 %*% vec(Khat) where newg2 and newGamma12 are g2 and Gamma12 returned by this function (to reuse variables for memory considerations)
void make_profile(int *pIn, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonals_with_multiplier){
	int p = *pIn, j, k, l;
	double *original_Gamma12 = (double*)malloc(p*p*sizeof(double));
	// Store Gamma12 in original_Gamma12 and replace Gamma12 <- Gamma12Gamma22inv = Gamma12 %*% Gamma22^(-1) flattened
	// Do this because we only need to return Gamma12Gamma22inv (to get eta) in the end, not Gamma12, but we still need the original Gamma12 to simplify calculations for Gamma below
	// Recall Gamma22 = d
	for (j=0; j<p; j++)
		for (k=0; k<p; k++){
			original_Gamma12[j*p+k] = Gamma12[j*p+k];
			Gamma12[j*p+k] /= d[j];
		}
	// Gamma = Gamma11 - Gamma12 %*% Gamma22^(-1) %*% Gamma12', all flattened
	for (j=0; j<p; j++)
		for (l=0; l<p; l++){
			diagonals_with_multiplier[j*p+l] -= Gamma12[j*p+l]*original_Gamma12[j*p+l];
			for (k=0; k<p; k++)
				Gamma[j*p*p+l*p+k] -= Gamma12[j*p+k]*original_Gamma12[j*p+l];
		}
	free(original_Gamma12); // No longer need original_Gamma12
	// g1 = g1 - Gamma12 %*% Gamma22^(-1) %*% g2
	for (j=0; j<p; j++)
		for (k=0; k<p; k++)
			g1[j*p+k] -= Gamma12[j*p+k] * g2[j];
	// g2 <- g2/d, used to be called t2 but renamed here to eliminate another variable
	for (j=0; j<p; j++)
		g2[j] = g2[j] / d[j];
	// Gamma12 and g2 needed to recover eta_hat
}


void elts_ab_c(int *nIn, int *pIn, double *a, double *hdx, double *hpdx, double *x, double *g1, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier){
	// NOTE THAT b IS NOT NEEDED
	// NOTE THAT UNLIKE elts_gauss_c for gaussian, d IS NOT IN THE ARGUMENT LIST because we cannot reuse it in g1
	int n = *nIn, p = *pIn, j, k, l, jn, lp, ln, jpp;
	double A = *a;
	for (j=0; j<p; j++){
		jpp = j*p*p; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=l; k<p; k++){
				Gamma[jpp+lp+k] = 0;
				for (int i=0; i<n; i++)
					Gamma[jpp+lp+k] += pow(x[k*n+i]*x[ln+i], A)*hdx[jn+i]*pow(x[jn+i], 2*A-2);
				Gamma[jpp+lp+k] /= n;
				Gamma[jpp+k*p+l] = Gamma[jpp+lp+k];
			}
			diagonals_with_multiplier[j*p+l] = *diagonal_multiplier * Gamma[jpp+lp+l]; //Gamma[jpp+lp+l] *= *diagonal_multiplier;
			g1[j*p+l] = (in_order_tri_dot_prod_pow(n, hpdx+j*n, x+j*n, x+l*n, 1, A-1, A)
						 +(A-1)*in_order_tri_dot_prod_pow(n, hdx+j*n, x+j*n, x+l*n, 1, A-2, A)) / n;
		}
		g1[j*p+j] += A*in_order_dot_prod_pow(n, hdx+j*n, x+j*n, 1, 2*A-2) / n;
	}
}

// elts for generalized a/b setting; this function is based on elts_gauss_np
void elts_ab_np(int *nIn, int *pIn, double *a, double *b, double *hdx, double *hpdx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	double A = *a, B = *b;
	elts_ab_c(nIn, pIn, a, hdx, hpdx, x, g1, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	for (j=0; j<p; j++){
		d[j] = in_order_dot_prod_pow(n, hdx+j*n, x+j*n, 1, 2*B-2) / n;
		g2[j] -= (in_order_dot_prod_pow(n, hpdx+j*n, x+j*n, 1, B-1) + (B-1)*in_order_dot_prod_pow(n, hdx+j*n, x+j*n, 1, B-2)) / n;
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = -in_order_tri_dot_prod_pow(n, x+k*n, hdx+j*n, x+j*n, A, 1, A+B-2) / n;
		}
	}
}

void elts_exp_c(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier){
	// Note that although we treat the centered case here, we include d here because it's needed for g1; g2 is also included b/c it's convenient
	// Also used for elts_gamma_np since it's the same; g2 and d are just ignored
	int n = *nIn, p = *pIn, j, k, l, jn, lp, ln, jpp;
	for (j=0; j<p; j++){ // Column
		double tmp, tmp2 = 0; d[j] = 0; g2[j] = 0;
		for (k=0; k<p; k++) // Column
			g1[j*p+k] = 0;
		for (int i=0; i<n; i++){
			tmp2 = hdx[j*n+i]/x[j*n+i];
			tmp = (hpdx[j*n+i]-0.5*tmp2)/sqrt(x[j*n+i]);
			d[j] += tmp2; g2[j] -= tmp;
			for (k=0; k<p; k++)
				g1[j*p+k] += tmp*sqrt(x[k*n+i]);
		}
		for (k=0; k<p; k++)
			g1[j*p+k] /= n;
		d[j] /= n; g2[j] /= n;
		g1[j*p+j] += d[j]/2;
	}
	for (j=0; j<p; j++){
		jpp = j*p*p; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=l; k<p; k++){
				Gamma[jpp+lp+k] = 0;
				for (int i=0; i<n; i++)
					Gamma[jpp+lp+k] += sqrt(x[k*n+i]*x[ln+i])*hdx[jn+i]/x[jn+i];
				Gamma[jpp+lp+k] /= n;
				Gamma[jpp+k*p+l] = Gamma[jpp+lp+k];
			}
			diagonals_with_multiplier[j*p+l] = *diagonal_multiplier * Gamma[jpp+lp+l]; //Gamma[jpp+lp+l] *= *diagonal_multiplier;
		}
	}
}

// elts for sqrt exp graphical models; corresponds to a = b = 1/2
void elts_exp_np(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	elts_exp_c(nIn, pIn, hdx, hpdx, x, g1, g2, d, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = 0;
			for (int i=0; i<n; i++)
				Gamma12[j*p+k] -= sqrt(x[k*n+i])*hdx[j*n+i]/x[j*n+i];
			Gamma12[j*p+k] /= n;
		}
	}
}

// elts for gamma graphical models; corresponds to a = 1/2, b = 0
void elts_gamma_np(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_exp_c(nIn, pIn, hdx, hpdx, x, g1, g2, d, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	int n = *nIn, p = *pIn, j, k;
	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = 0;
			for (int i=0; i<n; i++)
				Gamma12[j*p+k] -= sqrt(x[k*n+i]/x[j*n+i])*hdx[j*n+i]/x[j*n+i];
			Gamma12[j*p+k] /= n;
		}
	}
	for (j=0; j<p; j++){ // Column
		g2[j] = 0; d[j] = 0;
		for (int i=0; i<n; i++){
			g2[j] -= (hpdx[j*n+i]-hdx[j*n+i]/x[j*n+i])/x[j*n+i];
			d[j] += hdx[j*n+i]/x[j*n+i]/x[j*n+i];
		}
		g2[j] /= n; d[j] /= n;
	}
}

void elts_loglog_c(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier,
				   double *logx, double *h_over_xsq, double *hp_over_x){
	// NOTE THAT b IS NOT NEEDED
	// NOTE THAT UNLIKE elts_gauss_c for gaussian, d IS NOT IN THE ARGUMENT LIST because we cannot reuse it in g1
	int n = *nIn, p = *pIn, j, k, l, jn, lp, ln, jpp;
	for (int i=0; i<n; i++)
		for (int j=0; j<p; j++) {
			int idx = j*n+i;
			logx[idx] = log(x[idx]);
			h_over_xsq[idx] = hdx[idx] / x[idx] / x[idx];
			hp_over_x[idx] = hpdx[idx] / x[idx];
		}
	for (j=0; j<p; j++){
		jpp = j*p*p; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=l; k<p; k++)
				Gamma[jpp+lp+k] = Gamma[jpp+k*p+l] = in_order_tri_dot_prod(n, logx+k*n, logx+ln, h_over_xsq+jn) / n;
			g1[j*p+l] = (in_order_dot_prod(n, hp_over_x+j*n, logx+l*n) -
						 in_order_dot_prod(n, h_over_xsq+j*n, logx+l*n)) / n;
		}
		d[j] = sum(n, h_over_xsq+j*n) / n;
		g1[j*p+j] += d[j];
	}
	for (j=0; j<p; j++) {
		double *Gamma_j_diag = Gamma + j*p*p;
		for (l=0; l<p; l++) {
			*diagonals_with_multiplier++ = *diagonal_multiplier * (*Gamma_j_diag);
			Gamma_j_diag += p+1;
		}
	}
}

void elts_loglog_np(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	double *logx = (double*)malloc(n*p*sizeof(double));
	double *h_over_xsq = (double*)malloc(n*p*sizeof(double));
	double *hp_over_x = (double*)malloc(n*p*sizeof(double));
	elts_loglog_c(nIn, pIn, hdx, hpdx, x, g1, d, Gamma, diagonal_multiplier, diagonals_with_multiplier, logx, h_over_xsq, hp_over_x);
	for (j=0; j<p; j++){
		g2[j] = d[j] - sum(n, hp_over_x+j*n) / n;
		for (k=0; k<p; k++)
			Gamma12[j*p+k] = -in_order_dot_prod(n, logx+k*n, h_over_xsq+j*n) / n;
	}
	free(logx); free(h_over_xsq); free(hp_over_x);
}

// elts for non-centered and profiled; calls elts_gauss_np and make_profile
void elts_gauss_p(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_gauss_np(nIn, pIn, hdx, hpdx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_ab and make_profile
void elts_ab_p(int *nIn, int *pIn, double *a, double *b, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_ab_np(nIn, pIn, a, b, hdx, hpdx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_exp and make_profile
void elts_exp_p(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_exp_np(nIn, pIn, hdx, hpdx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_gamma and make_profile
void elts_gamma_p(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_gamma_np(nIn, pIn, hdx, hpdx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_gamma and make_profile
void elts_loglog_p(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_loglog_np(nIn, pIn, hdx, hpdx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

inline double shrink(double a, double b) {
	if (b < fabs(a)) {
		if (a > 0) return(a-b);
		else       return(a+b);
	} else {
		return(0.0);
	}
}

inline int lindx(int r, int c, int p) {
	int rr,cc;
	
	if (r < c){ rr = r; cc = c; }
	else      { rr = c; cc = r; }
	
	return( (p*rr)+cc - ((rr*(rr+1))/2) );
	
}

inline double set_KKT_bound(double bound, double tol){
	return ((bound * EPS >= 10*tol) ? (bound * (EPS)) : (bound - 10*tol));
}

// Loss for centered, or non-centered profiled estimator
// If used for refit eBIC, no penalty included
// Also works for an asymmetric K, assuming Gamma_K is "symmetric"
double loss_profiled (int p, double *Gamma_K, double *g_K, double *K, double *diagonals_with_multiplier, double lambda1){
	double crit1 = 0, crit2 = 0, crit3 = 0, crit4 = 0;
	int j,k;
	for (k=0;k<p;k++){ // Row
		crit1 -= in_order_dot_prod(p, K+k*p, g_K+k*p); // Sums over linear terms on Kji over j
	}
	if (diagonals_with_multiplier == NULL) // If refit, calculate the loss using the original diagonal entries instead of the l2-penalized loss
		for (k=0; k<p; k++)
			for (j=0; j<p; j++)
				crit2 += K[k*p+j]*K[k*p+j]*Gamma_K[k*p*p+j*p+j];
	else {
		for (k=0; k<p; k++){
			for (j=0; j<p; j++){
				crit2 += K[k*p+j]*K[k*p+j]*diagonals_with_multiplier[k*p+j]; // Quadratic terms on Kji
			}
			if (diagonals_with_multiplier != NULL) // If not refit (lambda1 = 0 if refit)
				crit4 += abs_sum(p, K+k*p) - fabs(K[k*p+k]); // sum(abs(K)) off diagonal
		}
	}
	for (k=0; k<p; k++){
		for (j=0; j<p-1; j++){
			crit3 += K[k*p+j]*in_order_dot_prod(p-j-1, K+k*p+j+1, Gamma_K+k*p*p+j*p+j+1);  // Cross term between Kjk and Kik, summed over i = (j+1):(p-1). Note that Sd2 is symmetric (Sd2[i,j,k]=Sd2[i,k,j]).
		}
	}
	return (crit1 + crit2/2 + crit3 + crit4 * lambda1);
}

void test_loss_profiled(double *crit, int *p, double *Gamma_K, double *g_K, double *K, double *diagonals_with_multiplier, double *lambda1) {
	*crit = loss_profiled (*p, Gamma_K, g_K, K, diagonals_with_multiplier, *lambda1);
}

// Loss for centered, or non-centered profiled estimator for UNTRUNCATED gaussian, assuming Gamma_K is of size p*p and g_K = c(diag(p))
// If used for refit eBIC, no penalty included
double loss_profiled_gauss (int p, double *Gamma_K, double *K, double *diagonals_with_multiplier, double lambda1){
	double crit1 = 0, crit2 = 0, crit3 = 0, crit4 = 0;
	int j,k;
	for (k=0;k<p;k++)
		crit1 -= K[k+k*p]; // g%*%c(K) = sum(diag(K))
	if (diagonals_with_multiplier == NULL) // If refit, calculate the loss using the original diagonal entries instead of the l2-penalized loss
		for (k=0; k<p; k++)
			for (j=0; j<p; j++)
				crit2 += K[k*p+j]*K[k*p+j]*Gamma_K[j*p+j];
	else {
		for (k=0; k<p; k++){
			for (j=0; j<p; j++){
				crit2 += K[k*p+j]*K[k*p+j]*diagonals_with_multiplier[j]; // Quadratic terms on Kji
			}
			if (diagonals_with_multiplier != NULL) // If not refit (lambda1 = 0 if refit)
				crit4 += abs_sum(p, K+k*p) - fabs(K[k*p+k]); // sum(abs(K)) off diagonal
		}
	}
	for (k=0; k<p; k++)
		for (j=0; j<p-1; j++)
			crit3 += K[k*p+j]*in_order_dot_prod(p-j-1, K+k*p+j+1, Gamma_K+j*p+j+1);  // Cross term between Kjk and Kik, summed over i = (j+1):(p-1). Note that Sd2 is symmetric (Sd2[i,j,k]=Sd2[i,k,j]).
	return (crit1 + crit2/2 + crit3 + crit4 * lambda1);
}

// Loss for non-centered non-profiled estimator
// If used for refit eBIC, no penalty included
double loss_full_penalized (int p, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *diagonals_with_multiplier, double lambda1, double lambda2){
	double crit = loss_profiled(p, Gamma_K, g_K, K, diagonals_with_multiplier, lambda1);
	for (int j=0; j<p; j++)
		crit += eta[j] * (in_order_dot_prod(p, Gamma_K_eta+j*p, K+j*p) - g_eta[j]);
	// Line above: -eta[j]*g_eta[j]: linear term on eta_j
	// eta[j] * in_order_dot_prod(p, Gamma_K_eta+j*p, K+j*p): Cross term between eta_j and Kkj, over k=0:(p-1)
	crit += in_order_tri_dot_prod(p, Gamma_eta, eta, eta) / 2; // Quadratic terms with eta_i; no interaction between different etas
	if (diagonals_with_multiplier != NULL) // If not refit (lambda2 = 0 if refit)
		crit += lambda2 * abs_sum(p, eta);
	return (crit);
}

void test_loss_full_penalized(double *crit, int *p, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *diagonals_with_multiplier, double *lambda1, double *lambda2) {
	*crit = loss_full_penalized(*p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, diagonals_with_multiplier, *lambda1, *lambda2);
}

// Loss for non-centered non-profiled estimator for UNTRUNCATED gaussian, assuming Gamma_K is of size p*p, g_K = c(diag(p)) and g_eta = numeric(p)
// If used for refit eBIC, no penalty included
double loss_full_penalized_gauss (int p, double *Gamma_K, double *Gamma_K_eta, double *K, double *eta, double *diagonals_with_multiplier, double lambda1, double lambda2){
	double crit = loss_profiled_gauss(p, Gamma_K, K, diagonals_with_multiplier, lambda1);
	for (int j=0; j<p; j++)
		crit += eta[j] * (in_order_dot_prod(p, Gamma_K_eta, K+j*p)); // Cross term between eta_j and Kkj, over k=0:(p-1); no linear term for eta for untruncated gaussian
	crit += in_order_dot_prod(p, eta, eta) / 2; // Since Gamma_eta=diag(p), the quadratic term one eta reduces to its squared l2 norm
	if (diagonals_with_multiplier != NULL) // If not refit
		crit += lambda2 * abs_sum(p, eta);
	return (crit);
}

void estimator_profiled (int *pIn, double *Gamma_K, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, int *exclude, double *diagonals_with_multiplier, int *gauss){
	// gauss: untruncated Gaussian if True
	int i,j;
	int p = *pIn;
	double lambda = *lambda1In;
	int tp=p*(p+1)/2;
	*converged=0;
	double sSum,s1,s2;
	double maxdiff=1.;
	
	double *oldK;
	oldK = (double*)malloc(tp*sizeof(double)); // Only need upper triangular to save space
	//double pen = 0;
	if (oldK == 0){
		Rprintf("Out of Memory!\n");
		return;
	}
	for (i = 0; i < p; i++){
		for (j = i; j < p; j++){
			oldK[lindx(i,j,p)] = K[i*p+j];
			//if (j != i) {pen += fabs(K[i*p+j]);} ////
		}
	}
	//*crit = loss_profiled(p, Gamma_K, g_K, K); ////
	//Rprintf("Initial crit %f \n", *crit+pen*2*lambda); ////
	
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		//pen = 0; ////
		
		// update diagonal elements of K
		for (i = 0; i < p; i++){
			int ip = ((*gauss) ? 0 : i*p), ipp = ip*p;
			s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+i*p);
			s1 += K[i*p+i] * Gamma_K[ipp+i*p+i] + ((*gauss) ? 1 : g_K[i*p+i]);
			if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
				sSum = Gamma_K[ipp+i*p+i];
			else
				sSum = diagonals_with_multiplier[ip+i];
			K[i*p+i]=s1/sSum; // Do not penalize diagonal terms, lambda = 0, no thresholding
			maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i,i,p)]-K[i*p+i]));
			oldK[lindx(i,i,p)] = K[i*p+i];
		}
		
		// update off-diagonal elements of K
		for (i = 0; i < (p-1); i++){
			for (j = i+1; j < p; j++){
				if (exclude != NULL && exclude[i*p+j])
					continue;
				int ip = ((*gauss) ? 0 : i*p), ipp = ip*p, jp = ((*gauss) ? 0 : j*p), jpp = jp*p;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+j*p); // -Sum over k of Gamma_i:ki,ji * K_ki -> (weight on Kji)
				s2 = -in_order_dot_prod(p, K+j*p, Gamma_K+jpp+i*p); // -Sum over k of Gamma_j:kj,ij * K_kj -> (weight on Kij)
				s1 += K[i*p+j]*Gamma_K[ipp+j*p+j]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j)
				s2 += K[j*p+i]*Gamma_K[jpp+i*p+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,j);
				if (!*gauss) {s1 += g_K[ip+j]; s2 += g_K[jp+i];}
				if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
					sSum = Gamma_K[ipp+j*p+j]+Gamma_K[jpp+i*p+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
				else
					sSum = diagonals_with_multiplier[ip+j]+diagonals_with_multiplier[jp+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
				K[i*p+j] = shrink((s1+s2)/sSum, 2*lambda/sSum); // Update Kij and Kji simultaneously by averaging; no need to divide by 2 since both numerator and denominator are doubled
				K[j*p+i] = K[i*p+j];
				maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, j, p)] - K[i*p+j])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[lindx(i,j,p)] = K[i*p+j]; // oldK contains upper triangular only
				//pen += fabs(K[i*p+j]); ////
			}
		}
		
		//*crit = loss_profiled(p, Gamma_K, g_K, K, lambda); ////
		//Rprintf("Iteration %d, crit = %f\n", iter, *crit+2*pen*lambda); ////
		
		if (maxdiff < *tol){
			*converged = 1;
			break;
		}
	}
	free(oldK);
	
}


void profiled(int *pIn, double *Gamma_K, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss){
	int p=*pIn;
	if (*is_refit){ // If refit, directly estimate with support restricted to exclude
		*lambda1In = 0;
		estimator_profiled (pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, NULL, gauss);
		if (*gauss)
			*crit = loss_profiled_gauss(p, Gamma_K, K, NULL, 0);
		else
			*crit = loss_profiled(p, Gamma_K, g_K, K, NULL, 0);
	}
	else {
		double KKT_bound = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol),
		KKT_bound_new = set_KKT_bound(*lambda1In, *tol);
		int first_time = 1, total_iters = 0, i,j;
		while (TRUE){
			if (!(first_time && KKT_bound > *lambda1In)){ // If first time and previous lambda is smaller, no need to check exclude since the support for the current lambda is necessarily a subset of the previous one
				int need_rerun = 0;
				// Only calculate for those not currently in the edge set
				for (i = 0; i < (p-1); i++){
					//grad[i*p+i] = g_K[i*p+i] - in_order_dot_prod(p, Gamma_K+i*p*p+i*p, K+i*p);
					for (j = i+1; j < p; j++){
						if (!exclude[i*p+j])
							continue;
						double grad = 0.0;
						if (*gauss)
							grad = (-in_order_dot_prod(p, Gamma_K+i*p, K+j*p) -
									in_order_dot_prod(p, Gamma_K+j*p, K+i*p) +
									(Gamma_K[i*p+i] - diagonals_with_multiplier[i]) * K[j*p+i] +
									(Gamma_K[j*p+j] - diagonals_with_multiplier[j]) * K[i*p+j]) / 2;
						else
							grad = (g_K[j*p+i] + g_K[i*p+j] -
									in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) -
									in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) +
									(Gamma_K[j*p*p+i*p+i] - diagonals_with_multiplier[j*p+i]) * K[j*p+i] +
									(Gamma_K[i*p*p+j*p+j]-diagonals_with_multiplier[i*p+j]) * K[i*p+j]) / 2;
						if (fabs(grad) > KKT_bound){
							need_rerun = 1; exclude[j*p+i] = 0; exclude[i*p+j] = 0;
						}
					}
				}//grad[p*p-1] = g_K[p*p-1] - in_order_dot_prod(p, Gamma_K+p*p*p+p*p, K+p*p);
				if (!first_time && !need_rerun)
					break;
			}
			estimator_profiled(pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, diagonals_with_multiplier, gauss);
			total_iters += *iters;
			first_time = 0; KKT_bound = KKT_bound_new;
		}
		estimator_profiled (pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, NULL, diagonals_with_multiplier, gauss);  // Run one more time to ensure correctness
		*iters = total_iters + *iters;
		if (*gauss)
			*crit = loss_profiled_gauss(p, Gamma_K, K, diagonals_with_multiplier, *lambda1In);
		else
			*crit = loss_profiled(p, Gamma_K, g_K, K, diagonals_with_multiplier, *lambda1In);
	}
}


void estimator_full_penalized (int *pIn, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, int *exclude, int *exclude_eta, double *diagonals_with_multiplier, int *gauss){
	int i,j;
	int p = *pIn;
	double lambda1 = *lambda1In, lambda2 = *lambda2In;
	int tp=p*(p+1)/2;
	*converged=0;
	double sSum,s1,s2;
	double maxdiff=1.;
	
	double *oldK, *oldeta;
	oldK = (double*)malloc(tp*sizeof(double)); // Only need upper triangular to save space
	oldeta = (double*)malloc(p*sizeof(double));
	if (oldK == 0 || oldeta == 0){
		Rprintf("Out of Memory!\n");
		return;
	}
	for (i = 0; i < p; i++){
		oldeta[i] = eta[i];
		for (j = i; j < p; j++){
			// Ignore diag if not estimating diagonals
			oldK[lindx(i,j,p)] = K[i*p+j];
		}
	}
	
	//*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta); ////
	//Rprintf("Crit %f \n", *crit); ////
	
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		
		// update diagonal elements of K
		for (i = 0; i < p; i++){
			int ip = ((*gauss) ? 0 : i*p), ipp = ip*p;
			s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+i*p);
			s1 += K[i*p+i] * Gamma_K[ipp+i*p+i] - Gamma_K_eta[ip+i]*eta[i] + ((*gauss) ? 1 : g_K[i*p+i]);
			if (diagonals_with_multiplier == NULL)
				sSum = Gamma_K[ipp+i*p+i];
			else
				sSum = diagonals_with_multiplier[ip+i];
			K[i*p+i] = s1/sSum; // Do not penalize diagonal terms, lambda1 = 0, no thresholding
			maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i,i,p)]-K[i*p+i]));
			oldK[lindx(i,i,p)] = K[i*p+i];
		}
		// update off-diagonal elements of K
		for (i = 0; i < (p-1); i++){
			for (j = i+1; j < p; j++){
				if (exclude != NULL && exclude[i*p+j])
					continue;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				int ip = ((*gauss) ? 0 : i*p), ipp = ip*p, jp = ((*gauss) ? 0 : j*p), jpp = jp*p;
				s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+j*p); // -Sum over k of Gamma_i:ki,ji * K_ki -> (weight on Kji)
				s2 = -in_order_dot_prod(p, K+j*p, Gamma_K+jpp+i*p); // -Sum over k of Gamma_j:kj,ij * K_kj -> (weight on Kij)
				s1 += K[i*p+j]*Gamma_K[ipp+j*p+j]-Gamma_K_eta[ip+j]*eta[i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j); additionally subtract the interaction btw Kji and etai
				s2 += K[j*p+i]*Gamma_K[jpp+i*p+i]-Gamma_K_eta[jp+i]*eta[j]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_s; consider s as (i,j)
				if (!*gauss) {s1 += g_K[i*p+j]; s2 += g_K[j*p+i];}
				if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
					sSum = Gamma_K[ipp+j*p+j]+Gamma_K[jpp+i*p+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
				else
					sSum = diagonals_with_multiplier[ip+j]+diagonals_with_multiplier[jp+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
				K[i*p+j] = shrink((s1+s2)/sSum, 2*lambda1/sSum); // Update Kij and Kji simultaneously by averaging; no need to divide by 2 since both numerator and denominator are doubled
				K[j*p+i] = K[i*p+j];
				maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, j, p)] - K[i*p+j])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[lindx(i,j,p)] = K[i*p+j]; // oldK contains upper triangular only
			}
		}
		for (i = 0; i < p; i++){
			if (exclude_eta != NULL && exclude_eta[i])
				continue;
			if (*gauss) {
				s1 = -in_order_dot_prod(p, K+i*p, Gamma_K_eta);
				sSum = 1;
			}
			else {
				s1 = g_eta[i] - in_order_dot_prod(p, K+i*p, Gamma_K_eta+i*p);
				sSum = Gamma_eta[i];
			}
			eta[i] = shrink(s1/sSum, lambda2/sSum);
			maxdiff = fmax2(maxdiff, fabs(oldeta[i] - eta[i]));
			oldeta[i] = eta[i];
		}
		
		//*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta); ////
		//Rprintf("Crit %f \n", *crit); ////
		
		if (maxdiff < *tol){
			*converged = 1;
			break;
		}
	}
	
	free(oldK);
	free(oldeta);
}


void full(int *pIn, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, int *exclude_eta, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss)
{
	int p=*pIn;
	if (*is_refit){ // If refit, directly estimate with support restricted to exclude
		*lambda1In = *lambda2In = 0;
		estimator_full_penalized (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta, NULL, gauss);
		if (*gauss)
			*crit = loss_full_penalized_gauss(p, Gamma_K, Gamma_K_eta, K, eta, NULL, 0, 0);
		else
			*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, NULL, 0, 0);
	}
	else {
		double KKT_bound1 = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol),
		KKT_bound1new = set_KKT_bound(*lambda1In, *tol);
		int first_time = 1, total_iters = 0, i,j;
		while (TRUE){
			if (!(first_time && KKT_bound1 > *lambda1In)){ // If first time and previous lambda is smaller, no need to check exclude since the support for the current lambda is necessarily a subset of the previous one
				int need_rerun = 0;
				// Only calculate for those not currently in the edge set
				for (i = 0; i < (p-1); i++){
					for (j = i+1; j < p; j++){
						if (!exclude[i*p+j])
							continue;
						double grad = 0.0;
						if (*gauss)
							grad = ((-in_order_dot_prod(p, Gamma_K+i*p, K+j*p) -
									 in_order_dot_prod(p, Gamma_K+j*p, K+i*p) +
									 (Gamma_K[i*p+i]-diagonals_with_multiplier[i]) * K[j*p+i] +
									 (Gamma_K[j*p+j] - diagonals_with_multiplier[j]) * K[i*p+j] -
									 Gamma_K_eta[j] * eta[i] - Gamma_K_eta[i]*eta[j]))/2;
						else
							grad = ((g_K[j*p+i] + g_K[i*p+j] -
									 in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) -
									 in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) +
									 (Gamma_K[j*p*p+i*p+i] - diagonals_with_multiplier[j*p+i]) * K[j*p+i] +
									 (Gamma_K[i*p*p+j*p+j]-diagonals_with_multiplier[i*p+j]) * K[i*p+j] -
									 Gamma_K_eta[i*p+j] * eta[i] - Gamma_K_eta[j*p+i] * eta[j]))/2;
						if (fabs(grad) > KKT_bound1){  // Currently does not check the KKT conditions for eta, as a final round will be run
							need_rerun += 1; exclude[j*p+i] = 0; exclude[i*p+j] = 0;
						}
					}
				}
				if (!first_time && !need_rerun)
					break;
			}
			estimator_full_penalized (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta, diagonals_with_multiplier, gauss);
			total_iters += *iters;
			first_time = 0; KKT_bound1 = KKT_bound1new;
		}
		estimator_full_penalized (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, NULL, NULL, diagonals_with_multiplier, gauss);  // Run one more time to ensure correctness
		*iters = total_iters + *iters;
		if (*gauss)
			*crit = loss_full_penalized_gauss(p, Gamma_K, Gamma_K_eta, K, eta, diagonals_with_multiplier, *lambda1In, *lambda2In);
		else
			*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, diagonals_with_multiplier, *lambda1In, *lambda2In);
	}
}

/// Asymmetric

void estimator_profiled_asymm (int *pIn, double *Gamma_K, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, int *exclude, double *diagonals_with_multiplier, int *gauss){
	int i,j;
	int p = *pIn;
	double lambda = *lambda1In;
	int tp=p*p;
	*converged=0;
	double sSum,s;
	double maxdiff=1.;
	double *oldK;
	oldK = (double*)malloc(tp*sizeof(double)); // Only need upper triangular to save space
	//double pen = 0;
	if (oldK == 0){
		Rprintf("Out of Memory!\n");
		return;
	}
	for (i = 0; i < p; i++){
		for (j = 0; j < p; j++){ // Kji
			oldK[i*p+j] = K[i*p+j];
			//if (j != i) {pen += fabs(K[i*p+j]);} ////
		}
	}
	
	//*crit = loss_profiled(p, Gamma_K, g_K, K); ////
	//Rprintf("Initial crit %f \n", *crit+pen*2*lambda); ////
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		// update all elements of K, off-diagonal and also diagonal
		for (i = 0; i < p; i++){ // Kji
			for (j = 0; j < p; j++){
				if (i != j && exclude != NULL && exclude[i*p+j])
					continue;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				int ip = ((*gauss) ? 0 : i*p), ipp = ip*p;
				s = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+j*p); // -Sum over k of Gamma_i:ki,ji * K_ki -> (weight on Kji)
				s += K[i*p+j]*Gamma_K[ipp+j*p+j]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i);
				s += ((*gauss) ? (i==j) : g_K[i*p+j]);
				if (diagonals_with_multiplier == NULL)  // For refitting without l2 penalty
					sSum = Gamma_K[ipp+j*p+j]; // Gamma_{ji,ji}
				else
					sSum = diagonals_with_multiplier[ip+j]; // Gamma_{ji,ji}
				if (i != j)
					K[i*p+j] = shrink(s/sSum, lambda/sSum); // shrink
				else
					K[i*p+j] = s/sSum; // no shrinkage for diagonal
				maxdiff = fmax2(maxdiff,fabs(oldK[i*p+j] - K[i*p+j])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[i*p+j] = K[i*p+j]; // oldK contains upper triangular only
				//pen += fabs(K[i*p+j]); ////
			}
		}
		//*crit = loss_profiled(p, Gamma_K, g_K, K); ////
		//Rprintf("Iteration %d, crit = %f\n", iter, *crit+2*pen*lambda); ////
		if (maxdiff < *tol){
			*converged = 1;
			break;
		}
	}
	free(oldK);
}


void profiled_asymm (int *pIn, double *Gamma_K, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss){
	int p=*pIn;
	if (*is_refit){ // If refit, directly estimate with support restricted to exclude
		*lambda1In = 0;
		estimator_profiled_asymm (pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, NULL, gauss);
		if (*gauss)
			*crit = loss_profiled_gauss(p, Gamma_K, K, NULL, 0);
		else
			*crit = loss_profiled(p, Gamma_K, g_K, K, NULL, 0);
	}
	else {
		double KKT_bound = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol),
		KKT_bound_new = set_KKT_bound(*lambda1In, *tol);
		int first_time = 1, total_iters = 0, i,j;
		while (TRUE){
			if (!(first_time && KKT_bound > *lambda1In)){
				int need_rerun = 0;
				// Only calculate for those not currently in the edge set
				for (i = 0; i < p; i++){ // Kji
					for (j = 0; j < p; j++){
						if (i==j || !exclude[i*p+j])
							continue;
						double grad = 0.0;
						if (*gauss)
							grad = (-in_order_dot_prod(p, Gamma_K+j*p, K+i*p) +
							(Gamma_K[j*p+j]-diagonals_with_multiplier[j])*K[i*p+j]);
						else
							grad = (g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) +
							(Gamma_K[i*p*p+j*p+j]-diagonals_with_multiplier[i*p+j]) * K[i*p+j]);
						if (fabs(grad) > KKT_bound){
							need_rerun = 1; exclude[i*p+j] = 0;
						}
					}
				}
				if (!first_time && !need_rerun)
					break;
			}
			estimator_profiled_asymm (pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, diagonals_with_multiplier, gauss);
			total_iters += *iters;
			first_time = 0; KKT_bound = KKT_bound_new;
		}
		estimator_profiled_asymm (pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, NULL, diagonals_with_multiplier, gauss);  // Run one more time to ensure correctness
		*iters = total_iters + *iters;
		if (*gauss)
			*crit = loss_profiled_gauss(p, Gamma_K, K, diagonals_with_multiplier, *lambda1In);
		else
			*crit = loss_profiled(p, Gamma_K, g_K, K, diagonals_with_multiplier, *lambda1In);
	}
}

void estimator_full_penalized_asymm (int *pIn, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, int *exclude, int *exclude_eta, double *diagonals_with_multiplier, int *gauss){
	int i,j;
	int p = *pIn;
	double lambda1 = *lambda1In, lambda2 = *lambda2In;
	int tp=p*p;
	*converged=0;
	double sSum,s;
	double maxdiff=1.;
	double *oldK, *oldeta;
	oldK = (double*)malloc(tp*sizeof(double)); // Only need upper triangular to save space
	oldeta = (double*)malloc(p*sizeof(double));
	if (oldK == 0 || oldeta == 0){
		Rprintf("Out of Memory!\n");
		return;
	}
	for (i = 0; i < p; i++){
		oldeta[i] = eta[i];
		for (j = 0; j < p; j++){ // Kji
			oldK[i*p+j] = K[i*p+j];
		}
	}

	//*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta); ////
	//Rprintf("Crit %f \n", *crit); ////
	
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		// update all elements of K, diagonal or off-diagonal
		for (i = 0; i < p; i++){
			for (j = 0; j < p; j++){
				if (i != j && exclude != NULL && exclude[i*p+j])
					continue;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				int ip = ((*gauss) ? 0 : i*p), ipp = ip*p;
				s = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+j*p); // -Sum over k of Gamma_i:ki,ji * K_ki -> (weight on Kji)
				s += K[i*p+j]*Gamma_K[ipp+j*p+j]-Gamma_K_eta[ip+j]*eta[i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j); additionally subtract the interaction btw Kji and etai
				s += ((*gauss) ? (i==j) : g_K[i*p+j]);
				if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
					sSum = Gamma_K[ipp+j*p+j]; // Gamma_{ji,ji}
				else
					sSum = diagonals_with_multiplier[ip+j]; // Gamma_{ji,ji}
				if (i != j)
					K[i*p+j] = shrink(s/sSum, lambda1/sSum); // shrinkage
				else
					K[i*p+j] = s/sSum; // no shrinkage for diagonal terms
				maxdiff = fmax2(maxdiff,fabs(oldK[i*p+j] - K[i*p+j])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[i*p+j] = K[i*p+j]; // oldK contains upper triangular only
			}
		}
		for (i = 0; i < p; i++){
			if (exclude_eta != NULL && exclude_eta[i])
				continue;
			if (*gauss) {
				s = - in_order_dot_prod(p, K+i*p, Gamma_K_eta);
				sSum = 1;
			}
			else {
				s = g_eta[i] - in_order_dot_prod(p, K+i*p, Gamma_K_eta+i*p);
				sSum = Gamma_eta[i];
			}
			eta[i] = shrink(s/sSum, lambda2/sSum);
			maxdiff = fmax2(maxdiff, fabs(oldeta[i] - eta[i]));
			oldeta[i] = eta[i];
		}
		//*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta); ////
		//Rprintf("Crit %f \n", *crit); ////
		
		if (maxdiff < *tol){
			*converged = 1;
			break;
		}
	}
	free(oldK);
	free(oldeta);
}

void full_asymm(int *pIn, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, int *exclude_eta, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss)
{
	int p=*pIn;
	if (*is_refit){ // If refit, directly estimate with support restricted to exclude
		*lambda1In = *lambda2In = 0;
		estimator_full_penalized_asymm (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta, NULL, gauss);
		if (*gauss)
			*crit = loss_full_penalized_gauss(p, Gamma_K, Gamma_K_eta, K, eta, NULL, 0, 0);
		else
			*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, NULL, 0, 0);
	}
	else {
		double KKT_bound1 = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol),
		KKT_bound1new = set_KKT_bound(*lambda1In, *tol);
		int first_time = 1, total_iters = 0, i,j;
		while (TRUE){
			if (!(first_time && KKT_bound1 > *lambda1In)){
				int need_rerun = 0;
				// Only calculate for those not currently in the edge set
				for (i = 0; i < p; i++){
					for (j = 0; j < p; j++){
						if (i==j || !exclude[i*p+j])
							continue;
						double grad = 0.0;
						if (*gauss)
							grad = (-in_order_dot_prod(p, Gamma_K+j*p, K+i*p) +
							(Gamma_K[j*p+j]-diagonals_with_multiplier[j]) * K[i*p+j] -
							Gamma_K_eta[j]*eta[i]);
						else
							grad = (g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) +
							(Gamma_K[i*p*p+j*p+j]-diagonals_with_multiplier[i*p+j]) * K[i*p+j] -
							Gamma_K_eta[i*p+j] * eta[i]);
						if (fabs(grad) > KKT_bound1){ // Currently does not check the KKT conditions for eta, as a final round will be run
							need_rerun += 1; exclude[i*p+j] = 0;
						}
					}
				}
				if (!first_time && !need_rerun)
					break;
			}
			estimator_full_penalized_asymm (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta, diagonals_with_multiplier, gauss);
			total_iters += *iters;
			first_time = 0; KKT_bound1 = KKT_bound1new;
		}
		estimator_full_penalized_asymm (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, NULL, NULL, diagonals_with_multiplier, gauss);  // Run one more time to ensure correctness
		*iters = total_iters + *iters;
		if (*gauss)
			*crit = loss_full_penalized_gauss(p, Gamma_K, Gamma_K_eta, K, eta, diagonals_with_multiplier, *lambda1In, *lambda2In);
		else
			*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, diagonals_with_multiplier, *lambda1In, *lambda2In);
	}
}




