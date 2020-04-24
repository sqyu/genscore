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
#include "simplex_genscore.h"

void elts_loglog_simplex_c(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, int *sum_to_zero, double *g_K, double *Gamma_K, double *Gamma_K_jp, double *Gamma_eta, double *Gamma_eta_jp, double *diagonal_multiplier, double *diagonals_with_multiplier, double *logx, double *h_over_xsq_nop, double *minus_h_over_x_xp_nop, double *sum_h_over_xmsq, double *hp_over_x_nop, double *sum_hp_over_xm, double *mean_sum_h_over_xmsq){
	// NOTE THAT b IS NOT NEEDED
	// NOTE THAT UNLIKE elts_gauss_c for gaussian, d IS NOT IN THE ARGUMENT LIST because we cannot reuse it in g1
	int n = *nIn, p = *pIn, j, k, l, jn, lp, ln, pp = p*p, jpp;
	for (int i=0; i<n; i++) {
		int pi = (p-1)*n+i;
		sum_h_over_xmsq[i] = sum_hp_over_xm[i] = 0;
		for (int j=0; j<p-1; j++) {
			int idx = j*n+i;
			logx[idx] = log(x[idx]);
			h_over_xsq_nop[idx] = hdx[idx] / x[idx] / x[idx];
			minus_h_over_x_xp_nop[idx] = -hdx[idx] / x[idx] / x[pi];
			hp_over_x_nop[idx] = hpdx[idx] / x[idx];
			sum_h_over_xmsq[i] += hdx[idx];
			sum_hp_over_xm[i] += hpdx[idx];
		}
		logx[pi] = log(x[pi]);
		sum_h_over_xmsq[i] /= x[pi] * x[pi];
		sum_hp_over_xm[i] /= x[pi];
	}
	for (j=0; j<p-1; j++){
		jpp = j*pp; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=l; k<p; k++) {
				Gamma_K[jpp+lp+k] = Gamma_K[jpp+k*p+l] = in_order_tri_dot_prod(n, logx+k*n, logx+ln, h_over_xsq_nop+jn) / n;
				Gamma_K_jp[jpp+lp+k] = Gamma_K_jp[jpp+k*p+l] = in_order_tri_dot_prod(n, logx+k*n, logx+ln, minus_h_over_x_xp_nop+jn) / n;
			}
			g_K[j*p+l] = (in_order_dot_prod(n, hp_over_x_nop+j*n, logx+l*n) -
						 in_order_dot_prod(n, h_over_xsq_nop+j*n, logx+l*n)) / n;
		}
		Gamma_eta[j] = sum(n, h_over_xsq_nop+j*n) / n;
		g_K[j*p+j] += Gamma_eta[j];
	}
	jpp = (p-1)*pp;
	for (l=0; l<p; l++){
		lp = l*p; ln = l*n;
		for (k=l; k<p; k++)
			Gamma_K[jpp+lp+k] = Gamma_K[jpp+k*p+l] = in_order_tri_dot_prod(n, logx+k*n, logx+ln, sum_h_over_xmsq) / n;
		g_K[(p-1)*p+l] = -(in_order_dot_prod(n, sum_h_over_xmsq, logx+l*n) +
					 in_order_dot_prod(n, sum_hp_over_xm, logx+l*n)) / n;
		if (l != (p-1)) {
			Gamma_eta_jp[l] = sum(n, minus_h_over_x_xp_nop+l*n) / n;
			g_K[(p-1)*p+l] += Gamma_eta_jp[l];
			g_K[l*p+(p-1)] += Gamma_eta_jp[l];
		}
	}
	*mean_sum_h_over_xmsq = sum(n, sum_h_over_xmsq) / n;
	g_K[(p-1)*p+(p-1)] += *mean_sum_h_over_xmsq;
	
	if (*sum_to_zero) {
		for (j=0; j<p; j++){
			eliminate_vec(pIn, g_K + j*p, j);
			eliminate_row_col(pIn, pIn, Gamma_K+j*pp, j, j);
		}
		for (j=0; j<p-1; j++)
			eliminate_row_col(pIn, pIn, Gamma_K_jp+j*pp, j, p-1);
	}
	
	for (j=0; j<p; j++) {
		double *Gamma_j_diag = Gamma_K + j*p*p;
		for (l=0; l<p; l++) {
			*diagonals_with_multiplier++ = *diagonal_multiplier * (*Gamma_j_diag);
			Gamma_j_diag += p+1;
		}
	}
}


void elts_loglog_simplex_np(int *nIn, int *pIn, double *hdx, double *hpdx, double *x, int *sum_to_zero, double *g_K, double *g_eta, double *Gamma_K, double *Gamma_K_eta, double *Gamma_K_jp, double *Gamma_Kj_etap, double *Gamma_Kp_etaj, double *Gamma_eta, double *Gamma_eta_jp, double *diagonal_multiplier, double *diagonals_with_multiplier) {
	int n = *nIn, p = *pIn, j, k;
	double *logx = (double*)malloc(n*p*sizeof(double));
	double *h_over_xsq_nop = (double*)malloc(n*(p-1)*sizeof(double));
	double *minus_h_over_x_xp_nop = (double*)malloc(n*(p-1)*sizeof(double));
	double *sum_h_over_xmsq = (double*)malloc(n*sizeof(double));
	double *hp_over_x_nop = (double*)malloc(n*p*sizeof(double));
	double *sum_hp_over_xm = (double*)malloc(n*sizeof(double));
	double *mean_sum_h_over_xmsq = (double*)malloc(sizeof(double));
	
	elts_loglog_simplex_c(nIn, pIn, hdx, hpdx, x, sum_to_zero, g_K, Gamma_K, Gamma_K_jp, Gamma_eta, Gamma_eta_jp, diagonal_multiplier, diagonals_with_multiplier, logx, h_over_xsq_nop, minus_h_over_x_xp_nop, sum_h_over_xmsq, hp_over_x_nop, sum_hp_over_xm, mean_sum_h_over_xmsq);
	
	for (j = 0; j < p-1; j++) {
		g_eta[j] = Gamma_eta[j] - sum(n, hp_over_x_nop+j*n) / n;
		for (k=0; k<p; k++) {
			Gamma_K_eta[j*p+k] = -in_order_dot_prod(n, logx+k*n, h_over_xsq_nop+j*n) / n;
			Gamma_Kj_etap[j*p+k] = Gamma_Kp_etaj[j*p+k] = -in_order_dot_prod(n, logx+k*n, minus_h_over_x_xp_nop+j*n) / n;
		}
	}
	g_eta[p-1] = *mean_sum_h_over_xmsq + sum(n, sum_hp_over_xm) / n;
	for (k=0; k<p; k++)
		Gamma_K_eta[(p-1)*p+k] = -in_order_dot_prod(n, logx+k*n, sum_h_over_xmsq) / n;
	Gamma_eta[p-1] = *mean_sum_h_over_xmsq;
	
	if (*sum_to_zero) {
		for (j=0; j<p-1; j++) {
			eliminate_vec(pIn, Gamma_Kj_etap + j*p, j);
			eliminate_vec(pIn, Gamma_Kp_etaj + j*p, p-1);
			eliminate_vec(pIn, Gamma_K_eta + j*p, j);
		}
		eliminate_vec(pIn, Gamma_K_eta + (p-1)*p, p-1);
	}
	free(logx); free(h_over_xsq_nop); free(minus_h_over_x_xp_nop);
	free(sum_h_over_xmsq); free(hp_over_x_nop); free(sum_hp_over_xm);
	free(mean_sum_h_over_xmsq);
}

double loss_loglog_simplex_profiled(int p, double *Gamma_K, double *g_K, double *Gamma_K_jp, double *K, double *diagonals_with_multiplier, double lambda1){
	// If sum_to_zero assumed FALSE, elements corresponding to diagonals should be 0 (or diagonals in K themselves should be 0) since the function calculates the entire Gamma_K*K - g_K
	double crit1 = 0, crit2 = 0, crit3 = 0, crit4 = 0, crit5 = 0;
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
			for (j=0; j<p; j++)
				crit2 += K[k*p+j]*K[k*p+j]*diagonals_with_multiplier[k*p+j]; // Quadratic terms on Kji
			if (diagonals_with_multiplier != NULL) { // If not refit (lambda1 = 0 if refit)
				crit4 += (abs_sum(p, K+k*p) - fabs(K[k*p+k])) *
				(k == p-1 ? (p-1) : 1); // sum(abs(K_{-k,k})); penalize last column by (p-1)*lambda1
				if (k != p-1)
					crit4 += fabs(K[k*p+p-1]) * ((p-1) - 1); // Penalize K_{p-1,k} by (p-1)*lambda1;
			}
		}
	}
	for (j=0; j<p-1; j++){
		for (k=0; k<p; k++){
			crit3 += K[k*p+j]*in_order_dot_prod(p-j-1, K+k*p+j+1, Gamma_K+k*p*p+j*p+j+1);  // Cross term between Kjk and Kik, summed over i = (j+1):(p-1). Note that Sd2 is symmetric (Sd2[i,j,k]=Sd2[i,k,j]).
			crit5 += K[(p-1)*p+k] * in_order_dot_prod(p, K+j*p, Gamma_K_jp+j*p*p+k*p); // Must be K[,j] %*% Gamma_K_jp_{j-th block} %*% K[,p-1] (order cannot be switched) since we subtracted the j-th row and the k-th column, not the j-th column and the k-th row. Gamma_K_jp not symmetric!!
		}
	}
	return (crit1 + crit2/2 + crit3 + crit4 * lambda1 + crit5);
}


void test_loss_loglog_simplex_profiled(double *crit, int *p, double *Gamma_K, double *g_K, double *Gamma_K_jp, double *K, double *diagonals_with_multiplier, double *lambda1){
	*crit = loss_loglog_simplex_profiled(*p, Gamma_K, g_K, Gamma_K_jp, K, diagonals_with_multiplier, *lambda1);
}

double loss_loglog_simplex_full_penalized (int p, double *Gamma_K, double *Gamma_K_eta, double *Gamma_K_jp, double *Gamma_Kj_etap, double *Gamma_Kp_etaj, double *Gamma_eta, double *Gamma_eta_jp, double *g_K, double *g_eta, double *K, double *eta, double *diagonals_with_multiplier, double lambda1, double lambda2){
	double crit1 = loss_loglog_simplex_profiled(p, Gamma_K, g_K, Gamma_K_jp, K, diagonals_with_multiplier, lambda1);
	double crit2 = 0, crit4 = 0, crit6 = 0;
	for (int j=0; j<p; j++)
		crit2 += eta[j] * (in_order_dot_prod(p, Gamma_K_eta+j*p, K+j*p) - g_eta[j]); // Cross term between eta_j and Kkj, over k=0:(p-1)
	double crit3 = eta[p-1] * in_order_dot_prod(p-1, eta, Gamma_eta_jp); // Interaction between eta_j and eta_p
	for (int j=0; j<p-1; j++)
		crit4 += eta[j] * in_order_dot_prod(p, Gamma_Kp_etaj+j*p, K+(p-1)*p) +
				eta[p-1] * in_order_dot_prod(p, Gamma_Kj_etap+j*p, K+j*p); // Interaction between eta_j and K_p, and between eta_p and K_j
	double crit5 = in_order_tri_dot_prod(p, Gamma_eta, eta, eta) / 2; // Quadratic terms with eta_i; no interaction between different etas (except for eta_p)
	if (diagonals_with_multiplier != NULL) { // If not refit (lambda2 = 0 if refit)
		crit6 = abs_sum(p, eta);
		crit6 += fabs(eta[p-1]) * ((p-1) - 1); // Penalize eta[p-1] by (p-1)*lambda2
	}
	return (crit1 + crit2 + crit3 + crit4 + crit5 + crit6 * lambda2);
}

void test_loss_loglog_simplex_full_penalized (double *crit, int *p, double *Gamma_K, double *Gamma_K_eta, double *Gamma_K_jp, double *Gamma_Kj_etap, double *Gamma_Kp_etaj, double *Gamma_eta, double *Gamma_eta_jp, double *g_K, double *g_eta, double *K, double *eta, double *diagonals_with_multiplier, double *lambda1, double *lambda2){
	*crit = loss_loglog_simplex_full_penalized (*p, Gamma_K, Gamma_K_eta, Gamma_K_jp, Gamma_Kj_etap, Gamma_Kp_etaj, Gamma_eta, Gamma_eta_jp, g_K, g_eta, K, eta, diagonals_with_multiplier, *lambda1, *lambda2);
}


void estimator_simplex_centered(int *pIn, int *sum_to_zero, double *Gamma_K, double *Gamma_K_jp, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, int *exclude, double *diagonals_with_multiplier){
	// Need to make sure all elements in Gamma_K and g_K corresponding to K are 0
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
	if (*sum_to_zero)
		for (int i = 0; i < p; i++)
			K[i*p+i] = 0;
	
	for (int i = 0; i < p; i++)
		for (int j = i; j < p; j++) // Off-diagonal
			oldK[lindx(i,j,p)] = K[i*p+j];
	
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		
		// Updates first p-1 diagonal entries, if not assuming rows and columns sum to 1
		if (!*sum_to_zero) {
			for (int i = 0; i < p - 1; i++){
				int ip = i*p, ipp = ip*p;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				s1 = -in_order_dot_prod(p, K+ip, Gamma_K+ipp+ip); // -Sum over k of Gamma_i:ki,ii * K_ki -> (weight on Kii)
				s1 += K[ip+i] * Gamma_K[ipp+ip+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,i); add back the (unnecessarily subtracted) term for (i,i), (i,i)
				s1 -= dot_prod_by_row(p, Gamma_K_jp+ipp+i, K+(p-1)*p); // Interaction of Kii with K(p-1) (K[,i] %*% Gamma_K_jp_{i-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
				s1 += g_K[ip+i];
				if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
					sSum = Gamma_K[ipp+ip+i]; // Gamma_{ii,ii}
				else
					sSum = diagonals_with_multiplier[ip+i]; // Gamma_{ii,ii}
				K[i*p+i] = s1/sSum;
				maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, i, p)] - K[i*p+i])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[lindx(i,i,p)] = K[i*p+i]; // oldK contains upper triangular only
			}
		}
		
		// Updates last diagonal entry, if not assuming rows and columns sum to 1
		if (!*sum_to_zero) {
			int ip = (p-1)*p, ipp = ip*p;
			s1 = -in_order_dot_prod(p, K+ip, Gamma_K+ipp+ip); // Gamma on K_{,p-1}
			s1 += K[ip+p-1]*Gamma_K[ipp+ip+p-1]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (p-1,p-1); add back the (unnecessarily subtracted) term for (p-1,p-1), (p-1,p-1)
			for (int k = 0; k < p-1; k++) // Interaction of K(p-1) with K0, ..., K(p-2)
				s1 -= in_order_dot_prod(p, K+k*p, Gamma_K_jp+k*p*p+ip); // Interaction of K_{,k} with K_{p-1,p-1} (K[,k] %*% Gamma_K_jp_{k-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
			s1 += g_K[ip+p-1];
			if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
				sSum = Gamma_K[ipp+ip+p-1]; // Gamma_{ji,ji} + Gamma_{ij,ij}
			else
				sSum = diagonals_with_multiplier[ip+p-1]; // Gamma_{ji,ji} + Gamma_{ij,ij}
			K[ip+(p-1)] = s1/sSum; // Update Kii
			maxdiff = fmax2(maxdiff, fabs(oldK[lindx(p-1, p-1, p)] - K[ip+(p-1)])); // ||theta^(t) - theta^(t-1)||_infinity
			oldK[lindx(p-1, p-1, p)] = K[ip+(p-1)];
		}
		
		
		// Updates off-diagonal elements of K, except for last row/column
		for (int i = 0; i < p - 2; i++){
			for (int j = i + 1; j < p - 1; j++){
				if (exclude != NULL && exclude[i*p+j])
					continue;
				int ip = i*p, ipp = ip*p, jp = j*p, jpp = jp*p;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+j*p); // -Sum over k of Gamma_i:ki,ji * K_ki -> (weight on Kji)
				s2 = -in_order_dot_prod(p, K+j*p, Gamma_K+jpp+i*p); // -Sum over k of Gamma_j:kj,ij * K_kj -> (weight on Kij)
				s1 += K[i*p+j]*Gamma_K[ipp+j*p+j]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j)
				s2 += K[j*p+i]*Gamma_K[jpp+i*p+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,j);
				s1 -= dot_prod_by_row(p, Gamma_K_jp+ipp+j, K+(p-1)*p); // Interaction of Kji with K(p-1) (K[,i] %*% Gamma_K_jp_{i-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
				s2 -= dot_prod_by_row(p, Gamma_K_jp+jpp+i, K+(p-1)*p); // Interaction of Kij with K(p-1) (K[,j] %*% Gamma_K_jp_{j-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
				s1 += g_K[ip+j];
				s2 += g_K[jp+i];
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
		
		// Updates off-diagonal elements of K, for last row and column
		for (int i = 0; i < p-1; i++){ // Row index of last column of Ki
			if (exclude != NULL && exclude[i*p+p-1])
				continue;
			int ip = i*p, ipp = ip*p, jp = (p-1)*p, jpp = jp*p;
			s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+(p-1)*p); // For K_{p-1,i}: Gamma on K_{,i}
			s2 = -in_order_dot_prod(p, K+(p-1)*p, Gamma_K+jpp+i*p); // For K_{i,p-1}: Gamma on K_{,p-1}
			s1 += K[i*p+p-1]*Gamma_K[ipp+(p-1)*p+p-1]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j)
			s2 += K[(p-1)*p+i]*Gamma_K[jpp+i*p+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,j);
			s1 -= dot_prod_by_row(p, Gamma_K_jp+ipp+p-1, K+(p-1)*p); // Interaction of K_{p-1,i} with K_{,p-1} (K[,i] %*% Gamma_K_jp_{i-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
			for (int k = 0; k < p-1; k++) // Interaction of K(p-1) with K0, ..., K(p-2)
				s2 -= in_order_dot_prod(p, K+k*p, Gamma_K_jp+k*p*p+i*p); // Interaction of K_{,k} with K_{i,p-1} (K[,k] %*% Gamma_K_jp_{k-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
			s1 += g_K[ip+p-1];
			s2 += g_K[jp+i];
			if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
				sSum = Gamma_K[ipp+(p-1)*p+p-1] + Gamma_K[jpp+i*p+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
			else
				sSum = diagonals_with_multiplier[ip+p-1] + diagonals_with_multiplier[jp+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
			K[i*p+(p-1)] = shrink((s1+s2)/sSum, 2 * (p-1) * lambda/sSum); // Update Kij and Kji simultaneously by averaging;
			K[(p-1)*p+i] = K[i*p+(p-1)];
			maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, p-1, p)] - K[i*p+(p-1)])); // ||theta^(t) - theta^(t-1)||_infinity
			oldK[lindx(i,p-1,p)] = K[i*p+(p-1)];
		}
		
		////Rprintf("%d: %f\n", *iters, loss_loglog_simplex_profiled(p, Gamma_K, g_K, Gamma_K_jp, K, NULL, 0));
		
		if (maxdiff < *tol){
			*converged = 1;
			break;
		}
	}
	free(oldK);
}

void simplex_centered(int *pIn, int *sum_to_zero, double *Gamma_K, double *Gamma_K_jp, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier){
	int p=*pIn;
	if (*is_refit){ // If refit, directly estimate with support restricted to exclude
		*lambda1In = 0;
		estimator_simplex_centered(pIn, sum_to_zero, Gamma_K, Gamma_K_jp, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, NULL);
		*crit = loss_loglog_simplex_profiled(p, Gamma_K, g_K, Gamma_K_jp, K, NULL, 0);
	}
	else {
		double KKT_bound = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol),
		KKT_bound_new = set_KKT_bound(*lambda1In, *tol);
		int first_time = 1, total_iters = 0;
		while (TRUE){
			if (!(first_time && KKT_bound > *lambda1In)){ // If first time and previous lambda is smaller, no need to check exclude since the support for the current lambda is necessarily a subset of the previous one
				int need_rerun = 0;
				// Off-diagonals except for last row/column
				for (int i = 0; i < p-2; i++){
					//grad[i*p+i] = g_K[i*p+i] - in_order_dot_prod(p, Gamma_K+i*p*p+i*p, K+i*p);
					for (int j = i+1; j < p-1; j++){
						if (!exclude[i*p+j])
							continue;
						double grad = (g_K[j*p+i] + g_K[i*p+j] -
									   in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) -
									   in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) -
									   dot_prod_by_row(p, Gamma_K_jp+i*p*p+j, K+(p-1)*p) -
									   dot_prod_by_row(p, Gamma_K_jp+j*p*p+i, K+(p-1)*p) +
									   (Gamma_K[j*p*p+i*p+i]-diagonals_with_multiplier[j*p+i]) * K[j*p+i] +
									   (Gamma_K[i*p*p+j*p+j]-diagonals_with_multiplier[i*p+j]) * K[i*p+j]) / 2;
						if (fabs(grad) > KKT_bound){
							need_rerun = 1; exclude[j*p+i] = 0; exclude[i*p+j] = 0;
						}
					}
				}
				// Off-diagonals in last row and column
				for (int i = 0; i < p-1; i++){
					if (!exclude[i*p+p-1])
						continue;
					double grad = (g_K[(p-1)*p+i] + g_K[i*p+p-1] -
								   in_order_dot_prod(p, Gamma_K+(p-1)*p*p+i*p, K+(p-1)*p) -
								   in_order_dot_prod(p, Gamma_K+i*p*p+(p-1)*p, K+i*p) -
								   dot_prod_by_row(p, Gamma_K_jp+i*p*p+p-1, K+(p-1)*p) +
								   (Gamma_K[(p-1)*p*p+i*p+i]-diagonals_with_multiplier[(p-1)*p+i]) * K[(p-1)*p+i] +
								   (Gamma_K[i*p*p+(p-1)*p+p-1]-diagonals_with_multiplier[i*p+p-1]) * K[i*p+p-1]) / 2;
					for (int k = 0; k < p-1; k++) // Interaction of K(p-1) with K0, ..., K(p-2)
						grad -= in_order_dot_prod(p, Gamma_K_jp+k*p*p+i*p, K+k*p) / 2; // Interaction of K_{,k} with K_{i,p-1} (K[,k] %*% Gamma_K_jp_{k-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
					if (fabs(grad) > KKT_bound * (p-1)){
						need_rerun = 1; exclude[(p-1)*p+i] = 0; exclude[i*p+p-1] = 0;
					}
				}
				if (!first_time && !need_rerun)
					break;
			}
			estimator_simplex_centered(pIn, sum_to_zero, Gamma_K, Gamma_K_jp, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, diagonals_with_multiplier);
			total_iters += *iters;
			first_time = 0; KKT_bound = KKT_bound_new;
		}
		estimator_simplex_centered (pIn, sum_to_zero, Gamma_K, Gamma_K_jp, g_K, K, lambda1In, tol, maxit, iters, converged, NULL, diagonals_with_multiplier);  // Run one more time to ensure correctness
		*iters = total_iters + *iters;
		*crit = loss_loglog_simplex_profiled(p, Gamma_K, g_K, Gamma_K_jp, K, diagonals_with_multiplier, *lambda1In);
	}
}



void estimator_simplex_full(int *pIn, int *sum_to_zero, double *Gamma_K, double *Gamma_K_eta, double *Gamma_K_jp, double *Gamma_Kj_etap, double *Gamma_Kp_etaj, double *Gamma_eta, double *Gamma_eta_jp, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, int *exclude, int *exclude_eta, double *diagonals_with_multiplier){
	// Need to make sure all elements in Gamma_K and g_K corresponding to K are 0
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
	if (*sum_to_zero)
		for (int i = 0; i < p; i++)
			K[i*p+i] = 0;
	
	for (int i = 0; i < p; i++)
		for (int j = i; j < p; j++) // Ignore diag
			oldK[lindx(i,j,p)] = K[i*p+j];

	for (int i = 0; i < p; i++)
		oldeta[i] = eta[i];
	
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		
		// Updates first p-1 diagonal entries, if not assuming rows and columns sum to 1
		if (!*sum_to_zero) {
			for (int i = 0; i < p - 1; i++){
				int ip = i*p, ipp = ip*p;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				s1 = -in_order_dot_prod(p, K+ip, Gamma_K+ipp+ip); // -Sum over k of Gamma_i:ki,ii * K_ki -> (weight on Kii)
				s1 += K[i*p+i]*Gamma_K[ipp+ip+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,i); add back the (unnecessarily subtracted) term for (i,i), (i,i)
				s1 -= dot_prod_by_row(p, Gamma_K_jp+ipp+i, K+(p-1)*p); // Interaction of Kji with K(p-1) (K[,i] %*% Gamma_K_jp_{i-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
				s1 -= Gamma_K_eta[ip+i] * eta[i]; // Interaction of Kii with etai
				s1 -= Gamma_Kj_etap[ip+i] * eta[p-1]; // Interaction of Kii with eta_{p-1}
				s1 += g_K[ip+i];
				if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
					sSum = Gamma_K[ipp+ip+i]; // Gamma_{ii,ji} + Gamma_{ij,ij}
				else
					sSum = diagonals_with_multiplier[ip+i]; // Gamma_{ii,ii}
				K[i*p+i] = s1/sSum; // Update Kii simultaneously by averaging; no need to divide by 2 since both numerator and denominator are doubled
				maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, i, p)] - K[i*p+i])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[lindx(i,i,p)] = K[i*p+i]; // oldK contains upper triangular only
			}
		}
		
		// Updates last diagonal entry, if not assuming rows and columns sum to 1
		if (!*sum_to_zero) {
			int ip = (p-1)*p, ipp = ip*p;
			s1 = -in_order_dot_prod(p, K+ip, Gamma_K+ipp+ip); // Gamma on K_{,p-1}
			s1 += K[ip+p-1]*Gamma_K[ipp+ip+p-1]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (p-1,p-1); add back the (unnecessarily subtracted) term for (p-1,p-1), (p-1,p-1)
			for (int k = 0; k < p-1; k++) // Interaction of K(p-1) with K0, ..., K(p-2)
				s1 -= in_order_dot_prod(p, K+k*p, Gamma_K_jp+k*p*p+ip); // Interaction of K_{,k} with K_{p-1,p-1} (K[,k] %*% Gamma_K_jp_{k-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
			s1 += g_K[ip+p-1];
			s1 -= Gamma_K_eta[ip+p-1] * eta[p-1]; // Interaction of K{p-1,p-1} with eta{p-1}
			for (int k=0; k <p-1;k++)
				s1 -= eta[k] * Gamma_Kp_etaj[k*p+p-1];  // Interaction of K{p-1,p-1} with eta0, ..., eta_{p-2}
			if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
				sSum = Gamma_K[ipp+ip+p-1]; // Gamma_{p-1,p-1}
			else
				sSum = diagonals_with_multiplier[ip+p-1]; // Gamma_{p-1,p-1}
			K[ip+(p-1)] = s1/sSum; // Update K{p-1,p-1}
			maxdiff = fmax2(maxdiff,fabs(oldK[lindx(p-1, p-1, p)] - K[ip+(p-1)])); // ||theta^(t) - theta^(t-1)||_infinity
			oldK[lindx(p-1, p-1, p)] = K[ip+(p-1)];
		}
		
		
		// update off-diagonal elements of K, except for the last row/column
		for (int i = 0; i < p - 2; i++){
			for (int j = i + 1; j < p - 1; j++){
				if (exclude != NULL && exclude[i*p+j])
					continue;
				int ip = i*p, ipp = ip*p, jp = j*p, jpp = jp*p;
				// Here: s1 = -[Gamma^T thetahat]_ji, s2 = ..._ij
				s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+j*p); // -Sum over k of Gamma_i:ki,ji * K_ki -> (weight on Kji)
				s2 = -in_order_dot_prod(p, K+j*p, Gamma_K+jpp+i*p); // -Sum over k of Gamma_j:kj,ij * K_kj -> (weight on Kij)
				s1 += K[i*p+j]*Gamma_K[ipp+j*p+j]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j)
				s2 += K[j*p+i]*Gamma_K[jpp+i*p+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,j);
				s1 -= dot_prod_by_row(p, Gamma_K_jp+ipp+j, K+(p-1)*p); // Interaction of Kji with K(p-1) (K[,i] %*% Gamma_K_jp_{i-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
				s2 -= dot_prod_by_row(p, Gamma_K_jp+jpp+i, K+(p-1)*p); // Interaction of Kij with K(p-1) (K[,j] %*% Gamma_K_jp_{j-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
				s1 -= Gamma_K_eta[ip+j] * eta[i]; // Interaction of Kji with etai
				s2 -= Gamma_K_eta[jp+i] * eta[j]; // Interaction of Kij with etaj
				s1 -= Gamma_Kj_etap[ip+j] * eta[p-1]; // Interaction of Kji with eta_{p-1}
				s2 -= Gamma_Kj_etap[jp+i] * eta[p-1]; // Interaction of Kij with eta_{p-1}
				s1 += g_K[ip+j];
				s2 += g_K[jp+i];
				if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
					sSum = Gamma_K[ipp+j*p+j]+Gamma_K[jpp+i*p+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
				else
					sSum = diagonals_with_multiplier[ip+j]+diagonals_with_multiplier[jp+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
				K[i*p+j] = shrink((s1+s2)/sSum, 2*lambda1/sSum); // Update Kij and Kji simultaneously by averaging; no need to divide by 2 since both numerator and denominator are doubled
				K[j*p+i] = K[i*p+j];
				maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, j, p)] - K[i*p+j])); // ||theta^(t) - theta^(t-1)||_infinity
				oldK[lindx(i,j,p)] = K[i*p+j]; // oldK contains upper triangular only
				//pen += fabs(K[i*p+j]); ////
			}
		}
		
		// Updates off-diagonal elements of K in the last row/column
		for (int i = 0; i < p-1; i++){ // Row index of last column of K
			if (exclude != NULL && exclude[i*p+p-1])
				continue;
			int ip = i*p, ipp = ip*p, jp = (p-1)*p, jpp = jp*p;
			s1 = -in_order_dot_prod(p, K+i*p, Gamma_K+ipp+(p-1)*p); // For K_{p-1,i}: Gamma on K_{,i}
			s2 = -in_order_dot_prod(p, K+(p-1)*p, Gamma_K+jpp+i*p); // For K_{i,p-1}: Gamma on K_{,p-1}
			s1 += K[i*p+p-1]*Gamma_K[ipp+(p-1)*p+p-1]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (j,i); add back the (unnecessarily subtracted) term for (i,j), (i,j)
			s2 += K[(p-1)*p+i]*Gamma_K[jpp+i*p+i]; // g(x)-Gamma_{:,s}^T theta + Gamma_ss theta_ss; consider s as (i,j);
			s1 -= dot_prod_by_row(p, Gamma_K_jp+ipp+p-1, K+(p-1)*p); // Interaction of K_{p-1,i} with K_{,p-1} (K[,i] %*% Gamma_K_jp_{i-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
			for (int k = 0; k < p-1; k++) // Interaction of K(p-1) with K0, ..., K(p-2)
				s2 -= in_order_dot_prod(p, K+k*p, Gamma_K_jp+k*p*p+i*p); // Interaction of K_{,k} with K_{i,p-1} (K[,k] %*% Gamma_K_jp_{k-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
			s1 += g_K[ip+p-1];
			s2 += g_K[jp+i];
			s1 -= Gamma_K_eta[ip+p-1] * eta[i]; // Interaction of K{p-1,i} with eta{i}
			s2 -= Gamma_K_eta[jp+i] * eta[p-1]; // Interaction of K{i,p-1} with eta{p-1}
			s1 -= Gamma_Kj_etap[ip+p-1] * eta[p-1]; // Interaction of K{p-1,i} with eta_{p-1}
			for (int k=0; k <p-1;k++)
				s2 -= eta[k] * Gamma_Kp_etaj[k*p+i];  // Interaction of K{i,p-1} with eta0, ..., eta_{p-2}
			if (diagonals_with_multiplier == NULL) // For refitting without l2 penalty
				sSum = Gamma_K[ipp+(p-1)*p+p-1] + Gamma_K[jpp+i*p+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
			else
				sSum = diagonals_with_multiplier[ip+p-1] + diagonals_with_multiplier[jp+i]; // Gamma_{ji,ji} + Gamma_{ij,ij}
			K[i*p+(p-1)] = shrink((s1+s2)/sSum, 2 * (p-1) * lambda1/sSum); // Update Kij and Kji simultaneously by averaging;
			K[(p-1)*p+i] = K[i*p+(p-1)];
			maxdiff = fmax2(maxdiff,fabs(oldK[lindx(i, p-1, p)] - K[i*p+(p-1)])); // ||theta^(t) - theta^(t-1)||_infinity
			oldK[lindx(i,p-1,p)] = K[i*p+(p-1)];
		}
		
		for (int i = 0; i < p-1; i++){
			if (exclude_eta != NULL && exclude_eta[i])
				continue;
			s1 = g_eta[i] - in_order_dot_prod(p, K+i*p, Gamma_K_eta+i*p);
			s1 -= eta[p-1] * Gamma_eta_jp[i]; // Interaction between eta[i] and eta[p-1]
			s1 -= in_order_dot_prod(p, Gamma_Kp_etaj+i*p, K+(p-1)*p); // Interaction between eta[i] and K_{,p-1}
			sSum = Gamma_eta[i];
			eta[i] = shrink(s1/sSum, lambda2/sSum);
			maxdiff = fmax2(maxdiff, fabs(oldeta[i] - eta[i]));
			oldeta[i] = eta[i];
		}
		
		if (exclude_eta == NULL || !exclude_eta[p-1]) { // eta[p-1]
			s1 = g_eta[p-1] - in_order_dot_prod(p, K+(p-1)*p, Gamma_K_eta+(p-1)*p);
			s1 -= in_order_dot_prod(p-1, eta, Gamma_eta_jp); // Interaction between eta[p-1] and eta[i]
			for (int k=0; k<p-1; k++)
				s1 -= in_order_dot_prod(p, Gamma_Kj_etap+k*p, K+k*p); // Interaction between eta[i] and K_{,p-1}
			sSum = Gamma_eta[p-1];
			eta[p-1] = shrink(s1/sSum, (p-1) * lambda2/sSum);
			maxdiff = fmax2(maxdiff, fabs(oldeta[p-1] - eta[p-1]));
			oldeta[p-1] = eta[p-1];
		}
		
		if (maxdiff < *tol){
			*converged = 1;
			break;
		}
	}
	free(oldK);
	free(oldeta);
}


void simplex_full(int *pIn, int *sum_to_zero, double *Gamma_K, double *Gamma_K_eta, double *Gamma_K_jp, double *Gamma_Kj_etap, double *Gamma_Kp_etaj, double *Gamma_eta, double *Gamma_eta_jp, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, int *exclude_eta, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier){
	int p=*pIn;
	if (*is_refit){ // If refit, directly estimate with support restricted to exclude
		*lambda1In = *lambda2In = 0;
		estimator_simplex_full(pIn, sum_to_zero, Gamma_K, Gamma_K_eta, Gamma_K_jp, Gamma_Kj_etap, Gamma_Kp_etaj, Gamma_eta, Gamma_eta_jp, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta, NULL);
		*crit = loss_loglog_simplex_full_penalized (p, Gamma_K, Gamma_K_eta, Gamma_K_jp, Gamma_Kj_etap, Gamma_Kp_etaj, Gamma_eta, Gamma_eta_jp, g_K, g_eta, K, eta, NULL, 0, 0);
	}
	else {
		double KKT_bound1 = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol),
		KKT_bound1_new = set_KKT_bound(*lambda1In, *tol);
		int first_time = 1, total_iters = 0;
		while (TRUE){
			if (!(first_time && KKT_bound1 > *lambda1In)){ // If first time and previous lambda is smaller, no need to check exclude since the support for the current lambda is necessarily a subset of the previous one
				int need_rerun = 0;
				// Only calculate for those not currently in the edge set
				for (int i = 0; i < p-2; i++){
					//grad[i*p+i] = g_K[i*p+i] - in_order_dot_prod(p, Gamma_K+i*p*p+i*p, K+i*p);
					for (int j = i+1; j < p-1; j++){
						if (!exclude[i*p+j])
							continue;
						double grad = (g_K[j*p+i] + g_K[i*p+j] -
									   in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) -
									   in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) -
									   dot_prod_by_row(p, Gamma_K_jp+i*p*p+j, K+(p-1)*p) -
									   dot_prod_by_row(p, Gamma_K_jp+j*p*p+i, K+(p-1)*p) +
									   (Gamma_K[j*p*p+i*p+i]-diagonals_with_multiplier[j*p+i]) * K[j*p+i] +
									   (Gamma_K[i*p*p+j*p+j]-diagonals_with_multiplier[i*p+j]) * K[i*p+j] -
									   Gamma_K_eta[i*p+j] * eta[i] - Gamma_K_eta[j*p+i] * eta[j] -
									   Gamma_Kj_etap[i*p+j] * eta[p-1] - Gamma_Kj_etap[j*p+i] * eta[p-1]
									   ) / 2;
						if (fabs(grad) > KKT_bound1){
							need_rerun = 1; exclude[j*p+i] = 0; exclude[i*p+j] = 0;
						}
					}
				}
				for (int i = 0; i < p-1; i++){ // Row index of last column of K
					if (!exclude[i*p+p-1])
						continue;
					double grad = (g_K[(p-1)*p+i] + g_K[i*p+p-1] -
								   in_order_dot_prod(p, Gamma_K+(p-1)*p*p+i*p, K+(p-1)*p) -
								   in_order_dot_prod(p, Gamma_K+i*p*p+(p-1)*p, K+i*p) -
								   dot_prod_by_row(p, Gamma_K_jp+i*p*p+p-1, K+(p-1)*p) +
								   (Gamma_K[(p-1)*p*p+i*p+i]-diagonals_with_multiplier[(p-1)*p+i]) * K[(p-1)*p+i] +
								   (Gamma_K[i*p*p+(p-1)*p+p-1]-diagonals_with_multiplier[i*p+p-1]) * K[i*p+p-1] -
								   Gamma_K_eta[i*p+p-1] * eta[i] - Gamma_K_eta[(p-1)*p+i] * eta[p-1] -
								   Gamma_Kj_etap[i*p+p-1] * eta[p-1] - dot_prod_by_row(p-1, Gamma_Kp_etaj+i, eta)
								   ) / 2;
					for (int k = 0; k < p-1; k++) // Interaction of K(p-1) with K0, ..., K(p-2)
						grad -= in_order_dot_prod(p, Gamma_K_jp+k*p*p+i*p, K+k*p) / 2; // Interaction of K_{,k} with K_{i,p-1} (K[,k] %*% Gamma_K_jp_{k-th block} %*% K[,p-1]); Gamma_K_jp not symmetric
					if (fabs(grad) > KKT_bound1 * (p-1)){
						need_rerun = 1; exclude[(p-1)*p+i] = 0; exclude[i*p+p-1] = 0;
					}
				}
				if (!first_time && !need_rerun)
					break;
			}
			estimator_simplex_full(pIn, sum_to_zero, Gamma_K, Gamma_K_eta, Gamma_K_jp, Gamma_Kj_etap, Gamma_Kp_etaj, Gamma_eta, Gamma_eta_jp, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta, diagonals_with_multiplier);
			total_iters += *iters;
			first_time = 0; KKT_bound1 = KKT_bound1_new;
		}
		estimator_simplex_full(pIn, sum_to_zero, Gamma_K, Gamma_K_eta, Gamma_K_jp, Gamma_Kj_etap, Gamma_Kp_etaj, Gamma_eta, Gamma_eta_jp, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, NULL, NULL, diagonals_with_multiplier);  // Run one more time to ensure correctness
		*iters = total_iters + *iters;
		*crit = loss_loglog_simplex_full_penalized (p, Gamma_K, Gamma_K_eta, Gamma_K_jp, Gamma_Kj_etap, Gamma_Kp_etaj, Gamma_eta, Gamma_eta_jp, g_K, g_eta, K, eta, diagonals_with_multiplier, *lambda1In, *lambda2In);
	}
}
