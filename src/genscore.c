/* Compile with R CMD SHLIB arms.c genscore.c  -o genscore.so
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
#include "arms.h"
#include "Score_gen.h"

int UNIT = 8;
double EPS = 1-1E-14;
#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */


/* ********************************************************************* */

struct ab_parm {
	double a, b, c, A, B, C;
};

/* ********************************************************************* */

double ab_density(double x, void *ab_data)
/* exp(A*x^a+B*x^(2a)+C*x^b)*/
{
	struct ab_parm *d;
	double y;
	
	/* cast voided pointer into pointer to struct normix_parm */
	if (x <= 0) return -INFINITY;
	d = ab_data;
	y = d->A*pow(x,d->a)+d->B*pow(x,2*d->a);
	y += d->b == 0.0 ? d->C*log(x) : d->C*pow(x,d->b);
	return y;
};

/* ********************************************************************* */

void samp_arms(int *n, double *samp, double *a, double *b, double *A, double *B, double *C)
{
	// In Gibbs sampling, this only works if samp is the point from the previous step, not when samp is fixed or completely random
	int neval, ninit = 5, npoint = 100, nsamp = 1, ncent = 0 ;
	double xinit[10]={1.0,3.0,10.0,17.0,20.0}, xl = 0, xr = 100.0;
	double xsamp[1], xcent[1], qcent[1];
	double convex = 1.;
	int dometrop = 1;
	
	/* initialise random number generator */
	
	double xprev = samp[0];//fabs(rnorm(1, 1));
	
	/* set up structures for each density function */
	
	if (*B >= 0){
		Rprintf("B needs to be negative! %4f provided.", *B);
		return;
	}
	
	struct ab_parm ab_data;
	
	/* initialise data needed by normal mixture density function */
	ab_data.a = *a;
	ab_data.b = *b;
	ab_data.A = *A;
	ab_data.B = *B;
	ab_data.C = *C;
	
	
	for(int i=0;i<*n;i++){
		int err = arms(xinit,ninit,&xl,&xr,ab_density,&ab_data,&convex,
					   npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
		if(err>0){
			Rprintf("error code in ARMS = %d\n",err);
			break;
			//exit(1);
		}
		if (xsamp[0] <= 0){
			Rprintf("Unable to generate ARMS sample. Got %f.\n", xsamp[0]);
			break;
		}
		samp[i] = xsamp[0];
	}
}


void one(int *n, double *x, double *A, double *B, double *C, int *max_iter){
	// Generates one sample from p \propto exp(A*sqrt(x)+B*x+C*log(x))
	// B must be < 0, C must be > -1; exponential if C = 0
	// M is so that exp(A*sqrt(x)+B*x) <= M*exp(-lambda*x)
	assert(*B < 0);
	double Aval = *A, Bval = *B, Cval = *C, M = 0, lambda = 0;
	if (Aval <= 0){
		M = 1; lambda = -Bval;
	} else {
		lambda = (Aval*Aval-8*Bval-Aval*sqrt(Aval*Aval-16*Bval))/8; M = exp(-Aval*Aval/4/(Bval+lambda));
	}
	//Rprintf("A=%4f, B=%4f, C=%4f, M=%4f, lambda=%4f", Aval, Bval, Cval, M, lambda);
	for (int i=0; i < *n; i++){
		x[i] = -1;
		for (int iter=0; iter < *max_iter; iter++){
			double y = (Cval == 0) ? -log(runif(0,1))/lambda : rgamma(Cval+1, 1/lambda);
			if (runif(0,1) < exp(Aval*sqrt(y)+(Bval+lambda)*y)/M){
				x[i] = y;
				break;
			}
		}
		assert(x[i] >= 0);
	}
}

inline double in_order_dot_prod(int len, double *l, double *r){
	// Computes dot product between vectors of l and r using unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0, UNIT = 8;
	while (i < len - len % UNIT) {
		total0 += l[0] * r[0]; total1 += l[1] * r[1];
		total2 += l[2] * r[2]; total3 += l[3] * r[3];
		total4 += l[4] * r[4]; total5 += l[5] * r[5];
		total6 += l[6] * r[6]; total7 += l[7] * r[7];
		i += UNIT; l += UNIT; r += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += (*(l++)) * (*(r++));
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

void rexp_gamma_reject(int *gamm, double *xinit, double *sqrtx, int *steps, int *m, double *Theta, double *Phi, int *max_iter, int *seed){
	// Runs Gibbs sampling for
	// sqrtx must be the square root of xinit, this is NOT checked
	int n = 1;
	if (*seed >= 0){
		srand(*seed);
		*seed += 1;
	}
	for (int iter = 0; iter < *steps; iter++){
		for (int j=0; j < *m; j++){
			double A = 2*(in_order_dot_prod(*m, Phi+j**m, sqrtx)-Phi[j+j**m]*sqrtx[j]);
			if (!*gamm)
				A += Theta[j]*2;
			double C = *gamm ? Theta[j] : 0;
			one(&n, xinit+j, &A, Phi+j**m+j, &C, max_iter);
			if (xinit[j] < 0){ // If failed to generate a sample
				xinit[0] = -1; // Signal R and stop generating
				return;
			}
			sqrtx[j] = sqrt(xinit[j]);
		}
	}
}


void rab_arms(double *a, double *b, double *xinit, double *xa, double *xb, int *steps, int *m, double *Theta, double *Phi, int *seed){
	//FILE *f;
	// Runs Gibbs sampling for exp(x^a%*%Phi%*%x^a/(2a) + Theta %*% (x^b-1)/b OR Theta %*% log(x))
	// xa = x^a, xb = x^b, not checked
	int n = 1;
	if (*seed >= 0){
		srand(*seed);
		*seed += 1;
	}
	for (int iter = 0; iter < *steps; iter++){
		for (int j=0; j < *m; j++){
			double A = (in_order_dot_prod(*m, Phi+j**m, xa)-Phi[j+j**m]*xa[j]) / *a;
			double B = Phi[j**m+j] / 2 / *a;
			double C = *b ? Theta[j] / *b : Theta[j];
			//Rprintf("A=%4f, B=%4f, C=%4f\n", A, B, C);
			samp_arms(&n, xinit+j, a, b, &A, &B, &C);
			//Rprintf("%4f, %4f, %4f\n", A, Phi[j**m+j], C);
			//Rprintf("%4f %4f\n", pow(xinit[j], *a), pow(xinit[j], *b), log(xinit[j]));
			xa[j] = pow(xinit[j], *a);
			xb[j] = pow(xinit[j], *b);
			//xb[j] = *b ? pow(xinit[j], *b) : log(xinit[j]);
		}
	}
}

inline double sum(int len, double *v){
	// Sums over v;  uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += v[0]; total1 += v[1]; total2 += v[2]; total3 += v[3];
		total4 += v[4]; total5 += v[5]; total6 += v[6]; total7 += v[7];
		i += UNIT; v += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += *(v++);
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

inline double abs_sum(int len, double *v){
	// Sums over |v|;  uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += fabs(v[0]); total1 += fabs(v[1]); total2 += fabs(v[2]); total3 += fabs(v[3]);
		total4 += fabs(v[4]); total5 += fabs(v[5]); total6 += fabs(v[6]); total7 += fabs(v[7]);
		i += UNIT; v += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += fabs(*(v++));
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}


inline double in_order_dot_prod_pow(int len, double *l, double *r, double lpow, double rpow){
	// Computes dot product between l^lpow and r^rpow; uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += pow(l[0], lpow) * pow(r[0], rpow); total1 += pow(l[1], lpow) * pow(r[1], rpow);
		total2 += pow(l[2], lpow) * pow(r[2], rpow); total3 += pow(l[3], lpow) * pow(r[3], rpow);
		total4 += pow(l[4], lpow) * pow(r[4], rpow); total5 += pow(l[5], lpow) * pow(r[5], rpow);
		total6 += pow(l[6], lpow) * pow(r[6], rpow); total7 += pow(l[7], lpow) * pow(r[7], rpow);
		i += UNIT; l += UNIT; r += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += (pow(*(l++), lpow)) * (pow(*(r++), rpow));
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

inline double in_order_tri_dot_prod(int len, double *l, double *m, double *r){
	// Sums over l*m*r; uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += l[0] * m[0] * r[0]; total1 += l[1] * m[1] * r[1];
		total2 += l[2] * m[2] * r[2]; total3 += l[3] * m[3] * r[3];
		total4 += l[4] * m[4] * r[4]; total5 += l[5] * m[5] * r[5];
		total6 += l[6] * m[6] * r[6]; total7 += l[7] * m[7] * r[7];
		i += UNIT; l += UNIT; m += UNIT; r += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += (*(l++)) * (*(m++)) * (*(r++));
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

inline double in_order_tri_dot_prod_pow(int len, double *l, double *m, double *r, double lpow, double mpow, double rpow){
	// Sums over l^lpow*m^mpow*r^rpow; uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += pow(l[0], lpow) * pow(m[0], mpow) * pow(r[0], rpow); total1 += pow(l[1], lpow) * pow(m[1], mpow) * pow(r[1], rpow);
		total2 += pow(l[2], lpow) * pow(m[2], mpow) * pow(r[2], rpow); total3 += pow(l[3], lpow) * pow(m[3], mpow) * pow(r[3], rpow);
		total4 += pow(l[4], lpow) * pow(m[4], mpow) * pow(r[4], rpow); total5 += pow(l[5], lpow) * pow(m[5], mpow) * pow(r[5], rpow);
		total6 += pow(l[6], lpow) * pow(m[6], mpow) * pow(r[6], rpow); total7 += pow(l[7], lpow) * pow(m[7], mpow) * pow(r[7], rpow);
		i += UNIT; l += UNIT; m += UNIT; r += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += (pow(*(l++), lpow)) * (pow(*(m++), mpow)) * (pow(*(r++), rpow));
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

// elts for centered
void elts_c(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier){
	//// d is here although center, because we can reuse its calculation in forming g1
	int n = *nIn, p = *pIn, j,k,l, jn, lp, ln, jpp;
	// g1 <- crossprod(hpx, x)/n; d <- colMeans(hx); diag(g1) <- diag(g1) + d
	for (j=0; j<p; j++){ // Column
		d[j] = sum(n, hx+j*n) / n;
		for (k=0; k<p; k++){ // Column
			g1[j*p+k] = in_order_dot_prod(n, hpx+(k*n), x+j*n) / n;
		}
		g1[j*p+j] += d[j];
	}
	// Gamma is the big matrix flattened (or its upper-left flattened for non-profiled non-centered)
	// Will be modified later for profiled non-centered
	for (j=0; j<p; j++){
		jpp = j*p*p; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=0; k<p; k++){
				Gamma[jpp+lp+k] = in_order_tri_dot_prod(n, x+k*n, x+ln, hx+jn) / n;
			}
			diagonals_with_multiplier[j*p+l] = *diagonal_multiplier * Gamma[jpp+lp+l]; //Gamma[jpp+lp+l] *= *diagonal_multiplier;
		}
	}
}

// elts for non-centered and non-profiled; calls elts_c
void elts_nc_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	elts_c(nIn, pIn, hx, hpx, x, g1, d, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	// Gamma12flat <- -crossprod(x,hx)/n
	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = -in_order_dot_prod(n, x+k*n, hx+j*n) / n;
		}
	}
	// g2 = -colMeans(hpx)
	for (j=0; j<p; j++){ // Column
		g2[j] = -sum(n, hpx+j*n) / n;
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


void elts_ab_c(int *nIn, int *pIn, double *a, double *hx, double *hpx, double *x, double *g1, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier){
	// NOTE THAT b IS NOT NEEDED
	// NOTE THAT UNLIKE elts_c for gaussian, d IS NOT IN THE ARGUMENT LIST because we cannot reuse it in g1
	int n = *nIn, p = *pIn, j, k, l, jn, lp, ln, jpp;
	double A = *a;
	for (j=0; j<p; j++){
		jpp = j*p*p; jn = j*n;
		for (l=0; l<p; l++){
			lp = l*p; ln = l*n;
			for (k=0; k<p; k++){
				Gamma[jpp+lp+k] = 0;
				for (int i=0; i<n; i++)
					Gamma[jpp+lp+k] += pow(x[k*n+i]*x[ln+i], A)*hx[jn+i]*pow(x[jn+i], 2*A-2);
				Gamma[jpp+lp+k] /= n;
			}
			diagonals_with_multiplier[j*p+l] = *diagonal_multiplier * Gamma[jpp+lp+l]; //Gamma[jpp+lp+l] *= *diagonal_multiplier;
			g1[j*p+l] = (in_order_tri_dot_prod_pow(n, hpx+j*n, x+j*n, x+l*n, 1, A-1, A)
						 +(A-1)*in_order_tri_dot_prod_pow(n, hx+j*n, x+j*n, x+l*n, 1, A-2, A)) / n;
		}
		g1[j*p+j] += A*in_order_dot_prod_pow(n, hx+j*n, x+j*n, 1, 2*A-2) / n;
	}
}

// elts for generalized a/b setting; this function is based on elts_nc_np
void elts_ab_np(int *nIn, int *pIn, double *a, double *b, double *hx, double *hpx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	double A = *a, B = *b;
	elts_ab_c(nIn, pIn, a, hx, hpx, x, g1, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	for (j=0; j<p; j++){
		d[j] = in_order_dot_prod_pow(n, hx+j*n, x+j*n, 1, 2*B-2) / n;
		g2[j] -= (in_order_dot_prod_pow(n, hpx+j*n, x+j*n, 1, B-1) + (B-1)*in_order_dot_prod_pow(n, hx+j*n, x+j*n, 1, B-2)) / n;
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = -in_order_tri_dot_prod_pow(n, x+k*n, hx+j*n, x+j*n, A, 1, A+B-2) / n;
		}
	}
}

void elts_exp_c(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier){
	// Note that although we treat the centered case here, we include d here because it's needed for g1; g2 is also included b/c it's convenient
	// Also used for elts_gamma_np since it's the same; g2 and d are just ignored
	int n = *nIn, p = *pIn, j, k, l, jn, lp, ln, jpp;
	for (j=0; j<p; j++){ // Column
		double tmp, tmp2 = 0; d[j] = 0; g2[j] = 0;
		for (k=0; k<p; k++) // Column
			g1[j*p+k] = 0;
		for (int i=0; i<n; i++){
			tmp2 = hx[j*n+i]/x[j*n+i];
			tmp = (hpx[j*n+i]-0.5*tmp2)/sqrt(x[j*n+i]);
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
			for (k=0; k<p; k++){
				Gamma[jpp+lp+k] = 0;
				for (int i=0; i<n; i++)
					Gamma[jpp+lp+k] += sqrt(x[k*n+i]*x[ln+i])*hx[jn+i]/x[jn+i];
				Gamma[jpp+lp+k] /= n;
			}
			diagonals_with_multiplier[j*p+l] = *diagonal_multiplier * Gamma[jpp+lp+l]; //Gamma[jpp+lp+l] *= *diagonal_multiplier;
		}
	}
}

// elts for sqrt exp graphical models; corresponds to a = b = 1/2
void elts_exp_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	int n = *nIn, p = *pIn, j, k;
	elts_exp_c(nIn, pIn, hx, hpx, x, g1, g2, d, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = 0;
			for (int i=0; i<n; i++)
				Gamma12[j*p+k] -= sqrt(x[k*n+i])*hx[j*n+i]/x[j*n+i];
			Gamma12[j*p+k] /= n;
		}
	}
}

// elts for gamma graphical models; corresponds to a = 1/2, b = 0
void elts_gamma_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_exp_c(nIn, pIn, hx, hpx, x, g1, g2, d, Gamma, diagonal_multiplier, diagonals_with_multiplier);
	int n = *nIn, p = *pIn, j, k;
	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			Gamma12[j*p+k] = 0;
			for (int i=0; i<n; i++)
				Gamma12[j*p+k] -= sqrt(x[k*n+i]/x[j*n+i])*hx[j*n+i]/x[j*n+i];
			Gamma12[j*p+k] /= n;
		}
	}
	for (j=0; j<p; j++){ // Column
		g2[j] = 0; d[j] = 0;
		for (int i=0; i<n; i++){
			g2[j] -= (hpx[j*n+i]-hx[j*n+i]/x[j*n+i])/x[j*n+i];
			d[j] += hx[j*n+i]/x[j*n+i]/x[j*n+i];
		}
		g2[j] /= n; d[j] /= n;
	}
}

// elts for non-centered and profiled; calls elts_nc_np and make_profile
void elts_nc_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_nc_np(nIn, pIn, hx, hpx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_ab and make_profile
void elts_ab_p(int *nIn, int *pIn, double *a, double *b, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_ab_np(nIn, pIn, a, b, hx, hpx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_exp and make_profile
void elts_exp_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_exp_np(nIn, pIn, hx, hpx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
	make_profile(pIn, g1, g2, d, Gamma, Gamma12, diagonals_with_multiplier);
}

// elts for non-centered and profiled; calls elts_gamma and make_profile
void elts_gamma_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier){
	elts_gamma_np(nIn, pIn, hx, hpx, x, g1, g2, d, Gamma, Gamma12, diagonal_multiplier, diagonals_with_multiplier);
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
	if (diagonals_with_multiplier != NULL) // If not refit
		crit += lambda2 * abs_sum(p, eta);
	return (crit);
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

double set_KKT_bound(double bound, double tol){
	return ((bound * EPS >= 10*tol) ? (bound * (EPS)) : (bound - 10*tol));
}

void estimator_profiled (int *pIn, double *Gamma_K, double *g_K, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, int *exclude, double *diagonals_with_multiplier, int *gauss){ // gauss: untruncated Gaussian if True
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
	if (oldK ==0){
		Rprintf("Out of Memory!\n");
		return;
	}
	for (i=0; i<p; i++){
		for (j=i; j<p; j++){
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
		
		for (i=0; i<p; i++){
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
							grad = (-in_order_dot_prod(p, Gamma_K+i*p, K+j*p) - in_order_dot_prod(p, Gamma_K+j*p, K+i*p)) / 2;
						else
							grad = (g_K[j*p+i] + g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p)) / 2;
						if (fabs(grad) > KKT_bound){
							need_rerun = 1; exclude[j*p+i] = 0; exclude[i*p+j] = 0;
						}
					}
				}//grad[p*p-1] = g_K[p*p-1] - in_order_dot_prod(p, Gamma_K+p*p*p+p*p, K+p*p);
				if (!first_time && !need_rerun)
					break;
			}
			estimator_profiled(pIn, Gamma_K, g_K, K, lambda1In, tol, maxit, iters, converged, exclude, diagonals_with_multiplier,  gauss);
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
	
	//double pen = 0;
	// Add penalty to loss
	//for (j=0;j<(p-1);j++) // Sum over upper triangle, since we only need l1 norm of off-diagonal entries
	//    pen += abs_sum(p-j-1, K+j*p+j+1);
	//*crit_pen = *crit + pen * 2 * (*lambda1In);
	
	/*
	 for (i=0; i<(p-1); i++){
	 for (j=i+1; j<p; j++){
	 exclude[i*p+j] = (fabs(K[i*p+j]) <= *tol); exclude[j*p+i] = exclude[i*p+j];
	 }
	 }
	 if (*refit){
	 double *K_throw = (double*)malloc(p*p*sizeof(double)), *lambda_throw = (double*)malloc(sizeof(double));
	 int *converged_throw = (int*)malloc(sizeof(int));
	 *converged_throw=0; *lambda_throw=0;
	 for (i=0; i<p; i++){
	 K_throw[i*p+i] = K[i*p+i];
	 for (j=i+1; j<p; j++){
	 K_throw[i*p+j] = K[i*p+j]; K_throw[j*p+i] = K_throw[i*p+j];
	 }
	 }
	 estimator_profiled (pIn, Gamma_K, g_K, K_throw, lambda_throw, tol, maxit, converged_throw, exclude);
	 *crit_refit = loss_profiled(p, Gamma_K, g_K, K_throw);
	 free(K_throw); free(converged_throw); free(lambda_throw);
	 }*/
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
	for (i=0; i<p; i++){
		oldeta[i] = eta[i];
		for (j=i; j<p; j++){
			oldK[lindx(i,j,p)] = K[i*p+j];
		}
	}
	
	//*crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta); ////
	//Rprintf("Crit %f \n", *crit); ////
	
	
	*iters = 0;
	while (*iters < *maxit){
		(*iters)++;
		maxdiff = 0;
		
		// update off-diagonal elements of K
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
		if (gauss)
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
							grad = ((-in_order_dot_prod(p, Gamma_K+i*p, K+j*p) - in_order_dot_prod(p, Gamma_K+j*p, K+i*p) - Gamma_K_eta[j]*eta[i] - Gamma_K_eta[i]*eta[j]))/2;
						else
							grad = ((g_K[j*p+i] + g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) - Gamma_K_eta[i*p+j]*eta[i] - Gamma_K_eta[j*p+i]*eta[j]))/2;
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




/*
 //// Screening on both K and eta: does not always work well
 void full_full(int *pIn, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, int *exclude_eta, double *previous_lambda1, double *previous_lambda2, int *is_refit){
 int p=*pIn;
 double KKT_bound1, KKT_bound2;
 if (*is_refit) // If refit, directly estimate with support restricted to exclude
 estimator_full_penalized (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta);
 else {
 KKT_bound1 = set_KKT_bound(2*(*lambda1In)-(*previous_lambda1), *tol);
 KKT_bound2 = set_KKT_bound(2*(*lambda2In)-(*previous_lambda2), *tol); //
 int first_time = 1, total_iters = 0, i,j;
 while (TRUE){
 int need_rerun = 0;
 // Only calculate for those not currently in the edge set
 for (i = 0; i < (p-1); i++){
 if (exclude_eta[i] && fabs(g_eta[i] - Gamma_eta[i]*eta[i] - in_order_dot_prod(p, K+i*p, Gamma_K_eta+i*p)) > KKT_bound2){
 need_rerun += 1; exclude_eta[i] = 0;
 }
 for (j = i+1; j < p; j++){
 if (!exclude[i*p+j])
 continue;
 double grad = ((g_K[j*p+i] + g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+j*p*p+i*p, K+j*p) - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) - Gamma_K_eta[i*p+j]*eta[i] - Gamma_K_eta[j*p+i]*eta[j]))/2;
 if (fabs(grad) > KKT_bound1){
 need_rerun += 1; exclude[j*p+i] = 0; exclude[i*p+j] = 0;
 }
 }
 }
 //if (need_rerun) {if (first_time) Rprintf("lambda=%f: %d added\n", *lambda1In, need_rerun); else Rprintf("lambda=%f: %d need rerun\n", *lambda1In, need_rerun);} ////
 if (!first_time && !need_rerun)
 break;
 estimator_full_penalized (pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, exclude, exclude_eta);
 //if (first_time) Rprintf("lambda=%f, %d iters\n", *lambda1In, *iters); else Rprintf("lambda=%f, %d MORE iters\n", *lambda1In, *iters);
 total_iters += *iters;
 first_time = 0; KKT_bound1 = set_KKT_bound(*lambda1In, *tol); KKT_bound2 = set_KKT_bound(*lambda2In, *tol);
 }
 estimator_full_penalized(pIn, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta, lambda1In, lambda2In, tol, maxit, iters, converged, NULL, NULL);
 *iters = total_iters + *iters;
 }
 *crit = loss_full_penalized(p, Gamma_K, Gamma_K_eta, Gamma_eta, g_K, g_eta, K, eta);
 
 //double pen = 0;
 // Add penalty to loss
 //for (j=0; j<(p-1); j++) // Sum over upper triangle, since K is symmetric by construction and we want l1 norm of off-diagonal entries
 //    pen += abs_sum(p-j-1, K+j*p+j+1);
 //pen *= 2 * lambda1;
 //pen += lambda2 * abs_sum(p, eta);
 // *crit_pen = *crit + pen;
 }*/



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
	if (oldK ==0){
		Rprintf("Out of Memory!\n");
		return;
	}
	for (i=0; i<p; i++){
		for (j=0; j<p; j++){ // Kji
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
		// update all elements of K, diagonal or off-diagonal
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
							grad = -in_order_dot_prod(p, Gamma_K+j*p, K+i*p);
						else
							grad = g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p);
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

void estimator_full_penalized_asymm (int *pIn, double *Gamma_K, double *Gamma_K_eta, double *Gamma_eta, double *g_K, double *g_eta, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, int *exclude, int *exclude_eta, double *diagonals_with_multiplier,  int *gauss){
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
	for (i=0; i<p; i++){
		oldeta[i] = eta[i];
		for (j=i; j<p; j++){ // Kji
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
							grad = -in_order_dot_prod(p, Gamma_K+j*p, K+i*p) - Gamma_K_eta[j]*eta[i];
						else
							grad = g_K[i*p+j] - in_order_dot_prod(p, Gamma_K+i*p*p+j*p, K+i*p) - Gamma_K_eta[i*p+j]*eta[i];
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




