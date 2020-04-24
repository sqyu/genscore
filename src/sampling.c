/* 
 Created by Shiqing Yu.
 */

#include <assert.h>
#include <math.h>
#include <R.h>
#include <R_ext/BLAS.h>
#include <Rembedded.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include "arms.h"
#include "sampling.h"
#include "utils.h"
#include "domain.h"


/* ********************************************************************* */

double ab_density(double x, void *ab_data)
/* Used for conditional densities on a non-null subset of R^p:
 log density of exp(A*x^a+B*(x^a)^2+C*x^b) where x^0 := log(x).*/
{
	/* cast voided pointer into pointer to struct normix_parm */
	struct ab_parm *d = ab_data;
	if (x < d->base.xl || x > d->base.xr) return -INFINITY;
	x = translate_unfuse(x, d->num_intervals, d->fused, d->disp);
	if (d->abs)
		x = fabs(x);
	double xa = frac_pow(x, d->base.a_numer, d->base.a_denom, false, TRUE);
	double fx = d->A * xa + d->B * pow(xa, 2) + d->C * frac_pow(x, d->base.b_numer, d->base.b_denom, false, TRUE);
	return fx;
}

/* ********************************************************************* */

struct ab_simplex_parm {
	// Elements for density on the simplex
	struct parm base;
	double Aj, Bj, Cj, Ap, Bp, Cp, Djm;
};

double ab_simplex_density(double xj, void *simplex_data)
/* Used for conditional densities on the simplex: log density of
 exp(Aj*xj^a + Bj*(xj^a)^2 + Cj*xj^b + Ap*xp^a + Bp*(xp^a)^2 + Cp*xp^b + Djm*xj^a*xp^a) where x^0 := log(x).
 Here xp := simplex_data -> base.xr - xj.*/
{
	struct ab_simplex_parm *d = simplex_data;
	if (xj <= 0 || xj >= d->base.xr) return -INFINITY;
	double xp = d->base.xr - xj;
	double xja = frac_pow(xj, d->base.a_numer, d->base.a_denom, false, TRUE);
	double fx = d->Aj * xja + d->Bj * pow(xja, 2) + d->Cj * frac_pow(xj, d->base.b_numer, d->base.b_denom, false, TRUE);
	double xpa = frac_pow(xp, d->base.a_numer, d->base.a_denom, false, TRUE);
	fx += d->Ap * xpa + d->Bp * pow(xpa, 2) + d->Cp * frac_pow(xp, d->base.b_numer, d->base.b_denom, false, TRUE);
	fx += d->Djm * xja * xpa;
	return fx;
}

/* ********************************************************************* */

void update_finite_infinity_and_finitify(const int *num_intervals, double *lefts, double *rights, double *finite_infinity) {
	// Updates finite_infinity to max(finite_infinity, 10*abs(lefts[0]), 10*abs(rights[0]), 10*abs(lefts[-1]), 10*abs(rights[-1])) (Ignoring the Infs in this list) and set lefts[0] and rights[-1] to -finite_infinity and finite_infinity, resp, if they are infinite.
	if (*num_intervals > 1)
		*finite_infinity = fmax(fmax(*finite_infinity, 10 * fabs(rights[0])), 10 * fabs(lefts[*num_intervals - 1]));
	if (!isinf(lefts[0]))
		*finite_infinity = fmax(*finite_infinity, 10 * fabs(lefts[0]));
	if (!isinf(rights[*num_intervals - 1]))
		*finite_infinity = fmax(*finite_infinity, 10 * fabs(rights[*num_intervals - 1]));
	else
		rights[*num_intervals - 1] = *finite_infinity;
	if (isinf(lefts[0]))
		lefts[0] = -*finite_infinity;
}


void form_density_elts1(const double *K, const double *eta, const int p, const int j, const double *xa, const int *a_numer, const int *a_denom, const int *b_numer, const int *b_denom, const int *abs, struct ab_parm *ab_data) {
	// Form objects for conditional log density of xj in exp(-(x^a)'K(x^a)/2 + eta'x^b), where x^0 := log(x).
	double A = K[j+j*p]*xa[j] - in_order_dot_prod(p, K+j*p, xa); // (Need to divide by a below) Coefficient on x^a in -1/2a*x^a'Kx^a for a != 0 and -1/2*log(x)'Klog(x) for a == 0
	double B = -K[j*p+j] / 2;  // (Need to divide by a below) Coefficient on x^(2a) in -1/2a*x^a'Kx^a for a != 0 and -1/2*log(x)'Klog(x) for a == 0
	if (B >= 0)
		error("In rab_arms(): K[%d,%d] needs to be positive for non-simplex domains! %4f provided.\n", j+1, j+1, K[j*p+j]);
	if (*a_denom != 0) {
		double a = (double)(*a_numer) / *a_denom;
		A /= a; // 1/a for a != 0 or 1 for a == 0
		B /= a; // 1/(2a) for a != 0 or 1 for a == 0
	}
	double C = eta[j]; // Coefficient on x^b, need to divide by b below
	if (*b_denom != 0)
		C /= (double)(*b_numer) / *b_denom; // Coefficient on x^b; if b = 0 (log(x)), use eta[j]
	
	/* set up structures for density function */
	ab_data -> base.a_numer = *a_numer;
	ab_data -> base.a_denom = *a_denom;
	reduce_gcd(&(ab_data -> base.a_numer), &(ab_data -> base.a_denom));
	ab_data -> base.b_numer = *b_numer;
	ab_data -> base.b_denom = *b_denom;
	reduce_gcd(&(ab_data -> base.b_numer), &(ab_data -> base.b_denom));
	ab_data -> A = A;
	ab_data -> B = B;
	ab_data -> C = C;
	ab_data -> abs = *abs;
}

void form_density_elts_bounds(const int *num_intervals, double *lefts, double *rights, double *finite_infinity, struct ab_parm *ab_data) {
	if (*num_intervals < 1)
		error("In form_density_elts(): number of intervals must be at least 1.\n");
	update_finite_infinity_and_finitify(num_intervals, lefts, rights, finite_infinity);
	double *fused = (double*)malloc((*num_intervals + 1) * sizeof(double));
	double *displacements = (double*)malloc(*num_intervals * sizeof(double));
	fuse_endpoints(num_intervals, lefts, rights, fused, displacements);
	ab_data -> fused = fused;
	ab_data -> disp = displacements;
	ab_data -> num_intervals = *num_intervals;
	ab_data -> lefts = lefts;
	ab_data -> rights = rights;
	ab_data -> base.xl = fused[0];
	ab_data -> base.xr = fused[*num_intervals];
}

/* ********************************************************************* */

void form_simplex_density_elts(const double *K, const double *eta, const int p, const int j, const double *xa, const double xj_add_xp, const int *a_numer, const int *a_denom, const int *b_numer, const int *b_denom, struct ab_simplex_parm *ab_data) {
	// Form objects for conditional log density of xj in exp(-(x^a)'K(x^a)/2 + eta'x^b) on the simplex, where x^0 := log(x) and
	// where it is assumed that x sum up to 1. Only the first p-1 components are free to vary.
	// xj_add_xp is x[j] + x[p-1], which is the boundary for x[j] (equals to 1-x[1]-...-x[j-1]-x[j-1]-...-x[p-2]).
	double Aj = K[j+j*p] * xa[j] + K[p-1+j*p] * xa[p-1] - in_order_dot_prod(p, K+j*p, xa); // (Need to divide by a below) Coefficient on x^a in -1/2a*x^a'Kx^a for a != 0 and -1/2*log(x)'Klog(x) for a == 0
	double Ap = K[j+(p-1)*p] * xa[j] + K[p-1+(p-1)*p] * xa[p-1] - in_order_dot_prod(p, K+(p-1)*p, xa);
	double Bj = -K[j*p+j] / 2, Bp = -K[p*p-1] / 2; // (Need to divide by a below) Coefficient on x^(2a) in -1/2a*x^a'Kx^a for a != 0 and -1/2*log(x)'Klog(x) for a == 0
	if (*a_denom != 0) {
		double a = (double)(*a_numer) / *a_denom;
		Aj /= a; // 1/a for a != 0 or 1 for a == 0
		Ap /= a;
		Bj /= a; // 1/(2a) for a != 0 or 1 for a == 0
		Bp /= a;
	}
	double Cj = eta[j], Cp = eta[p-1]; // Coefficient on x^b, need to divide by b below
	if (*b_denom != 0) {
		double b = (double)(*b_numer) / *b_denom;
		Cj /= b; // Coefficient on x^b; if b = 0 (log(x)), use eta[j]
		Cp /= b;
	}
	double Djm = -K[j*p+p-1];
	
	ab_data -> base.a_numer = *a_numer;
	ab_data -> base.a_denom = *a_denom;
	reduce_gcd(&(ab_data -> base.a_numer), &(ab_data -> base.a_denom));
	ab_data -> base.b_numer = *b_numer;
	ab_data -> base.b_denom = *b_denom;
	reduce_gcd(&(ab_data -> base.b_numer), &(ab_data -> base.b_denom));
	ab_data -> Aj = Aj;
	ab_data -> Bj = Bj;
	ab_data -> Cj = Cj;
	ab_data -> Ap = Ap;
	ab_data -> Bp = Bp;
	ab_data -> Cp = Cp;
	ab_data -> Djm = Djm;
	ab_data -> base.xl = 0;
	ab_data -> base.xr = xj_add_xp;
}

/* ********************************************************************* */

double translate_unfuse_unified(double x, void *density_elts) {
	if (density_elts) {
		struct ab_parm *d = density_elts;
		return translate_unfuse(x, d->num_intervals, d->fused, d->disp);
	} else return x;
}

double translate_fuse_unified(double x, void *density_elts) {
	if (density_elts) {
		struct ab_parm *d = density_elts;
		return translate_fuse(x, d->num_intervals, d->lefts, d->rights, d->disp);
	} else return x;
}

/* ********************************************************************* */

void samp_arms(const int not_simplex, const int *n, const int *every, double *samp, double (*den)(double x, void *den_elts), void *den_elts)
/*
 *n: number of samples to return
 *every: how many samples to draw for generating each sample; only the last one will be kept;
 good practice to set every > 1
 samp: array for storing the samples; must have n doubles allocated
 */
{
	// In Gibbs sampling, this only works if samp is the point from the previous step, not when samp is fixed or completely random
	
	int neval, ninit = 17, npoint = 100, nsamp = *n, ncent = 0, dometrop = 1;
	double xinit[17]={1.0,5.0,10.0,20.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,80.0,90.0,95.0,99.0};
	double *xsamp = (double*)malloc((*every) * sizeof(double));
	double xl = ((struct parm*)den_elts) -> xl, xr = ((struct parm*)den_elts) -> xr;
	double xcent[1], qcent[1], convex = 1.;
	for (int i = 0; i < ninit; i++)
		xinit[i] = xinit[i] / 100 * (xr - xl) + xl;
	
	if (xr - xl < TOL) { // If degenerate domain (singleton)
		double xmid = (xl + xr) / 2;
		(*den)(xmid, den_elts);
		/*if (*(((struct parm*)den_elts) -> errno_status) || *errno_status) {
			Rprintf("!!!Fused left (%f) and right (%f) endpoints are within tolerance (%f), but error occurred when evaluating conditional density at their midpoint (%f)!!!\n", xl, xr, TOL, xmid);
		} else {*/
		xmid = translate_unfuse_unified(xmid, den_elts);
		for (int i = 0; i < nsamp; i++)
			samp[i] = xmid;
		//}
		return;
	}
	
	/* initialise random number generator */
	// Using samp[0] as the previous value
	double xprev = samp[0]; //fabs(rnorm(1, 1));
	if (not_simplex)
		xprev = translate_fuse_unified(xprev, den_elts);

	if (xprev < xl || xprev > xr)
		error("In samp_arms(): translated xprev = %f, but fused domain is [%f, %f].\n", xprev, xl, xr);
	
	//if (xprev < xl) // If xprev provided not in the domain
	//	exit(1);
	
	for (int i = 0; i < nsamp; i++) { // nsamp = 1
		int err = arms(xinit, ninit, &xl, &xr, den, den_elts, &convex,
					   npoint, dometrop, &xprev, xsamp, *every, qcent, xcent, ncent, &neval);
		if (err > 0)
			error("In samp_arms(): error code in ARMS = %d.\n", err);
		if (isnan(xsamp[*every-1]))
			error("In samp_arms(): NaN generated, possibly due to overflow in density (e.g. with densities involving exp(exp(...))).\n");
		if (xsamp[*every-1] < xl || xsamp[*every-1] > xr)
			error("In samp_arms(): %d-th sample out of range [%f, %f] (fused domain). Got %f.\n", i, xl, xr, xsamp[*every-1]);
		samp[i] = xsamp[*every - 1];
		if (not_simplex)
			samp[i] = translate_unfuse_unified(samp[i], den_elts);
		xprev = xsamp[*every - 1];
	}
	free(xsamp);
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


void rexp_gamma_reject(int *gamm, double *xinit, double *sqrtx, int *steps, int *p, double *eta, double *K, int *max_iter//, int *seed
){
	// Runs Gibbs sampling for
	// sqrtx must be the square root of xinit, this is NOT checked
	int n = 1;
	/*if (*seed >= 0){
		srand(*seed);
		*seed += 1;
	}*/
	for (int iter = 0; iter < *steps; iter++){
		for (int j=0; j < *p; j++){
			double A = 2*(K[j+j**p]*sqrtx[j]-in_order_dot_prod(*p, K+j**p, sqrtx));
			if (!*gamm)
				A += eta[j]*2;
			double C = *gamm ? eta[j] : 0;
			double B = -K[j**p+j];
			one(&n, xinit+j, &A, &B, &C, max_iter);
			if (xinit[j] < 0){ // If failed to generate a sample
				xinit[0] = -1; // Signal R and stop generating
				return;
			}
			sqrtx[j] = sqrt(xinit[j]);
		}
	}
}


int int_runif(int lo, int hi) {
	// One integer from {lo, lo+1, ..., hi-1}
	return (int)(runif(lo, hi));
}

double rexp_truncated(double lo, double hi) {
	// Generates X ~ exp(-x) on (lo, hi)
	return (-log(runif(exp(lo - hi), 1)) + lo);
}

double rlaplace_truncated(double lo, double hi) {
	// Generates X ~ exp(-|x|)
	if (lo >= 0)
		return rexp_truncated(lo, hi);
	if (hi <= 0)
		return -rexp_truncated(-hi, -lo);
	if (runif(0, 1) < (1 - exp(lo)) / (2 - exp(lo) - exp(-hi)))
		return -rexp_truncated(0, -lo);
	else
		return rexp_truncated(0, hi);
}

double rand_init(const int *num_intervals, const double *lefts, const double *rights){
	int bin = int_runif(0, *num_intervals);
	if (lefts[bin] == -INFINITY || rights[bin] == INFINITY) // left or right end infinite
		return (rlaplace_truncated(lefts[bin], rights[bin])); // Generates laplace
	else // both ends finite
		return (runif(lefts[bin], rights[bin])); // Generates uniformly
}

double laplace_center(struct ab_parm *ab_data){
	// Lazily find a center for random_init_laplace, preferrably the x that maximizes exp(A*x^a+B*(x^a)^2+C*x^b) with a=a_numer/a_denom and b=b_numer/b_denom and x^(0/0)=log(x), x^(n/0)=exp(n*x) for n != 0
	// If abs is TRUE, return 0.0 (since the density will be centered around 0 anyways)
	// If a_numer / a_denom != b_numer / b_denom, return 0.0 as the maximizer has no closed form solution
	// Otherwise, A*x^a+B*(x^a)^2+C*x^b is maximized at x^a = -(A+C)/2/B. Return such an x if possible, otherwise 0.0.
	if (ab_data -> abs) return 0.0;
	if (ab_data -> base.a_numer == ab_data -> base.b_numer &&
		ab_data -> base.a_denom == ab_data -> base.b_denom) {
		double xa = -(ab_data -> A + ab_data -> C)/2/ab_data -> B, res = 0.0;
		if (ab_data -> base.a_denom != 0){
			if (xa < 0.0 && (ab_data -> base.a_numer % 2 == 0 || ab_data -> base.a_denom % 2 == 0))
				return 0.0;
			res = frac_pow(xa, ab_data -> base.a_denom, ab_data -> base.a_numer, FALSE, FALSE); // Try xa^(a_denom/a_numer), but silent any errors
		} else if (ab_data -> base.a_numer == 0) // xa = log(x)
			res = exp(xa);
		else // xa = exp(a_numer * x)
			res = xa > 0 ? (log(xa) / ab_data -> base.a_numer) : 0;
		if (isinf(res) || isnan(res))
			return 0.0;
		else
			return res;
	} else
		return 0.0;
}

double random_init_laplace(const int *num_intervals, const double *lefts, const double *rights, double *center){
	if (*num_intervals == 1)
		return (rlaplace_truncated(lefts[0] - *center, rights[0] - *center) + *center);
	double *norm_consts = (double*)malloc((*num_intervals + 1) * sizeof(double));
	norm_consts[0] = 0;
	if (*center > rights[*num_intervals - 1]) // In case center is very large
		*center = rights[*num_intervals - 1] + 1;
	else if (*center < lefts[0]) // In case center is very negative
		*center = lefts[0] - 1;
	for (int i = 0; i < *num_intervals; i++) {
		if (lefts[i] < *center) {
			if (*center < rights[i]) // left < center < right: exp(-|x-center|) on [left, right]
				norm_consts[i+1] = 2 - exp(lefts[i] - *center) - exp(*center - rights[i]);
			else // left < right < center: exp(x-center) on [left, right]
				norm_consts[i+1] = exp(rights[i] - *center) - exp(lefts[i] - *center);
		} else // center < left < right: exp(center-x) on [left, right]
			norm_consts[i+1] = exp(*center - lefts[i]) - exp(*center - rights[i]);
		norm_consts[i+1] += norm_consts[i]; // Cumulative sum of normalizing constants
	}
	for (int i = 1; i < *num_intervals + 1; i++)
		norm_consts[i] /= norm_consts[*num_intervals];
	int bin = search_fused(norm_consts, *num_intervals+1, runif(0, 1));
	free(norm_consts);
	return (rlaplace_truncated(lefts[bin] - *center, rights[bin] - *center) + *center);
}

void rab_arms(const int *nsamp, const int *burnin, const int *p, const int *every, const int *a_numer, const int *a_denom, const int *b_numer, const int *b_denom, const int *abs, double *xinit, double *xres, const double *eta, const double *K, //int *seed,
			  double *finite_infinity, const int *num_char_params, const char **char_params, const int *num_int_params, int *int_params, int *num_double_params, double *double_params, int *verbose){
	// Runs Gibbs sampling for exp(-x^a %*% K %*% x^a/(2a) + eta %*% (x^b-1)/b OR eta %*% log(x))
	int n = 1;//, every = 1;
	
	if (*finite_infinity <= TOL || isinf(*finite_infinity))
		error("In rab_arms(): finite_infinity must be finite and > TOL=%f! Got %f.\n", TOL, *finite_infinity);

	if ((((*a_numer >= 0) ^ (*a_denom >= 0)) && *a_denom != 0)  ||
		(((*b_numer >= 0) ^ (*b_denom >= 0)) && *b_denom != 0))
		error("In rab_arms(): if the denominators are non-zero, a (a_numer/a_denom) and b (b_numer/b_denom) must both be positive.\n");
	/*if (*seed >= 0) {
		srand(*seed);
		*seed += 1;
	}*/

	double *lefts, *rights, *xa = (double*)malloc(*p * sizeof(double));
	for (int j = 0; j < *p; j++)
		xa[j] = frac_pow(xinit[j], *a_numer, *a_denom, *abs, TRUE);
	
	int num_intervals, counter = 0, total_iters = *burnin + *nsamp, progress_pointer = 0;
	double *print_checkpoints;
	if (*verbose) print_progress_setup(&print_checkpoints, total_iters);

	if (strcmp(char_params[0], "simplex") == 0) { // Simplex
		if (fabs(sum(*p, xinit) - 1) > TOL)
			error("In rab_arms(): sum(xinit) must be close to 1 for simplex.\n");
		for (int iter = 0; iter < total_iters; iter++){
			for (int j = 0; j < *p-1; j++){
				double xj_add_xp = xinit[j] + xinit[*p-1];
				struct ab_simplex_parm *ab_simplex_data = (struct ab_simplex_parm*)malloc(sizeof(struct ab_simplex_parm));
				form_simplex_density_elts(K, eta, *p, j, xa, xj_add_xp, a_numer, a_denom, b_numer, b_denom, ab_simplex_data); // These elements depend on all x except for xj and xp-1
				
				// Randomly generates an initial point for x[j] since the domain might have changed
				xinit[j] = runif(0, xj_add_xp);
				samp_arms(FALSE, &n, every, xinit+j, ab_simplex_density, ab_simplex_data);

				xinit[*p-1] = xj_add_xp - xinit[j];
				// Update x^a
				xa[j] = frac_pow(xinit[j], *a_numer, *a_denom, *abs, TRUE);
				xa[*p-1] = frac_pow(xinit[*p-1], *a_numer, *a_denom, *abs, TRUE);
				if (iter >= *burnin)
					xres[counter++] = xinit[j];
				free(ab_simplex_data);
			}
			if (iter >= *burnin)
				xres[counter++] = xinit[*p-1];
			if (*verbose)
				print_progress(print_checkpoints, &progress_pointer, iter, total_iters);
		}
	} else { // Non-simplex
		for (int iter = 0; iter < total_iters; iter++){
			for (int j = 0; j < *p; j++){
				// Calculate domain for x[j]
				domain_1d(&j, p, xinit, num_char_params, char_params,
						 num_int_params, int_params, num_double_params, double_params,
						 &num_intervals, &lefts, &rights, NULL);
				struct ab_parm *ab_data = (struct ab_parm*)malloc(sizeof(struct ab_parm));
				// Construct data for calculating densities
				form_density_elts1(K, eta, *p, j, xa, a_numer, a_denom, b_numer, b_denom, abs, ab_data);
				// Randomly generates an initial point for x[j] since the domain might have changed
				//xinit[j] = rand_init(&num_intervals, lefts, rights); // Allows INFINITY in lefts and rights
				double center = laplace_center(ab_data);
				xinit[j] = random_init_laplace(&num_intervals, lefts, rights, &center);
				if (*finite_infinity < 10 * fabs(xinit[j]))
					*finite_infinity = 10 * fabs(xinit[j]);
				// Fuse the multiple intervals in the domain into one and construct data for calculating densities; note that INFINITY is changed to FINITE_INFINITY in this call
				form_density_elts_bounds(&num_intervals, lefts, rights, finite_infinity, ab_data);
				samp_arms(TRUE, &n, every, xinit+j, ab_density, ab_data);
				// Update x^a
				xa[j] = frac_pow(xinit[j], *a_numer, *a_denom, *abs, TRUE);
				if (iter >= *burnin)
					xres[counter++] = xinit[j];
				free(ab_data->fused);
				free(ab_data->disp);
				free(ab_data);
				free(lefts); free(rights);
			}
			if (*verbose)
				print_progress(print_checkpoints, &progress_pointer, iter, total_iters);
		}
	}
	free(xa);
	if (*verbose) free(print_checkpoints);
}
