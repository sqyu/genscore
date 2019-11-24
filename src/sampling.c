/* Compile with R CMD SHLIB arms.c utils.c domain.c sampling.c genscore.c -o genscore.so
 Created by Shiqing Yu.
 */

#include <assert.h>
//#include <errno.h>
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

struct ab_parm {
	double a, b, c, A, B, C, xl, xr;
	int num_intervals;
	double *fused, *disp;
	int *errno;
};

/* ********************************************************************* */

double ab_density(double x, void *ab_data)
/* log density of exp(A*x^a+B*x^(2a)+C*x^b) (or C*log(b) if b == 0)*/
{
	struct ab_parm *d;
	double y;
	/* cast voided pointer into pointer to struct normix_parm */
	d = ab_data;
	
	if (x < d->xl || x > d->xr) return -INFINITY;
	x = translate_unfuse(x, d->num_intervals, d->fused, d->disp, d->errno);
	
	if (*(d->errno))
		return NAN;
	
	y = d->A*pow(x,d->a)+d->B*pow(x,2*d->a);
	y += d->b == 0.0 ? d->C*log(x) : d->C*pow(x,d->b);
	return y;
};

/* ********************************************************************* */

void samp_arms(const int *n, const int *every, double *samp, const double *a, const double *b, const double *A, const double *B, const double *C, const int *num_intervals, double *lefts, double *rights, const double *finite_infinity, int *errno)
/*
 *n: number of samples to return
 *every: how many samples to draw for generating each sample; only the last one will be kept;
 good practice to set every > 1
 samp: array for storing the samples; must have n doubles allocated
 *a, *b, *A, *B, *C: parameters for the density exp(A*x^a+B*x^(2a)+C*x^b) (or C*log(b) if b == 0)
 *num_intervals: number of subintervals in the domain
 *lefts: left endpoints of the subintervals; must be FINITE
 *rights: right endpoints of the subintervals; must be FINITE
 lefts and rights must satisfy lefts[i] < rights[i], lefts[i] < lefts[j], rights[i] < rights[j] for all i < j.
 */
{
	// In Gibbs sampling, this only works if samp is the point from the previous step, not when samp is fixed or completely random
	assert(*num_intervals >= 1);
	int neval, ninit = 17, npoint = 100, nsamp = *n, ncent = 0 ;
	double xinit[17]={1.0,5.0,10.0,20.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,80.0,90.0,95.0,99.0};
	
	lefts[0] = fmax(lefts[0], -*finite_infinity);
	rights[*num_intervals - 1] = fmin(rights[*num_intervals - 1], *finite_infinity);
	
	double *fused = (double*)malloc((*num_intervals + 1) * sizeof(double));
	double *displacements = (double*)malloc(*num_intervals * sizeof(double));
	fuse_endpoints(num_intervals, lefts, rights, fused, displacements, errno);
	
	if (*errno) {
		Rprintf("!!!Error occurred in fuse_endpoints.!!!\n");
		return;
	}
	
	double *xsamp = (double*)malloc((*every) * sizeof(double));
	double xl = fused[0], xr = fused[*num_intervals];
	
	for (int i = 0; i < ninit; i++)
		xinit[i] = xinit[i] / 100 * (xr - xl) + xl;
	double xcent[1], qcent[1];
	double convex = 1.;
	int dometrop = 1;
	
	/* initialise random number generator */
	// Using samp[0] as the previous value
	double xprev = translate_fuse(samp[0], *num_intervals, lefts, rights, displacements, errno);//fabs(rnorm(1, 1));
	
	if (*errno) return;
	
	//if (xprev < xl) // If xprev provided not in the domain
	//	exit(1);
	
	/* set up structures for each density function */
	
	if (*B >= 0){
		*errno = 1;
		Rprintf("!!!B needs to be negative! %4f provided!!!\n", *B);
		return;
	}
	
	struct ab_parm ab_data;
	
	/* initialise data needed by normal mixture density function */
	ab_data.a = *a;
	ab_data.b = *b;
	ab_data.A = *A;
	ab_data.B = *B;
	ab_data.C = *C;
	ab_data.fused = fused;
	ab_data.disp = displacements;
	ab_data.xl = xl;
	ab_data.xr = xr;
	ab_data.num_intervals = *num_intervals;
	ab_data.errno = errno;
	
	for (int i = 0; i < nsamp; i++) {
		int err = arms(xinit,ninit,&xl,&xr,ab_density,&ab_data,&convex,
					   npoint,dometrop,&xprev,xsamp,*every,qcent,xcent,ncent,&neval);
		if (err > 0) {
			Rprintf("!!!error code in ARMS = %d\n", err);
			*errno = 1;
			break;
		}
		if (xsamp[*every-1] < xl || xsamp[*every-1] > xr){
			Rprintf("!!!%d-th sample out of range [%f, %f] (fused domain). Got %f!!!\n", i, xl, xr, xsamp[*every-1]);
			*errno  = 1;
			break;
		}
		samp[i] = translate_unfuse(xsamp[*every-1], *num_intervals, fused, displacements, errno);
		if (samp[i] < fused[0]) {
			Rprintf("!!!Generated sample out of bound!!!\n");
			*errno = 1;
			break;
		}
		xprev = xsamp[*every-1];
	}
	free(fused);
	free(displacements);
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

double rand_init(const int *num_intervals, const double *lefts, const double *rights, const double *finite_infinity){
	int bin = int_runif(0, *num_intervals);
	if (lefts[bin] == -INFINITY || rights[bin] == INFINITY) // left or right end infinite
		return (rlaplace_truncated(fmax(lefts[bin], -*finite_infinity),
								   fmin(rights[bin], *finite_infinity))); // Generates laplace
	else // both ends finite
		return (runif(lefts[bin], rights[bin])); // Generates uniformly
}

void rand_init_test(int *num_intervals, double *lefts, double *rights, double *res) {
	double finite_infinity = 100;
	*res = rand_init(num_intervals, lefts, rights, &finite_infinity);
}

void rab_arms(const int *nsamp, const int *burnin, const int *m, const int *every, const int *a_numer, const int *a_denom, const int *b_numer, const int *b_denom, double *xinit, double *xres, const double *Theta, const double *Phi, int *seed, const double *finite_infinity, const int *num_char_params, const char **char_params, const int *num_int_params, int *int_params, int *num_double_params, double *double_params, int *errno){
	// Runs Gibbs sampling for exp(x^a%*%Phi%*%x^a/(2a) + Theta %*% (x^b-1)/b OR Theta %*% log(x))
	int n = 1;//, every = 1;
	if (*finite_infinity <= 0) {
		*errno = 1;
		Rprintf("!!!finite_infinity must be > 0! Got %f!!!\n", *finite_infinity);
		return;
	}
	if (*seed >= 0){
		srand(*seed);
		*seed += 1;
	}
	double a = (double)(*a_numer) / *a_denom, b = (double)(*b_numer) / *b_denom;
	int num_intervals, counter = 0;
	double *lefts, *rights, *xa = (double*)malloc((*m) * sizeof(double));
	for (int j = 0; j < *m; j++) {
		xa[j] = frac_pow(xinit[j], *a_numer, *a_denom, false, errno);
		if (*errno) {
			Rprintf("!!!Error occurred when computing power for xinit[%d]!!!\n", j);
			return;
		}
	}

	for (int iter = 0; iter < *burnin + *nsamp; iter++){
		for (int j = 0; j < *m; j++){
			double A = (in_order_dot_prod(*m, Phi+j**m, xa)-Phi[j+j**m]*xa[j]) / a;
			double B = Phi[j**m+j] / 2 / a;
			double C = *b_numer ? Theta[j] / b : Theta[j]; // if b = 0, use Theta[j]
			//Rprintf("A=%4f, B=%4f, C=%4f\n", A, B, C);
			domain_1d(&j, m, xinit, num_char_params, char_params,
					 num_int_params, int_params, num_double_params, double_params,
					 &num_intervals, &lefts, &rights, NULL, errno);
			if (*errno) {
				Rprintf("!!!Error occurred when calculating domain!!!\n");
				return;
			}
			// Randomly generates an initial point since the domain might have changed
			xinit[j] = rand_init(&num_intervals, lefts, rights, finite_infinity);
			samp_arms(&n, every, xinit+j, &a, &b, &A, &B, &C, &num_intervals, lefts, rights, finite_infinity, errno); // Note that in samp_arms INFINITY is changed to FINITE_INFINITY
			if (*errno) {
				Rprintf("!!!Error occurred when calculating domain!!!\n");
				return;
			}
			xa[j] = frac_pow(xinit[j], *a_numer, *a_denom, false, errno);
			if (*errno) {
				Rprintf("!!!Error occurred when calculating domain!!!\n");
				return;
			}
			if (iter >= *burnin)
				xres[counter++] = xinit[j];
		}
	}
}
