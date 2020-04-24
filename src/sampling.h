//
//  sampling.h
//  
//
//  Created by Shiqing Yu.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif

struct parm {
	double xl, xr;
	int a_numer, a_denom, b_numer, b_denom;
};

struct ab_parm {
	// Elements for density on R^p
	struct parm base;
	double A, B, C;
	int abs, num_intervals;
	double *fused, *disp, *lefts, *rights;
};

int int_runif(int lo, int hi);
double rexp_truncated(double lo, double hi);
double rlaplace_truncated(double lo, double hi);
double rand_init(const int *num_intervals, const double *lefts, const double *rights);
double random_init_laplace(const int *num_intervals, const double *lefts, const double *rights, double *center);
double laplace_center(struct ab_parm *ab_data);

void samp_arms(const int not_simplex, const int *n, const int *every, double *samp, double (*den)(double x, void *den_elts), void *den_elts);

void rexp_gamma_reject(int *gamm, double *xinit, double *sqrtx, int *steps, int *p, double *eta, double *K, int *max_iter//, int *seed
);

void rab_arms(const int *nsamp, const int *burnin, const int *p, const int *every, const int *a_numer, const int *a_denom, const int *b_numer, const int *b_denom, const int *abs, double *xinit, double *xres, const double *eta, const double *K, //int *seed,
			  double *finite_infinity, const int *num_char_params, const char **char_params, const int *num_int_params, int *int_params, int *num_double_params, double *double_params, int *verbose);
