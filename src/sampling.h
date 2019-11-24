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

int int_runif(int lo, int hi);
double rexp_truncated(double lo, double hi);
double rlaplace_truncated(double lo, double hi);
double rand_init(const int *num_intervals, const double *lefts, const double *rights,  const double *finite_infinity);


void samp_arms(const int *n, const int *every, double *samp, const double *a, const double *b, const double *A, const double *B, const double *C, const int *num_intervals, double *lefts, double *rights, const double *finite_infinity, int *errno);

void rexp_gamma_reject(int *gamm, double *xinit, double *sqrtx, int *steps, int *m, double *Theta, double *Phi, int *max_iter, int *seed);

void rab_arms(const int *nsamp, const int *burnin, const int *m, const int *every, const int *a_numer, const int *a_denom, const int *b_numer, const int *b_denom, double *xinit, double *xres, const double *Theta, const double *Phi, int *seed, const double *finite_infinity, const int *num_char_params, const char **char_params, const int *num_int_params, int *int_params, int *num_double_params, double *double_params, int *errno);
