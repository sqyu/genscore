//
//  genscore.h
//  
//
//  Created by Shiqing Yu on 2017-09-10.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif

void elts_gauss_c(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_gauss_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);

void make_profile(int *pIn, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonals_with_multiplier);

void elts_ab_c(int *nIn, int *pIn, double *a, double *hx, double *hpx, double *x, double *g1, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_ab_np(int *nIn, int *pIn, double *a, double *b, double *hx, double *hpx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);

void elts_exp_c(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_exp_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);

void elts_gamma_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);

void elts_loglog_c(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *d, double *Gamma, double *diagonal_multiplier, double *diagonals_with_multiplier,
				   double *logx, double *h_over_xsq, double *hp_over_x);
void elts_loglog_np(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2,  double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);


void elts_gauss_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_ab_p(int *nIn, int *pIn, double *a, double *b, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_exp_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_gamma_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);
void elts_loglog_p(int *nIn, int *pIn, double *hx, double *hpx, double *x, double *g1, double *g2, double *d, double *Gamma, double *Gamma12, double *diagonal_multiplier, double *diagonals_with_multiplier);


double shrink(double a, double b);

int lindx(int r, int c, int p);

double set_KKT_bound(double bound, double tol);


double loss_profiled (int p, double *Sd20, double *Sd11, double *K, double *diagonals_with_multiplier, double lambda1);
double loss_profiled_gauss (int p, double *Gamma_K, double *K, double *diagonals_with_multiplier, double lambda1);
double loss_full_penalized(int p, double *Sd20, double *Sd21, double *Sd22, double *Sd11, double *Sd12, double *K, double *eta, double *diagonals_with_multiplier, double lambda1, double lambda2);
double loss_full_penalized_gauss (int p, double *Gamma_K, double *Gamma_K_eta, double *K, double *eta, double *diagonals_with_multiplier, double lambda1, double lambda2);

// For non-centered case with penalty on both K and eta = K*mu
void estimator_full_penalized (int *pIn, double *Sd20, double *Sd21, double *Sd22, double *Sd11, double *Sd12, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, int *exclude, int *exclude_eta, double *diagonals_with_multiplier, int *gauss);

// For centered case, or non-centered case with profiled likelihood (i.e. no penalty on eta)
void estimator_profiled (int *pIn, double *Sd20, double *Sd11, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, int *exclude, double *diagonals_with_multiplier, int *gauss);

void profiled(int *pIn, double *Sd20, double *Sd11, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss);

void full(int *pIn, double *Sd20, double *Sd21, double *Sd22, double *Sd11, double *Sd12, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, int *exclude_eta, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss);


void estimator_full_penalized_asymm (int *pIn, double *Sd20, double *Sd21, double *Sd22, double *Sd11, double *Sd12, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, int *exclude, int *exclude_eta, double *diagonals_with_multiplier, int *gauss);

// For centered case, or non-centered case with profiled likelihood (i.e. no penalty on eta)
void estimator_profiled_asymm (int *pIn, double *Sd20, double *Sd11, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, int *exclude, double *diagonals_with_multiplier, int *gauss);

void profiled_asymm(int *pIn, double *Sd20, double *Sd11, double *K, double *lambda1In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss);

void full_asymm(int *pIn, double *Sd20, double *Sd21, double *Sd22, double *Sd11, double *Sd12, double *K, double *eta, double *lambda1In, double *lambda2In, double *tol, int *maxit, int *iters, int *converged, double *crit, int *exclude, int *exclude_eta, double *previous_lambda1, int *is_refit, double *diagonals_with_multiplier, int *gauss);

