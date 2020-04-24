//
//  utils.h
//
//
//  Created by Shiqing Yu.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif

double sum(int len, const double *v);
double abs_sum(int len, const double *v);
double in_order_dot_prod(int len, const double *l, const double *r);

double in_order_tri_dot_prod(int len, const double *l, const double *m, double *r);
double in_order_dot_prod_pow(int len, const double *l, const double *r, double lpow, double rpow);
double in_order_tri_dot_prod_pow(int len, const double *l, const double *m, const double *r, double lpow, double mpow, double rpow);

int gcd(int a, int b);
void reduce_gcd(int *a, int *b);
double frac_pow(double num, int power_numer, int power_denom, int abs, int print_error);

double in_order_sum_uniform_pow(int len, const double *arr, int power_numer, int power_denom, int abs);
double in_order_sum_different_pow(int len, const double *arr, int *power_numers, int *power_denoms, int abs);
double dot_prod_by_row(int len, const double *m, const double *v);


void eliminate_vec(const int *p, double *vec, const int j);
void eliminate_col(const int *n, const int *p, double *mat, const int j);
void eliminate_row(const int *n, const int *p, double *mat, const int j);
void eliminate_row_col(const int *n, const int *p, double *mat, const int j, const int k);


void print_progress_setup(double **checkpoints, int total_iters);
void print_progress(double *checkpoints, int *pointer, int iter, int total_iters);
