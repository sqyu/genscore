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
double frac_pow(double num, int power_numer, int power_denom, int abs, int* errno);

double in_order_sum_uniform_pow(int len, const double *arr, int power_numer, int power_denom, int abs, int* errno);
double in_order_sum_different_pow(int len, const double *arr, int *power_numers, int *power_denoms, int abs, int* errno);
