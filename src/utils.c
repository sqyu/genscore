#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <sys/param.h>
#include <math.h>
#include <stdio.h>
#include <R_ext/BLAS.h>
#include "utils.h"

#ifdef NAN
/* NAN is supported */
#endif

#define UNIT 8

inline double sum(int len, const double *v){
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

inline double abs_sum(int len, const double *v){
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

inline double in_order_dot_prod(int len, const double *l, const double *r){
	// Computes dot product between vectors of l and r using unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
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

inline double dot_prod_by_row(int len, const double *m, const double *v){
	// Computes dot product between a row of square matrix m and a vector v
	// Specifically, computes m[0]*v[0] + m[len]*v[1] + ... + m[len*(len-1)]*v[len-1]
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += m[0] * v[0]; total1 += m[len] * v[1];
		total2 += m[2*len] * v[2]; total3 += m[3*len] * v[3];
		total4 += m[4*len] * v[4]; total5 += m[5*len] * v[5];
		total6 += m[6*len] * v[6]; total7 += m[7*len] * v[7];
		i += UNIT; m += UNIT * len; v += UNIT;
	}
	for (; i < len; i++){ // Leftovers
		total7 += (*m) * (*(v++));
		m += len;
	}
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

/* ********************************************************************* */
inline int gcd(int a, int b) {
	return b ? gcd(b, a % b) : a;
}

void x(double *c, double *a, double *b){
	*c = INFINITY;
}

inline void reduce_gcd(int *a, int *b){
	if (*b == 0) return; // If log or exp, do nothing
	int c = gcd(*a, *b); // c returned may be negative
	*a = *a / c; *b = *b / c;
	if (*b < 0) {
		*a = -*a; *b = -*b; // To ensure b is > 0.
	}
}

inline double frac_pow(double num, int power_numer, int power_denom, int abs, int print_error) {
/*
 Returns |num| ^ (power_numer / power_denom) if abs == 1 or else num ^ (power_numer / power_denom).
 When one of power_numer and power_denom is 0:
 1. x^(0/0) is understood as log(x).
 2. x^(n/0) is treated as exp(n*x) for n != 0.
 3. x^(0/n) returns 1 for any x and n != 0. Note that 0^0 = 1.
 Otherwise, if num is negative, power_denom must be odd, and num ^ (power_numer / power_denom) = (-1) ^ power_numer * (-num)^(power_numer / power_denom).
 If num is 0 and power_numer / power_denom is negative, NAN will be returned.
 
 print_error should be set to FALSE when frac_pow is used in calculating the sum of other coordinates in calculating the boundary, since
 the current x may be out of domain for the current equation (but should at least be in domain for one equation) and it is thus not guaranteed
 that frac_pow is well-defined.
 */
	if (power_denom == 0) {
		if (power_numer == 0) {
			if (abs) return log(fabs(num));
			if (num > 0) return (log(num));
			if (print_error)
				error("In frac_pow(): x^(0/0) is treated as log(x), but x = %f < 0 is provided.\n", num);
			return NAN;
		} else
			return abs ? exp(power_numer * fabs(num)) : exp(power_numer * num);
	}
	if (power_numer == 0)
		return 1;
	double power = (double)(power_numer) / power_denom;
	if (num > 0) // If num > 0
		return (pow(num, power));
	else if (num == 0) {
		if (power <= 0) {
			if (print_error)
				error("In frac_pow(): 0^(%d/%d) encountered.\n", power_numer, power_denom);
			return NAN;
		}
		return 0;
	}
	if (abs) // If use abs(num)
		return (pow(-num, power));
	if (power_denom % 2 == 0) { // If num < 0 but need to take an even root, error
		if (print_error)
			error("In frac_pow(): A negative number (%f) cannot be raised to a power with even denominator (%d/%d). Returning NAN.\n", num, power_numer, power_denom);
		return NAN;
	}
	if (power_numer % 2) // If power numerator is odd, returned value should be negative
		return (-pow(-num, power));
	else
		return (pow(-num, power));
}

inline double in_order_sum_uniform_pow(int len, const double *arr, int power_numer, int power_denom, int abs){
	// Computes dot product between l^lpow and r^rpow; uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += frac_pow(arr[0], power_numer, power_denom, abs, TRUE);
		total1 += frac_pow(arr[1], power_numer, power_denom, abs, TRUE);
		total2 += frac_pow(arr[2], power_numer, power_denom, abs, TRUE);
		total3 += frac_pow(arr[3], power_numer, power_denom, abs, TRUE);
		total4 += frac_pow(arr[4], power_numer, power_denom, abs, TRUE);
		total5 += frac_pow(arr[5], power_numer, power_denom, abs, TRUE);
		total6 += frac_pow(arr[6], power_numer, power_denom, abs, TRUE);
		total7 += frac_pow(arr[7], power_numer, power_denom, abs, TRUE);
		i += UNIT; arr += UNIT;
	}
	for (; i < len; i++) // Leftovers
		total7 += frac_pow(*(arr++), power_numer, power_denom, abs, TRUE);
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

inline double in_order_sum_different_pow(int len, const double *arr, int *power_numers, int *power_denoms, int abs){
	// Computes dot product between l^lpow and r^rpow; uses unrolled loops
	double total0 = 0, total1 = 0, total2 = 0, total3 = 0, total4 = 0, total5 = 0, total6 = 0, total7 = 0;
	int i = 0;
	while (i < len - len % UNIT) {
		total0 += frac_pow(arr[0], power_numers[0], power_denoms[0], abs, TRUE);
		total1 += frac_pow(arr[1], power_numers[1], power_denoms[1], abs, TRUE);
		total2 += frac_pow(arr[2], power_numers[2], power_denoms[2], abs, TRUE);
		total3 += frac_pow(arr[3], power_numers[3], power_denoms[3], abs, TRUE);
		total4 += frac_pow(arr[4], power_numers[4], power_denoms[4], abs, TRUE);
		total5 += frac_pow(arr[5], power_numers[5], power_denoms[5], abs, TRUE);
		total6 += frac_pow(arr[6], power_numers[6], power_denoms[6], abs, TRUE);
		total7 += frac_pow(arr[7], power_numers[7], power_denoms[7], abs, TRUE);
		i += UNIT; arr += UNIT;
	}
	for (; i < len; i++) // Leftovers
		total7 += frac_pow(*(arr++), *(power_numers++), *(power_denoms++), abs, TRUE);
	return (total0 + total1 + total2 + total3 + total4 + total5 + total6 + total7);
}

/* ********************************************************************* */

inline double in_order_tri_dot_prod(int len, const double *l, const double *m, double *r){
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

inline double in_order_dot_prod_pow(int len, const double *l, const double *r, double lpow, double rpow){
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

inline double in_order_tri_dot_prod_pow(int len, const double *l, const double *m, const double *r, double lpow, double mpow, double rpow){
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


/* ********************************************************************* */
inline void eliminate_vec(const int *p, double *vec, const int j) {
	// Subtracts the j-th component from all components
	double tmp = vec[j];
	for (int k=0; k<*p; k++)
		vec[k] -= tmp;
}

inline void eliminate_col(const int *n, const int *p, double *mat, const int j) {
	// Subtracts the j-th column from all columns
	for (int i=0; i<*n; i++) {
		double tmp = mat[i+(*n)*j];
		for (int k=0; k<*p; k++) {
			mat[i+k*(*n)] -= tmp;
		}
	}
}

inline void eliminate_row(const int *n, const int *p, double *mat, const int j) {
	// Subtracts the j-th row from all rows
	for (int k=0; k<*p; k++)
		eliminate_vec(n, mat+k*(*n), j);
}

inline void eliminate_row_col(const int *n, const int *p, double *mat, const int j, const int k) {
	// Subtracts the j-th row from all rows and then the k-th column from all columns
	eliminate_col(n, p, mat, k);
	eliminate_row(n, p, mat, j);
}

/* ********************************************************************* */

void print_progress_setup(double **checkpoints, int total_iters) {
	double print_percentages[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	*checkpoints = (double*)malloc(14 * sizeof(double));
	int pointer = -1;
	for (int i = 0; i < 13; i++) {
		int print_iter = floor(print_percentages[i] * total_iters);
		if (print_iter >= 1 && (i == 0 || print_iter > (*checkpoints)[pointer]))
			(*checkpoints)[++pointer] = print_iter;
	}
	(*checkpoints)[pointer + 1] = INFINITY; // So that checkpoints[pointer] will never go out of bound when accessed in printing
}

void print_progress(double *checkpoints, int *pointer, int iter, int total_iters) {
	if (iter + 1 >= checkpoints[*pointer]) {
		Rprintf("Sampling %.1f%% done.\n", (iter + 1.0) / total_iters * 100);
		(*pointer)++;
	}
}
