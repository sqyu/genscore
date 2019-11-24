/* Compile with R CMD SHLIB arms.c utils.c set_ops.c domain.c sampling.c genscore.c tests.c -o genscore.so
 Written by Shiqing Yu.
 */

#include <ctype.h>
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
#include "set_ops.h"

#define TOL 1e-6
#define MAXOPSTACK 100
#define MAXNUMSTACK 100


void poly_domain_nonnegative_1d(int a, int b, double c, int larger,
								int *num_intervals, double **lefts_pt, double **rights_pt, int *errno){
	// If requires x to be non-negative.
	*num_intervals = 0; // Default
	*lefts_pt = (double*)malloc(sizeof(double)); // Max num of intervals is 1
	*rights_pt = (double*)malloc(sizeof(double));
	if (b == 0) { // Special functions log and exp for b == 0
		if (a == 0) { // a == 0, b == 0 -> log
			*num_intervals = 1;
			if (larger) {
				**lefts_pt = exp(c); **rights_pt = INFINITY;
			} else {
				**lefts_pt = 0; **rights_pt = exp(c);
			}
		} else if (a == 1) { // a == 1, b == 0 -> exp
			if (c <= 1) {
				if (larger) {
					*num_intervals = 1;
					**lefts_pt = 0; **rights_pt = INFINITY;
				} else  // exp(x) < c <= 1 is never true for x >= 0
					return;
			} else {
				*num_intervals = 1;
				if (larger) {
					**lefts_pt = log(c); **rights_pt = INFINITY;
				} else {
					**lefts_pt = 0; **rights_pt = log(c);
				}
			}
		} else { // Undefined
			*errno = 1; *num_intervals = 0;
			Rprintf("!!!x^(%d/0) not defined!!!\n", a);
			return;
		}
	} else if (a == 0) { // If b != 0 and a == 0, x^(a/b) always 1.0
		if ((c >= 1.0 && !larger) || (c <= 1.0 && larger)) {
			*num_intervals = 1;
			**lefts_pt = 0; **rights_pt = INFINITY;
		}
	}
	else if (c <= 0) { // b != 0 and a != 0 and c <= 0
		if (larger) { // x^(a/b) >= 0 >= c is always true for x >= 0; otherwise never true
			*num_intervals = 1;
			**lefts_pt = 0; **rights_pt = INFINITY;
		}
	} else { // b != 0 and a != 0 and c > 0
		*num_intervals = 1;
		if ((a > 0) ^ larger) { // x^(a/b) < c for a > 0 or x^(a/b) > c for a < 0
			**lefts_pt = 0; **rights_pt = frac_pow(c, b, a, false, errno);
		} else { // x^(a/b) > c for a > 0 or x^(a/b) < c for a < 0
			**lefts_pt = frac_pow(c, b, a, false, errno); **rights_pt = INFINITY;
		}
	}
	return;
}

void log_exp_domain_1d(int a, double c, int larger, int abs,
					   int *num_intervals, double **lefts_pt, double **rights_pt, int *errno) {
	if (a == 0) { // a == 0, b == 0 -> log
		if (abs) { // log(|x|)
			if (larger) { // log(|x|) > c
				*num_intervals = 2;
				(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = exp(c);
				(*rights_pt)[0] = -exp(c); (*rights_pt)[1] = INFINITY;
			} else { // log(|x|) < c
				*num_intervals = 1;
				**lefts_pt = -exp(c); **rights_pt = exp(c);
			}
		} else { // log(x)
			*num_intervals = 1;
			if (larger) { // log(x) > c
				**lefts_pt = exp(c); **rights_pt = INFINITY;
			} else { // log(x) < c
				**lefts_pt = 0; **rights_pt = exp(c);
			}
		}
	} else if (a == 1) { // a == 1, b == 0 -> exp
		if (c > 0 && (!abs)) { // exp(x) compared to c > 0
			*num_intervals = 1;
			if (larger) { // exp(x) > c
				**lefts_pt = log(c); **rights_pt = INFINITY;
			} else { // exp(x) < c
				**lefts_pt = -INFINITY; **rights_pt = log(c);
			}
		} else if (c > 1 && abs) { // exp(|x|) compared to c > 1
			if (larger) { // exp(|x|) > c
				*num_intervals = 2;
				(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = log(c);
				(*rights_pt)[0] = -log(c); (*rights_pt)[1] = INFINITY;
			} else { // exp(|x|) < c
				*num_intervals = 1;
				**lefts_pt = -log(c); **rights_pt = log(c);
			}
		} else {
			if (larger) { // exp(x) > 0 >= c or exp(|x|) > 1 >= c always true
				*num_intervals = 1;
				**lefts_pt = -INFINITY; **rights_pt = INFINITY;
			} else // exp(x) < c <= 0 or exp(|x|) < c <= 1 never true
				return;
		}
	} else {
		*errno = 1; *num_intervals = 0;
		Rprintf("!!!x^(%d/0) not defined!!!\n", a);
		return;
	}
}

void poly_domain_1d(int a, int b, double c, int larger, int abs, int nonnegative,
					  int *num_intervals, double **lefts_pt, double **rights_pt, int *errno){
	/* Returns the intervals on which x^(a/b) >/< c (larger == 1 or 0), or abs(x) if abs == 1.
	 Assumes b > 0, and a and b coprime. Ignores Lebesgue null sets and
	 joins intervals whenever possible (i.e. (-INF,0),(0,INF) will be joined to (-INF,INF)).
	 */
	*num_intervals = 0; // default setting; if error occurs returns with 0
	*lefts_pt = (double*)malloc(2 * sizeof(double)); // Max num of intervals is 2, so preset for convenience
	*rights_pt = (double*)malloc(2 * sizeof(double));
	if (b < 0) {
		*errno = 1;
		Rprintf("!!!power denominator must be > 0!! Got %d!!!\n", b);
		return;
	}
	if (nonnegative) {
		poly_domain_nonnegative_1d(a, b, c, larger, num_intervals, lefts_pt, rights_pt, errno);
		return;
	}
	if (b == 0) { // Special functions log and exp for b == 0
		log_exp_domain_1d(a, c, larger, abs, num_intervals, lefts_pt, rights_pt, errno);
		return;
	}
	if (a == 0) { // For b != 0, x^(a/b) always 1.0
		if ((c <= 1.0 && larger) || (c >= 1.0 && !larger)) {
			*num_intervals = 1;
			**lefts_pt = nonnegative ? 0 : -INFINITY; **rights_pt = INFINITY;
		}
		return;
	}
	if (b % 2 == 0) { // x (or |x|) must be >= 0
		if (c <= 0) { // b even, x >= 0, so x^(a/b) >= c always holds true for c <= 0.
			if (larger) { // Always >; never true if <
				*num_intervals = 1;
				**lefts_pt = abs ? -INFINITY : 0; // x: [0,INF); abs(x): (-INF,INF)
				**rights_pt = INFINITY;
			}
		} else { // b even, c > 0
			double cba = frac_pow(c, b, a, false, errno); // c^(b/a)
			if (abs) {
				if ((a > 0) ^ larger) { // |x| < C^(b/a): [-c^(b/a), c^(b/a)]
					*num_intervals = 1;
					**lefts_pt = -cba; **rights_pt = cba;
				} else { // |x| > C^(b/a): (-INF, -c^(b/a)], [c^(b/a), INF)
					*num_intervals = 2;
					(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = cba;
					(*rights_pt)[0] = -cba; (*rights_pt)[1] = INFINITY;
				}
			} else { // x^a | C^b
				*num_intervals = 1;
				if ((a > 0) ^ larger) { // 0 < x < C^(b/a): [0, c^(b/a)]
					**lefts_pt = 0; **rights_pt = cba;
				} else { // x > C^(b/a): [c^(b/a), INF)
					**lefts_pt = cba; **rights_pt = INFINITY;
				}
			}
		} // b even and c > 0
	} else { // b odd
		abs = abs || (a % 2 == 0); // If a even, x^a is equivalent to |x|^a
		if (abs){
			if (c <= 0) { // |x|^(a/b) >= c always holds true for c <= 0
				if (larger) { // (-INF, INF)
					*num_intervals = 1;
					**lefts_pt = -INFINITY; **rights_pt = INFINITY;
				}
			} else { // b odd, abs or a even, c > 0
				double cba = frac_pow(c, b, a, false, errno);
				if ((a > 0) ^ larger) { // [-c^(b/a), c^(b/a)]
					*num_intervals = 1;
					**lefts_pt = -cba; **rights_pt = cba;
				} else { // (-INF, -c^(b/a)], [c^(b/a), INF)
					*num_intervals = 2;
					(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = cba;
					(*rights_pt)[0] = -cba; (*rights_pt)[1] = INFINITY;
				}
			} // b odd, abs or a even, c > 0
		} else { // b odd, no abs and a odd
			if (c == 0) { // b odd, no abs and a odd, c == 0: x^(a/b) | 0 same as x | 0
				if (larger) { // x^(a/b) > 0: [0, INF)
					*num_intervals = 1; **lefts_pt = 0; **rights_pt = INFINITY;
				} else { // x^(a/b) < 0: (-INF, 0]
					*num_intervals = 1; **lefts_pt = -INFINITY; **rights_pt = 0;
				}
			} else { // b odd, no abs and a odd, c != 0
				double cba = frac_pow(c, b, a, false, errno);
				if (a > 0) { // b odd, no abs and a odd, c != 0, a > 0
					*num_intervals = 1;
					if (larger) { // x^(a/b) | c is simply x | c^(b/a): [c^(b/a), INF)
						**lefts_pt = cba; **rights_pt = INFINITY;
					} else { // (-INF, c^(b/a)]
						**lefts_pt = -INFINITY; **rights_pt = cba;
					}
				} else { // b odd, no abs and a odd, c != 0, a < 0
					if (c > 0) { // b odd, no abs and a odd, c > 0, a < 0
						if (larger) { // [0, c^(b/a)]
							*num_intervals = 1;
							**lefts_pt = 0; **rights_pt = cba;
						} else { // (-INF,0], [c^(b/a),+INF)
							*num_intervals = 2;
							(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = cba;
							(*rights_pt)[0] = 0; (*rights_pt)[1] = INFINITY;
						}
					} else { // b odd, no abs and a odd, c < 0, a < 0
						if (larger) { // (-INF,c^(b/a)], [0,+INF)
							*num_intervals = 2;
							(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = 0;
							(*rights_pt)[0] = cba; (*rights_pt)[1] = INFINITY;
						} else { // [c^(b/a), 0]
							*num_intervals = 1;
							**lefts_pt = cba; **rights_pt = 0;
						}
					} // b odd, no abs and a odd, c < 0, a < 0
				}  // b odd, no abs and a odd, c != 0, a < 0
			} // b odd, no abs and a odd, c != 0
		}  // b odd, no abs and a odd
	} // b odd
}

void domain_1d(const int *idx, const int *m, const double *x,
			  const int *num_char_params, const char **char_params,
			  const int *num_int_params, int *int_params,
			  int *num_double_params, double *double_params,
			  int *num_intervals, double **lefts_pt, double **rights_pt,
			  double **cache, int *errno){
	/*
	 cache: For caching results common to all *idx to speed up calculation of h. Boundary type-specific. Defined as double** since each boundary type may need to allocate different size of memory.
	 	uniform: Not used.
	 	polynomial: If not NULL, (*cache) would cache the sum of powers (after being calculated for *idx = 0); used to avoid repeated calculations of powers when calculating boundaries for all dimensions for a single fixed x (for calculating h(x)).
	 
	 char_params[0] must specify the type.
	 uniform: boundaries uniform and independent of values of the other components. double_params must either be empty (entire R), or contain the left end-points followed by the right ones. Thus, *num_double_params must be even.
	 polynomial:
	 	char_params[1] must specify the logic operation on the domains defined by the polynomials, in the postfix notation. E.g. "1 2 | 3 &" returns the region which is the union of those defined by equations 0 and 1, then intersected with that defined by equation 2.
	 	The first four elements of int_params must be [larger or smaller than (1 or 0), abs in x (1 or 0), uniform powers (1 or 0), restrict to non-negative orthant (1 or 0)].
	 	double_params should contain 1 element -- the constant to compare to.
	 	If int_params[0] is 1, the domain will be double_params[0]*x0^pow0 + ... + double_params[m-1]*x(m-1)^pow(m-1) >= double_params[m], otherwise <= will be used.
	 	If int_params[1] is 1, |x0|, ..., |x(m-1)| will be used.
	 	If int_params[2] is 1, pow0 = ... = pow(m-1) = int_params[4] / int_params[5]; otherwise, pow0 = int_params[4] / int_params[4+m], ... , pow(m-1) = int_params[3+m] / int_params[3+2*m].
	 	If int_params[3] is 1, x will be restricted to the non-negative orthant.
	 */
	if (*num_char_params < 1) {
		*errno = 1;
		Rprintf("!!!Number of string parameters must be at least 1 and the first string must be the domain type!!!\n");
		return;
	}
	if (*idx < 0 || *idx >= *m) {
		*errno = 1;
		Rprintf("!!!Invalid index %d. Should be between [0, %d]!!!\n", *idx, *m-1);
		return;
	}
	if (strcmp(char_params[0], "uniform") == 0) {
		if (*num_double_params % 2) {
			*errno = 1;
			Rprintf("!!!For uniform boundaries, num_double_params must be an even number and double params should contain the left endpoints followed by the right endpoints!!!\n");
			return;
		}
		if (*num_double_params == 0) {
			*num_double_params = 2;
			double_params = (double*)malloc(2 * sizeof(double));
			double_params[0] = -INFINITY; double_params[1] = INFINITY;
		}
		*num_intervals = *num_double_params / 2;
		*lefts_pt = double_params; *rights_pt = double_params + *num_intervals;
		return;
	} else if (strcmp(char_params[0], "polynomial") == 0) {
		if (*num_char_params != 2) {
			*errno = 1;
			Rprintf("!!!Number of string parameters must be 2 for polynomial type boundary, where the second string parameter is the rule for combining the boundaries defined by each polynomial!!!\n");
			return;
		}
		int num_eqs = int_params[0];
		int_params++;
		const char *postfix = char_params[1];
		int *num_intervals_list = (int*)malloc(num_eqs * sizeof(int));
		double **lefts_list = (double**)malloc(num_eqs * sizeof(double *));
		double **rights_list = (double**)malloc(num_eqs * sizeof(double *));
		// Calculates the domain for each equation
		for (int eq_i = 0; eq_i < num_eqs; eq_i++) {
			int larger = int_params[0], abs = int_params[1];
			int unif_pow = int_params[2], nonneg = int_params[3];
			int power_numer_this, power_denom_this;
			double sum_other = 0;
			// Calculate total sum of powers, or read from cache if available
			if (unif_pow) { // Uniform power
				reduce_gcd(int_params + 4, int_params + 5); // Make the numer and denom coprime
				if (cache && *idx != 0) // If pow sum already calculated (m != 0) and cached
					sum_other = (*cache)[eq_i];
				else
					for (int j = 0; j < *m; j++)
						sum_other += double_params[j] * frac_pow(x[j], int_params[4], int_params[5], abs, errno);
				power_numer_this = int_params[4]; power_denom_this = int_params[5];
				int_params += 6;
			} else { // Different powers
				if (cache && *idx != 0) // If pow sum already calculated (m != 0) and cached
					sum_other = (*cache)[eq_i];
				else
					for (int j = 0; j < *m; j++) {
						reduce_gcd(int_params + 4 + j, int_params + 4 + *m + j); // Make the numer and denom coprime
						sum_other += double_params[j] * frac_pow(x[j], int_params[4 + j], int_params[4 + *m + j], abs, errno);
					}
				power_numer_this = int_params[*idx + 4]; power_denom_this = int_params[*idx + 4 + *m];
				int_params += 4 + 2 * (*m);
			}
			// Cache if needed
			if (cache && *idx == 0) { // If need to cache pow sum
				if (eq_i == 0) // Allocate memory and initialize
					*cache = (double*)malloc(num_eqs * sizeof(double));
				(*cache)[eq_i] = sum_other;
			}
			// Calculate sum of powers of OTHER variables only
			sum_other -= double_params[*idx] * frac_pow(x[*idx], power_numer_this, power_denom_this, abs, errno);
			// Budget for this variable
			double budget = double_params[*m] - sum_other;
			// If coefficient is 0
			if (fabs(double_params[*idx]) < TOL) {
				if ((budget > -TOL && larger) || (budget < TOL && !larger))  // If want larger than bound but sum_other already < bound OR want smaller than bound but sum_other already > bound -> x always out of bound, bound for *idx is empty set
					num_intervals_list[eq_i] = 0;
				else {
					num_intervals_list[eq_i] = 1;
					lefts_list[eq_i] = (double *)malloc(sizeof(double));
					*(lefts_list[eq_i]) = nonneg ? 0 : -INFINITY;
					rights_list[eq_i] = (double *)malloc(sizeof(double));
					*(rights_list[eq_i]) = INFINITY;
				}
			} else {
				if (double_params[*idx] < 0) {
					larger = !larger;
					budget = -budget / double_params[*idx];
				} else
					budget = budget / double_params[*idx];
				poly_domain_1d(power_numer_this, power_denom_this,
						   budget, larger, abs, nonneg,
						   num_intervals_list + eq_i, lefts_list + eq_i, rights_list + eq_i, errno);
				if (*errno) {
					Rprintf("!!!Error occurred when calling poly_domain_1d() in domain_1d()!!!\n");
					return;
				}
			}
			double_params += *m + 1;
		}
		// Evaluates the domain using domains for each equation with intersections/unions
		evaluate_logic(&num_eqs, postfix, num_intervals_list, lefts_list,
					   rights_list, num_intervals, lefts_pt, rights_pt, errno);
		free(num_intervals_list); free(lefts_list); free(rights_list);
		return;
	}
}


void poly_domain_1d_for_R(int *a, int *b, double *c, int *larger, int *abs, int *nonnegative, int *num_intervals, double *lefts, double *rights, int *print, int *errno){
	// Prints and returns 1d domain for a domain specified by a polynomial x^(a/b) >/< c, used for calls from R. Assumes lefts and rights have been allocated at least 2 elements since the max number of intervals is 2.
	double *lefts0, *rights0;
	poly_domain_1d(*a, *b, *c, *larger, *abs, *nonnegative, num_intervals, &lefts0, &rights0, errno);
	if (*errno) {
		Rprintf("!!!Error occurred when calling poly_domain_1d in poly_domain_1d_for_R()!!!\n");
		return;
	}
	for (int i = 0; i < *num_intervals; i++) {
		if (*print)
			Rprintf("Interval %d: [%f, %f]\n", i, lefts0[i], rights0[i]);
		lefts[i] = lefts0[i]; rights[i] = rights0[i];
	}
}

void h(const int *m, const double *x, double *hx, double *hpx, const double *h_pow, const double *h_trunc,
	   const int *num_char_params, const char **char_params,
	   const int *num_int_params, int *int_params,
	   int *num_double_params, double *double_params, int *errno){
	double *lefts, *rights;
	double **cache = (double**)malloc(sizeof(double*)); // Caching results common to all idx
	int num_intervals;
	if (*h_pow <= 0) {*errno = 1; Rprintf("!!!h_pow must be > 0!!!\n"); return;}
	if (*h_trunc <= 0) {*errno = 1; Rprintf("!!!h_trunc must be > 0!!!\n"); return;}
	for (int idx = 0; idx < *m; idx++) {
		domain_1d(&idx, m, x, num_char_params, char_params,
			   num_int_params, int_params, num_double_params, double_params,
			   &num_intervals, &lefts, &rights, cache, errno);
		if (*errno) {
			Rprintf("!!!Error ocurred when calling domain_1d() in h()!!!\n");
			return;
		}
		int bin = search_unfused(lefts, rights, num_intervals, x[idx], errno);
		if (*errno) {
			Rprintf("!!!Error occurred when calling search_unfused() in h()!!!\n");
			return;
		}
		if (bin == -1) {
			*errno = 1;
			Rprintf("!!!h: x given not in domain!!!\n");
			return;
		}
		hx[idx] = *h_trunc; hpx[idx] = 0;
		if (lefts[bin] != -INFINITY && x[idx] - lefts[bin] < hx[idx]) {
			hx[idx] = x[idx] - lefts[bin];
			hpx[idx] = 1;
		}
		if (rights[bin] != INFINITY && rights[bin] - x[idx] < hx[idx]) {
			hx[idx] = rights[bin] - x[idx];
			hpx[idx] = -1;
		}
		if (hpx[idx] != 0) {
			if (hx[idx] > 0) // If x not at boundary: a(x-x0)^(a-1) or -a(x0-x)^(a-1)
				hpx[idx] *= *h_pow * pow(hx[idx], *h_pow - 1);
			else
				hpx[idx] = 0; // at boundary
		}
		hx[idx] = pow(hx[idx], *h_pow); // |x-x0|^(a-1)
	}
	free(cache);
}

void fuse_endpoints(const int *num_intervals, const double *lefts, const double *rights,
					double *fused, double *disp, int *errno)
/* to fuse intervals with endpoints specified by lefts and rights, left-aligned;
 [x1, x2], [x3, x4], [x5, x6] becomes {x1, x2, x2+x4-x3, x2+x4-x3+x6-x5}
 *num_intervals: number of intervals
 *lefts: left bounds, an array of length num_intervals
 *rights: right bounds, an array of length num_intervals
 *fused: the fused intervals, an array of length num_intervals + 1
 *disp: displacements of intervals after fusion, an array of length num_intervals,
 e.g. {0, x3-x2, x5-x4+x3-x2}
 */
{
	fused[0] = lefts[0];
	fused[1] = rights[0];
	disp[0] = 0.0;
	if (*num_intervals < 1) {
		*errno = 1;
		Rprintf("!!!In fuse_endpoints: number of intervals < 1!!!\n");
	}
	for (int interval = 1; interval < *num_intervals; interval++){
		fused[1 + interval] = fused[interval] + rights[interval] - lefts[interval];
		disp[interval] = disp[interval - 1] + lefts[interval] - rights[interval - 1];
	}
}

/* ********************************************************************* */

int binarySearch_fused(const double *arr, int l, int r, double x, int *errno)
/* Binary search for the left endpoint of the FUSED interval x belongs to.
 Assumes arr has at least 2 elements and x is inside the range of arr.
 *arr: an array
 l: the left index to search (inclusive)
 r: the right index to search (inclusive)
 x: the number to be searched
 */
{
	if (r > l + 1) {
		int mid = (l + r) / 2;
		if (arr[mid] >= x)
			return binarySearch_fused(arr, l, mid, x, errno);
		return binarySearch_fused(arr, mid, r, x, errno);
	}
	return l;
}

int naiveSearch_fused(const double *arr, int length, double x, int *errno)
/* Naive search for the left endpoint of the FUSED interval x belongs to.
 Assumes arr has at least 2 elements and x is inside the range of arr.
 *arr: an array
 length: the length of arr
 x: the number to be searched
 */
{
	for (int i = 1; i <= length - 1; i++)
		if (arr[i] >= x)
			return i - 1;
	*errno = 1;
	Rprintf("!!!search_fused: %f not in fused domain!!!\n", x);
	return -1;
}

int search_fused(const double *arr, int length, double x, int *errno)
/* to return the index of left endpoint of the FUSED interval x belongs to using
 naive search; i.e. find smallest i s.t. arr[i] <= x <= arr[i+1]
 (if x is at the boundaries the smaller index is picked).
 Calls binarySearch_fused() when length > 8 and naiveSearch_fused() otherwise.
 *arr: an array
 length: the length of arr
 x: the number to be searched
 */
{
	if (length < 2 || x < arr[0] || x > arr[length - 1]){
		*errno = 1;
		Rprintf("!!!search_fused: %f not in fused domain!!!\n", x);
		return -1;
	}
	return (length > 8) ? binarySearch_fused(arr, 0, length - 1, x, errno) : naiveSearch_fused(arr, length, x, errno);
}

double translate_unfuse(double x, int num_intervals, const double *fused,
						const double *disp, int *errno)
/* translates x (relative the fused intervals) back to the x in the original domain.
 Example: [0, 1], [3, 6], [10, 15], [21, 28] is fused into {0,1,4,9,16};
 given x = 5, search for the bin it belongs to ([4,9]) and add back 6 (since
 [10,15] was translated to [4,9]), so in the original space it becomes x = 11.
 x: the number relative to the fused intervals
 num_intervals: the number of intervals
 *fused: the fused intervals, an array of length num_intervals + 1
 *disp: displacements of intervals after fusion, an array of length num_intervals
 */
{	if (num_intervals == 1) return x;  // No translation required
	int bucket = search_fused(fused, num_intervals + 1, x, errno);
	if (bucket == -1) {
		*errno = 1;
		Rprintf("!!!Unable to find a bin for %f in translate_unfuse!!!\n", x);
		return fused[0] - 1; // Return an invalid value
	}
	return x + disp[bucket];
}


int search_unfused(const double *lefts, const double *rights, int length, double x, int *errno)
/* to return the index of left endpoint of the ORIGINAL interval x belongs to using
 naive search; i.e. find the unique i s.t. lefts[i] <= x <= rights[i]. Uses naive search.
 *lefts: left bounds, an array of length num_intervals
 *rights: right bounds, an array of length num_intervals
 length: the length of lefts and rights
 x: the number to be searched
 */
{
	if (length < 1 || x < lefts[0] || x > rights[length - 1]) {
		*errno = 1;
		Rprintf("!!!search_unfused: %f not in domain!!!\n", x);
		return -1;
	}
	for (int i = length - 1; i >= 0; i--)
		if (lefts[i] <= x) {
			if (rights[i] >= x)
				return i;
			else {
				*errno = 1;
				Rprintf("!!!search_unfused: %f not in domain!!!\n", x);
				return -1;
			}
		}
	*errno = 1;
	Rprintf("!!!search_unfused: %f not in domain!!!\n", x);
	return -1;
}

double translate_fuse(double x, int num_intervals, const double *lefts,
					  const double *rights, const double *disp, int *errno)
/* translates x in the original domain to its position relative to the fused intervals.
 Example: [0, 1], [3, 6], [10, 15], [21, 28] is fused into {0,1,4,9,16};
 given x = 11, search for the bin it belongs to ([10,15]) and subtract 6 (since
 [10,15] should be translated to [4,9]), so in the fused space it becomes x = 5.
 x: the number in the original domain
 num_intervals: the number of intervals
 *lefts: left bounds, an array of length num_intervals
 *rights: right bounds, an array of length num_intervals
 *disp: displacements of intervals after fusion, an array of length num_intervals
 */
{	if (num_intervals == 1) return x;  // No translation required
	int bucket = search_unfused(lefts, rights, num_intervals, x, errno);
	if (*errno || bucket == -1) {
		Rprintf("!!!Error occurred when calling search_unfused() in translate_fuse()!!!\n");
		return lefts[0] - 1; // Return an invalid value
	}
	return x - disp[bucket];
}
