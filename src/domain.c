/*
 Written by Shiqing Yu.
 */

#include <ctype.h>
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

#define MAXOPSTACK 100
#define MAXNUMSTACK 100


void poly_domain_nonnegative_1d(int a, int b, double c, int larger,
								int *num_intervals, double **lefts_pt, double **rights_pt){
	// If requires x to be non-negative.
	*num_intervals = 0; // Default
	if (b == 0) { // Special functions log and exp for b == 0
		if (a == 0) { // a == 0, b == 0 -> log
			*num_intervals = 1;
			if (larger) {
				**lefts_pt = exp(c); **rights_pt = INFINITY;
			} else {
				**lefts_pt = 0; **rights_pt = exp(c);
			}
		} else if (a > 0) { // a > 0, b == 0 -> exp(a*x) >/< c -> x >/< log(c)/a
			if (c <= 1) { // log(c) < 0
				if (larger) {
					*num_intervals = 1;
					**lefts_pt = 0; **rights_pt = INFINITY;
				} else  // exp(a * x) < c <= 1 is never true for x >= 0 and a > 0
					return;
			} else {
				*num_intervals = 1;
				if (larger) {
					**lefts_pt = log(c) / a; **rights_pt = INFINITY;
				} else {
					**lefts_pt = 0; **rights_pt = log(c) / a;
				}
			}
		} else { // a < 0, b == 0 -> exp(a*x) >/< c -> x </> log(c)/a
			if ((c <= 0 && larger) || (c >= 1 && !larger)) { //
				// exp(a*x) > 0 >= c and exp(a*x) < 1 <= c always true for a < 0 and x >= 0
				*num_intervals = 1; **lefts_pt = 0; **rights_pt = INFINITY;
			} else if ((c <= 0 && !larger) || (c >= 1 && larger))
				// exp(a*x) < c <= 0 and exp(a*x) > c >= 1 never true for a < 0 and x >= 0
				return;
			else {
				*num_intervals = 1;
				if (larger) {
					**lefts_pt = 0; **rights_pt = log(c) / a;
				} else {
					**lefts_pt = log(c) / a; **rights_pt = INFINITY;
				}
			}
		}
	} else if (a == 0) { // If b != 0 and a == 0, x^(a/b) always 1.0
		if ((c >= 1.0 && !larger) || (c <= 1.0 && larger)) {
			*num_intervals = 1;
			**lefts_pt = 0; **rights_pt = INFINITY;
		}
	} else if (c <= 0) { // b != 0 and a != 0 and c <= 0
		if (larger) { // x^(a/b) >= 0 >= c is always true for x >= 0; otherwise never true
			*num_intervals = 1;
			**lefts_pt = 0; **rights_pt = INFINITY;
		}
	} else { // b != 0 and a != 0 and c > 0
		*num_intervals = 1;
		if ((a > 0) ^ larger) { // x^(a/b) < c for a > 0 or x^(a/b) > c for a < 0
			**lefts_pt = 0; **rights_pt = frac_pow(c, b, a, false, TRUE);
		} else { // x^(a/b) > c for a > 0 or x^(a/b) < c for a < 0
			**lefts_pt = frac_pow(c, b, a, false, TRUE); **rights_pt = INFINITY;
		}
	}
	return;
}

void log_exp_domain_1d(int a, double c, int larger, int abs,
					   int *num_intervals, double **lefts_pt, double **rights_pt) {
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
	} else if (a > 0) { // a > 0, b == 0 -> exp(a*x) >/< c -> x >/< log(c) / a
		if (c > 0 && (!abs)) { // x compared to log(c) / a
			*num_intervals = 1;
			if (larger) { // exp(x) > c
				**lefts_pt = log(c) / a; **rights_pt = INFINITY;
			} else { // exp(x) < c
				**lefts_pt = -INFINITY; **rights_pt = log(c) / a;
			}
		} else if (c > 1 && abs) { // |x| compared to log(c) / a > 0
			if (larger) { // exp(|x|) > c
				*num_intervals = 2;
				(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = log(c) / a;
				(*rights_pt)[0] = -log(c) / a; (*rights_pt)[1] = INFINITY;
			} else { // exp(|x|) < c
				*num_intervals = 1;
				**lefts_pt = -log(c) / a; **rights_pt = log(c) / a;
			}
		} else { // c <= 0 OR c <= 1 && abs
			if (larger) { // exp(a*x) > 0 >= c or exp(a*|x|) > 1 >= c always true
				*num_intervals = 1;
				**lefts_pt = -INFINITY; **rights_pt = INFINITY;
			} else // exp(a*x) < c <= 0 or exp(a*|x|) < c <= 1 never true
				return;
		}
	} else { // a < 0, b == 0 -> exp(a*x) >/< c -> x </> log(c) / a
		if (c > 0 && (!abs)) { // x compared to log(c) / a (>/< reversed)
			*num_intervals = 1;
			if (larger) { // exp(a*x) > c -> x < log(c) / a
				**lefts_pt = -INFINITY; **rights_pt = log(c) / a;
			} else { // exp(a*x) < c -> x > log(c) / a
				**lefts_pt = log(c) / a; **rights_pt = INFINITY;
			}
		} else if (c > 0 && c < 1 && abs) { // |x| compared to log(c) / a
			if (larger) { // exp(a*|x|) > c -> |x| < log(c) / a
				*num_intervals = 1;
				**lefts_pt = -log(c) / a; **rights_pt = log(c) / a;
			} else { // exp(a*|x|) < c -> |x| > log(c) / a
				*num_intervals = 2;
				(*lefts_pt)[0] = -INFINITY; (*lefts_pt)[1] = log(c) / a;
				(*rights_pt)[0] = -log(c) / a; (*rights_pt)[1] = INFINITY;
			}
		} else if (c <= 0) {
			if (larger) { // exp(a*x) > 0 >= c always true
				*num_intervals = 1;
				**lefts_pt = -INFINITY; **rights_pt = INFINITY;
			} else // exp(a*x) < c <= 0 never true
				return;
		} else { // c >= 1 && abs
			if (larger) // exp(a|x|) > c >= 1 never true for a < 0
				return;
			else { // exp(a|x|) < c <= 1 always true for a < 0
				*num_intervals = 1;
				**lefts_pt = -INFINITY; **rights_pt = INFINITY;
			}
		}
	}
}

void poly_domain_1d(int a, int b, double c, int larger, int abs, int nonnegative,
					  int *num_intervals, double **lefts_pt, double **rights_pt){
	/* Returns the intervals on which x^(a/b) >/< c (larger == 1 or 0), or abs(x) if abs == 1.
	 Assumes b > 0, and a and b coprime. Ignores Lebesgue null sets and
	 joins intervals whenever possible (i.e. (-INF,0),(0,INF) will be joined to (-INF,INF)).
	 */
	*num_intervals = 0; // default setting; if error occurs returns with 0
	*lefts_pt = (double*)malloc(2 * sizeof(double)); // Max num of intervals is 2, so preset for convenience
	*rights_pt = (double*)malloc(2 * sizeof(double));
	if (b < 0)
		error("In poly_domain_1d(): Power denominator must be > 0. Got %d.\n", b);
	if (nonnegative) {
		poly_domain_nonnegative_1d(a, b, c, larger, num_intervals, lefts_pt, rights_pt);
		return;
	}
	if (b == 0) { // Special functions log and exp for b == 0
		log_exp_domain_1d(a, c, larger, abs, num_intervals, lefts_pt, rights_pt);
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
			double cba = frac_pow(c, b, a, false, TRUE); // c^(b/a)
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
				double cba = frac_pow(c, b, a, false, TRUE);
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
				double cba = frac_pow(c, b, a, false, TRUE);
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

void domain_1d(const int *idx, const int *p, const double *x,
			  const int *num_char_params, const char **char_params,
			  const int *num_int_params, int *int_params,
			  const int *num_double_params, double *double_params,
			  int *num_intervals, double **lefts_pt, double **rights_pt,
			  double **cache){
	/*
	 cache: For caching results common to all *idx to speed up calculation of h. Boundary type-specific. Defined as double** since each boundary type may need to allocate different size of memory.
	 	uniform: Not used.
	 	polynomial: If not NULL, (*cache) would cache the sum of powers (after being calculated for *idx = 0); used to avoid repeated calculations of powers when calculating boundaries for all dimensions for a single fixed x (for calculating h(x)).
	 
	 char_params[0] must specify the type.
	 R: entire R
	 uniform: boundaries uniform and independent of values of the other components. double_params must contain the left end-points followed by the right ones. Thus, *num_double_params must be even. int_params must contain two 0/1s, indicating whether lefts[0] is -Inf and whether rights[-1] is Inf.
	 polynomial:
	 	char_params[1] must specify the logic operation on the domains defined by the polynomials, in the postfix notation. E.g. "1 2 | 3 &" returns the region which is the union of those defined by inequalities 0 and 1, then intersected with that defined by inequality 2.
	 	The first four elements of int_params must be [larger or smaller than (1 or 0), abs in x (1 or 0), uniform powers (1 or 0), restrict to non-negative orthant (1 or 0)].
	 	double_params should contain 1 element -- the constant to compare to.
	 	If int_params[0] is 1, the domain will be double_params[0]*x0^pow0 + ... + double_params[p-1]*x(p-1)^pow(p-1) >= double_params[p], otherwise <= will be used.
	 	If int_params[1] is 1, |x0|, ..., |x(p-1)| will be used.
	 	If int_params[2] is 1, pow0 = ... = pow(p-1) = int_params[4] / int_params[5]; otherwise, pow0 = int_params[4] / int_params[4+p], ... , pow(p-1) = int_params[3+p] / int_params[3+2*p].
	 	If int_params[3] is 1, x will be restricted to the non-negative orthant.
	 */
	if (*num_char_params < 1)
		error("In domain_1d(): Number of string parameters must be at least 1 and the first string must be the domain type.\n");
	if (*idx < 0 || *idx >= *p)
		error("In domain_1d(): Invalid index %d. Should be between [0, %d].\n", *idx, *p-1);
	if (strcmp(char_params[0], "R") == 0 || strcmp(char_params[0], "R+") == 0) {
		*num_intervals = 1;
		*lefts_pt = (double*)malloc(sizeof(double));
		*rights_pt = (double*)malloc(sizeof(double));
		**lefts_pt = strcmp(char_params[0], "R") == 0 ? -INFINITY : 0;
		**rights_pt = INFINITY;
		return;
	} else if (strcmp(char_params[0], "uniform") == 0) {
		if (*num_double_params % 2 || *num_double_params <= 0 || *num_int_params != 2)
			error("In domain_1d(): For uniform boundaries, num_double_params must be an even positive number and double params should contain the left endpoints followed by the right endpoints, and num_int_params must be 2 with int_params[0] indicating whether lefts[0] is -Inf and int_params[1] whether rights[-1] is Inf.\n");
		*num_intervals = *num_double_params / 2;
		*lefts_pt = (double*)malloc(*num_intervals * sizeof(double));
		*rights_pt = (double*)malloc(*num_intervals * sizeof(double));
		for (int j = 0; j < *num_intervals; j++) {
			(*lefts_pt)[j] = double_params[j];
			(*rights_pt)[j] = double_params[j + *num_intervals];
		}
		// *lefts_pt = double_params; *rights_pt = double_params + *num_intervals;
		if (int_params[0])
			(*lefts_pt)[0] = -INFINITY;
		if (int_params[1])
			(*rights_pt)[*num_intervals - 1] = INFINITY;
		return;
	} else if (strcmp(char_params[0], "simplex") == 0) {
		double sum_all = 0;
		if (cache && *idx != 0) // If sum of x already calculated (p != 0) and cached
			sum_all = **cache;
		else {
			for (int j = 0; j < *p; j++) {
				if (x[j] <= -TOL) //// Or <= TOL?
					error("In domain_1d(): x[%d] = %f <= 0 not allowed for the simplex.\n", j, x[j]);
				sum_all += x[j];
			}
			if (sum_all >= 1) //// Or >= 1 - TOL?
				error("In domain_1d(): sum(x) = %f >= 1 not allowed for the simplex.\n", sum_all);
		}
		// Cache if needed
		if (cache && *idx == 0) { // If need to cache sum of x
			*cache = (double*)malloc(sizeof(double));
			**cache = sum_all;
		}
		*num_intervals = 1;
		*lefts_pt = (double*)malloc(sizeof(double));
		*rights_pt = (double*)malloc(sizeof(double));
		**lefts_pt = 0;
		**rights_pt = 1 - (sum_all - x[*idx]);
		return;
	} else if (strcmp(char_params[0], "polynomial") == 0) {
		if (*num_char_params != 2)
			error("In domain_1d(): Number of string parameters must be 2 for polynomial type boundary, where the second string parameter is the rule for combining the boundaries defined by each polynomial.\n");
		int num_eqs = *int_params++;
		const char *postfix = char_params[1];
		int *num_intervals_list = (int*)malloc(num_eqs * sizeof(int));
		double **lefts_list = (double**)malloc(num_eqs * sizeof(double *));
		double **rights_list = (double**)malloc(num_eqs * sizeof(double *));
		
		/*for (int i = 0; i < *p; i++) ////
			Rprintf("x%d=%f, ", i, x[i]); ////
		Rprintf("\n"); ////*/
		
		// Calculates the domain for each inequality
		for (int eq_i = 0; eq_i < num_eqs; eq_i++) {
			int uniform = int_params[0], larger = int_params[1], abs = int_params[2];
			int unif_pow = int_params[3], nonneg = int_params[4];
			
			if (uniform) { // Uniform inequalities
				int power_numer = int_params[5], power_denom = int_params[6];
				double budget = *double_params++;
				int_params += 7;
				poly_domain_1d(power_numer, power_denom, budget, larger, abs, nonneg,
							   num_intervals_list + eq_i, lefts_list + eq_i, rights_list + eq_i);
				/*if (*errno_status) {
					error("In domain_1d(): Error occurred when calling poly_domain_1d() for (uniform) inequality %d. Trying to find domain for %sx%d%s^(%d/%d)%s%f (nonneg=%d).\n", eq_i + 1, abs?"|":"", *idx + 1, abs?"|":"", power_numer, power_denom, larger?">":"<", budget, nonneg);
					return;
				}*/
				for (int i = 0; i < *p; i++)
					if (i != *idx) {
						if (nonneg && x[i] < 0) { // If other components do not satisfy nonnegative requirement, domain for this component should also be empty
							num_intervals_list[eq_i] = 0;
							free(lefts_list[eq_i]); free(rights_list[eq_i]);
							break;
						}
						double tmp = frac_pow(x[i], power_numer, power_denom, abs, FALSE);
						if (isnan(tmp) || (tmp < budget && larger) || (tmp > budget && !larger)) {
							num_intervals_list[eq_i] = 0;
							free(lefts_list[eq_i]); free(rights_list[eq_i]);
							break;
						}
					}
			} else { // Non-uniform inequalities
				int power_numer_this, power_denom_this;
				double sum_all = 0, sum_other = 0;
				
				// If need to cache pow sum
				if (cache && *idx == 0 && eq_i == 0) {
					// Allocate memory and initialize
					*cache = (double*)malloc(num_eqs * sizeof(double));
					for (int i = 0; i < num_eqs; i++)
						(*cache)[i] = NAN;
				}
				
				if (nonneg) {
					int other_components_oob = 0;
					for (int i = 0; i < *p; i++)
						if (i != *idx && x[i] < 0) { // If other components do not satisfy nonnegative requirement, domain for this component should also be empty
							other_components_oob = 1;
							num_intervals_list[eq_i] = 0;
							int_params += unif_pow ? 7 : (5 + 2 * (*p));
							double_params += *p + 1;
							break;
						}
					if (other_components_oob)
						continue;
				}
				
				// Calculate total sum of powers, or read from cache if available
				if (unif_pow) { // Uniform power
					reduce_gcd(int_params + 5, int_params + 6); // Make the numer and denom coprime
					power_numer_this = int_params[5]; power_denom_this = int_params[6];
					if (cache && *idx != 0 && !isnan((*cache)[eq_i])) {// If pow sum already calculated (p != 0) and cached
						sum_all = (*cache)[eq_i];
						sum_other = sum_all - double_params[*idx] * frac_pow(x[*idx], power_numer_this, power_denom_this, abs, FALSE);
					} else {
						for (int j = 0; j < *p; j++)
							if (j != *idx && double_params[j] != 0.0)
								// Calculate sum_other and then sum_all so that even when x[*idx] is out of bound so that frac_pow returns NAN, we can still calculate its domain as long as frac_pow is valid for all other components; calculating sum_all and then sum_other will make both NAN in that case
								sum_other += double_params[j] * frac_pow(x[j], int_params[5], int_params[6], abs, FALSE); // FALSE: do not produce warning if power not well-defined
						sum_all = sum_other + double_params[*idx] * frac_pow(x[*idx], power_numer_this, power_denom_this, abs, FALSE);
					}
					int_params += 7;
				} else { // Different powers
					reduce_gcd(int_params + 5 + *idx, int_params + 5 + *p + *idx);
					power_numer_this = int_params[*idx + 5]; power_denom_this = int_params[*idx + 5 + *p];
					if (cache && *idx != 0 && !isnan((*cache)[eq_i])) {// If pow sum already calculated (p != 0) and cached
						sum_all = (*cache)[eq_i];
						sum_other = sum_all - double_params[*idx] * frac_pow(x[*idx], power_numer_this, power_denom_this, abs, FALSE);
					}
					else {
						for (int j = 0; j < *p; j++)
							if (j != *idx && double_params[j] != 0.0){
								reduce_gcd(int_params + 5 + j, int_params + 5 + *p + j); // Make the numer and denom coprime
								// Calculate sum_other and then sum_all so that even when x[*idx] is out of bound so that frac_pow returns NAN, we can still calculate its domain as long as frac_pow is valid for all other components; calculating sum_all and then sum_other will make both NAN in that case
								sum_other += double_params[j] * frac_pow(x[j], int_params[5 + j], int_params[5 + *p + j], abs, FALSE); // FALSE: do not produce warning if power not well-defined
							}
						sum_all = sum_other + double_params[*idx] * frac_pow(x[*idx], power_numer_this, power_denom_this, abs, FALSE);
					}
					int_params += 5 + 2 * (*p);
				}
				if (cache && *idx == 0) { // If need to cache pow sum
					(*cache)[eq_i] = sum_all; // Cache even if NaN otherwise other *idx might not know the current x is out of domain for inequality eq_i
				}
				if (isnan(sum_other)) { // It is possible that error occurs in frac_pow for other some x[j] at the current x; this x may be generated from other inequalities with an OR rule, and becomes out of bound for the current inequality; simply set the domain for the current inequality to empty set. (E.g. x = [-1,-1,-1] generated from "x^(1/2)>1 || x < 0". Then the domain for x^(1/2)>1 for all components will be set to empty at x = [-1,-1,-1].)
					num_intervals_list[eq_i] = 0;
					double_params += *p + 1;
					continue;
				}
				// Rprintf("eq_i=%d, double_params[*p]=%f, sum_all=%f, sum_other=%f, larger=%d\n", eq_i, double_params[*p], sum_all, sum_other, larger);
				// Budget for this variable
				double budget = double_params[*p] - sum_other;
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
							   num_intervals_list + eq_i, lefts_list + eq_i, rights_list + eq_i);
					/*if (*errno_status) {
						Rprintf("!!!In domain_1d(): Error occurred when calling poly_domain_1d() for (non-uniform) inequality %d. Trying to find domain for %sx%d%s^(%d/%d)%s%f (nonneg=%d)!!!\n", eq_i + 1, abs?"|":"", *idx + 1, abs?"|":"", power_numer_this, power_denom_this, larger?">":"<", budget, nonneg);
						return;
					}*/
				}
				double_params += *p + 1;
			}
		}
		/*for (int j = 0; j < *p; j++)
			Rprintf("idx=%d: x[%d] = %f\n", *idx, j, x[j]);
		for (int eq_i = 0; eq_i < num_eqs; eq_i++) {
			Rprintf("Pre dm(): num_intervals_list[%d]: %d\n", eq_i, num_intervals_list[eq_i]);
			if (num_intervals_list[eq_i])
				Rprintf("Pre dm(): lefts_list[%d]: %p  rights_list[%d]: %p\n", eq_i, (lefts_list[eq_i]), eq_i, (rights_list[eq_i]));
		}*/
		
		// Evaluates the domain using domains for each inequality with intersections/unions
		evaluate_logic(&num_eqs, postfix, num_intervals_list, lefts_list,
					   rights_list, num_intervals, lefts_pt, rights_pt);
		
		/*for (int i = 0; i < num_eqs; i++) {
			Rprintf("\nSubinterval for x%d in eq %d: ", *idx, i);
			for (int j = 0; j < num_intervals_list[i]; j++)
				Rprintf("%d, [%f, %f]; ", j, lefts_list[i][j], rights_list[i][j]);
		}*/
		
		for (int eq_i = 0; eq_i < num_eqs; eq_i++)
			if (num_intervals_list[eq_i]) {
				//Rprintf("dm(): lefts_list[%d]: %p  rights_list[%d]: %p\n", eq_i, (lefts_list[eq_i]), eq_i, (rights_list[eq_i])); ////
				free(lefts_list[eq_i]); free(rights_list[eq_i]);
			}
		free(num_intervals_list); free(lefts_list); free(rights_list);
		return;
	} else
		error("In domain_1d(): domain type %s not supported.\n", char_params[0]);
}


void poly_domain_1d_for_R(int *a, int *b, double *c, int *larger, int *abs, int *nonnegative, int *num_intervals, double *lefts, double *rights, int *print){
	// Prints and returns 1d domain for a domain specified by a polynomial x^(a/b) >/< c, used for calls from R. Assumes lefts and rights have been allocated at least 2 elements since the max number of intervals is 2.
	double *lefts0, *rights0;
	poly_domain_1d(*a, *b, *c, *larger, *abs, *nonnegative, num_intervals, &lefts0, &rights0);
	for (int i = 0; i < *num_intervals; i++) {
		if (*print)
			Rprintf("In poly_domain_1d_for_R(): Interval %d: [%f, %f]\n", i, lefts0[i], rights0[i]);
		lefts[i] = lefts0[i]; rights[i] = rights0[i];
	}
	free(lefts0); free(rights0);
}

void dist(const int *n, const int *p, const double *x, double *dists, int *dist_ps,
	   const int *num_char_params, const char **char_params,
	   const int *num_int_params, int *int_params,
	   int *num_double_params, double *double_params){
	/*
	 x: of size p (num parameters) * n (num samples)
	 dists: of size p * n
	 dist_ps: of size p * n
	 */
	double *lefts, *rights;
	int num_intervals;
	for (int samp_i = 0; samp_i < *n; samp_i++) {
		double **cache = (double**)malloc(sizeof(double*)); // Caches an intermediate result common to all idx for one x (e.g. sum of all components raised to some powers)
		for (int idx = 0; idx < *p; idx++) {
			domain_1d(&idx, p, x, num_char_params, char_params,
					  num_int_params, int_params, num_double_params, double_params,
					  &num_intervals, &lefts, &rights, cache);

			int bin = search_unfused(lefts, rights, num_intervals, x[idx]);
			dist_ps[idx] = 0; // If domain is R, will return 0 derivative
			dists[idx] = INFINITY; // If domain is R, will return INF distance
			if (lefts[bin] != -INFINITY && x[idx] - lefts[bin] < dists[idx]) {
				dists[idx] = x[idx] - lefts[bin];
				dist_ps[idx] = 1;
			}
			if (rights[bin] != INFINITY) {
				double new_dist_minus_old = rights[bin] - x[idx] - dists[idx];
				if (new_dist_minus_old < 0) { // Update distance if shorter to the right
					dists[idx] = rights[bin] - x[idx];
					dist_ps[idx] = -1;
				}
				if (fabs(new_dist_minus_old) < TOL) // If equal distance to left and right
					dist_ps[idx] = 0;
			}
			if (dists[idx] < TOL) // If at boundary, derivative 0
				dist_ps[idx] = 0;
			free(lefts); free(rights);
		}
		if (strcmp(char_params[0], "simplex") == 0 || strcmp(char_params[0], "polynomial") == 0)
			free(*cache);
		free(cache);
		x += *p; dists += *p; dist_ps += *p;
	}
}


void fuse_endpoints(const int *num_intervals, const double *lefts, const double *rights,
					double *fused, double *disp)
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
	if (*num_intervals < 1)
		error("In fuse_endpoints: number of intervals < 1.\n");
	for (int interval = 1; interval < *num_intervals; interval++){
		fused[1 + interval] = fused[interval] + rights[interval] - lefts[interval];
		disp[interval] = disp[interval - 1] + lefts[interval] - rights[interval - 1];
	}
}

/* ********************************************************************* */

int binarySearch_fused(const double *arr, int l, int r, double x)
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
			return binarySearch_fused(arr, l, mid, x);
		return binarySearch_fused(arr, mid, r, x);
	}
	return l;
}

int naiveSearch_fused(const double *arr, int length, double x)
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
	error("In search_fused(): %f not in fused domain.\n", x);
	return -1;
}

int search_fused(const double *arr, int length, double x)
/* to return the index of left endpoint of the FUSED interval x belongs to using
 naive search; i.e. find smallest i s.t. arr[i] <= x <= arr[i+1]
 (if x is at the boundaries the smaller index is picked).
 Calls binarySearch_fused() when length > 8 and naiveSearch_fused() otherwise.
 *arr: an array
 length: the length of arr
 x: the number to be searched
 */
{
	if (length < 2 || x < arr[0] || x > arr[length - 1])
		error("In search_fused(): %f not in fused domain.\n", x);
	return (length > 8) ? binarySearch_fused(arr, 0, length - 1, x) : naiveSearch_fused(arr, length, x);
}

double translate_unfuse(double x, int num_intervals, const double *fused,
						const double *disp)
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
	int bucket = search_fused(fused, num_intervals + 1, x);
	return x + disp[bucket];
}

int search_unfused(const double *lefts, const double *rights, int length, double x)
/* to return the index of left endpoint of the ORIGINAL interval x belongs to using
 naive search; i.e. find the unique i s.t. lefts[i] <= x <= rights[i]. Uses naive search.
 *lefts: left bounds, an array of length num_intervals
 *rights: right bounds, an array of length num_intervals
 length: the length of lefts and rights
 x: the number to be searched
 */
{
	if (length < 1 || x < lefts[0] || x > rights[length - 1])
		error("In search_unfused(): %f not in domain. lefts[0] = %f, rights[%d] = %f. If %f is the finite_infinity you set in make_domain(), please increase the value.\n", x, lefts[0], length, rights[length-1], rights[length-1]);
	for (int i = length - 1; i >= 0; i--)
		if (lefts[i] <= x) {
			if (rights[i] >= x)
				return i;
			else
				error("In search_unfused(): %f not in domain.\n", x);
		}
	error("In search_unfused(): %f not in domain.\n", x);
	return -1;
}

double translate_fuse(double x, int num_intervals, const double *lefts,
					  const double *rights, const double *disp)
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
	int bucket = search_unfused(lefts, rights, num_intervals, x);
	return x - disp[bucket];
}
