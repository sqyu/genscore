#include <ctype.h>
#include <math.h>
#include <R.h>
#include <R_ext/BLAS.h>
#include <Rembedded.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include "arms.h"
#include "sampling.h"
#include "utils.h"
#include "domain.h"
#include "set_ops.h"

void shunting_yard_test(int *num_eqs, const char **infix) {
	char *postfix = (char*)malloc((eq_nums_len(num_eqs) + 2*(*num_eqs)) * sizeof(char));
	int postfix_length;
	shunting_yard(num_eqs, infix, &postfix, &postfix_length);
	Rprintf("Length %d: %s\n", postfix_length, postfix);
	free(postfix);
}


void logic_domain_test(int *num_eqs, const char **infix, int *num_intervals_list,
					   double *lefts, double *rights) {
	double **lefts_list = (double**)malloc(*num_eqs * sizeof(double*));
	double **rights_list = (double**)malloc(*num_eqs * sizeof(double*));
	for (int i = 0; i < *num_eqs; i++) {
		Rprintf("Equation %d:\n", i + 1);
		Rprintf("\t%d intervals: ", num_intervals_list[i]);
		lefts_list[i] = lefts;
		rights_list[i] = rights;
		for (int j = 0; j < num_intervals_list[i]; j++)
			Rprintf("[%f, %f] ", lefts_list[i][j], rights_list[i][j]);
		lefts += num_intervals_list[i];
		rights += num_intervals_list[i];
		Rprintf("\n");
	}
	char *postfix = (char*)malloc((eq_nums_len(num_eqs) + 2*(*num_eqs)) * sizeof(char));
	int postfix_length;
	shunting_yard(num_eqs, infix, &postfix, &postfix_length);
	Rprintf("In-fix notation: %s\n", *infix);
	Rprintf("Post-fix notation: %s, length %d.\n", postfix, postfix_length);
	int res_num_intervals;
	double *res_lefts, *res_rights;
	evaluate_logic(num_eqs, postfix, num_intervals_list, lefts_list, rights_list,
				   &res_num_intervals, &res_lefts, &res_rights);
	Rprintf("\n Resulting %d intervals:\n\t", res_num_intervals);
	for (int j = 0; j < res_num_intervals; j++)
		Rprintf("[%f, %f] ", res_lefts[j], res_rights[j]);
	Rprintf("\n");
	free(lefts_list); free(rights_list);
}

void frac_pow_test(double *num, int *power_numer, int *power_denom, int *abs) {
	double res = frac_pow(*num, *power_numer, *power_denom, *abs, TRUE);
	Rprintf("Res = %f.\n", res);
}

void intersection_test(const int *A_num_intervals, const double *A_lefts, const double *A_rights, const int *B_num_intervals, const double *B_lefts, const double *B_rights){
	int res_num_intervals = -1;
	double *res_lefts, *res_rights;
	intersection(A_num_intervals, A_lefts, A_rights, B_num_intervals, B_lefts, B_rights, &res_num_intervals, &res_lefts, &res_rights);
	Rprintf("Num intervals: %d\n", res_num_intervals);
	for (int i = 0; i < res_num_intervals; i++)
		Rprintf("Interval %d: [%f, %f]\n", i, res_lefts[i], res_rights[i]);
}

void merge_sorted_test(const int *A_length, double *A, const int *B_length, double *B){
	double *res;
	merge_sorted_arrays(A_length, A, B_length, B, &res);
	for (int i = 0; i < *A_length + *B_length; i++)
		Rprintf("%f, ", res[i]);
	Rprintf("\n");
}

void setunion_test(const int *A_num_intervals, double *A_lefts, double *A_rights, const int *B_num_intervals, double *B_lefts, double *B_rights){
	int res_num_intervals = -1;
	double *res_lefts, *res_rights;
	setunion(A_num_intervals, A_lefts, A_rights, B_num_intervals, B_lefts, B_rights, &res_num_intervals, &res_lefts, &res_rights);
	Rprintf("Num intervals: %d\n", res_num_intervals);
	for (int i = 0; i < res_num_intervals; i++)
		Rprintf("Interval %d: [%f, %f]\n", i, res_lefts[i], res_rights[i]);
}

void rand_init_test(int *num_intervals, double *lefts, double *rights, int *left_inf, int *right_inf, double *res) {
	if (*left_inf) lefts[0] = -INFINITY;
	if (*right_inf) rights[*num_intervals - 1] = INFINITY;
	*res = rand_init(num_intervals, lefts, rights);
}

void laplace_center_test(double *A, double *B, double *C, int *a_numer, int *a_denom, int *b_numer, int *b_denom, int *abs, double *res){
	struct ab_parm *ab_data = (struct ab_parm*)malloc(sizeof(struct ab_parm));
	ab_data -> A = *A;
	ab_data -> B = *B;
	ab_data -> C = *C;
	ab_data -> abs = *abs;
	ab_data -> base.a_numer = *a_numer;
	ab_data -> base.a_denom = *a_denom;
	reduce_gcd(&(ab_data -> base.a_numer), &(ab_data -> base.a_denom));
	ab_data -> base.b_numer = *b_numer;
	ab_data -> base.b_denom = *b_denom;
	reduce_gcd(&(ab_data -> base.b_numer), &(ab_data -> base.b_denom));
	*res = laplace_center(ab_data);
	free(ab_data);
}


void random_init_laplace_test(int *num_intervals, double *lefts, double *rights, int *left_inf, int *right_inf, double *center, double *res) {
	if (*left_inf) lefts[0] = -INFINITY;
	if (*right_inf) rights[*num_intervals - 1] = INFINITY;
	*res = random_init_laplace(num_intervals, lefts, rights, center);
}

void domain_1d_for_R_test(int *idx, int *m, double *x,
						  const int *num_char_params, const char **char_params,
						  const int *num_int_params, int *int_params,
						  int *num_double_params, double *double_params){
	double *lefts, *rights;
	int num_intervals;
	domain_1d(idx, m, x, num_char_params, char_params,
			  num_int_params, int_params, num_double_params, double_params,
			  &num_intervals, &lefts, &rights, NULL);
	if (num_intervals == 0)
		error("In domain_1d_for_R_test(): No feasible point found using domain_1d().\n");
	for (int i = 0; i < num_intervals; i++)
		Rprintf("Interval %d: [%f, %f].\n", i, lefts[i], rights[i]);
}


void search_fused_test(const double *arr, const int *length, const double *x, int *res) {
	*res = (search_fused(arr, *length, *x));
}

void translate_unfuse_test(const double *x, const  int *num_intervals, const double *fused, const double *disp, double *res) {
	*res = translate_unfuse(*x, *num_intervals, fused, disp);
}

void search_unfused_test(const double *lefts, const double *rights, const int *length, const double *x, int *res) {
	*res = search_unfused(lefts, rights, *length, *x);
}

void translate_fuse_test(const double *x, const int *num_intervals, const double *lefts, const double *rights, const double *disp, double *res) {
	*res = translate_fuse(*x, *num_intervals, lefts, rights, disp);
}


