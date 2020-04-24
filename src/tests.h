//
//  tests.h
//
//
//  Written by Shiqing Yu.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif

void shunting_yard_test(int *num_eqs, const char **infix);
void logic_domain_test(int *num_eqs, const char **infix, int *num_intervals_list,
					   double *lefts, double *rights);
void frac_pow_test(double *num, int *power_numer, int *power_denom, int *abs);
void intersection_test(const int *A_num_intervals, const double *A_lefts, const double *A_rights, const int *B_num_intervals, const double *B_lefts, const double *B_rights);
void merge_sorted_test(const int *A_length, double *A, const int *B_length, double *B);
void setunion_test(const int *A_num_intervals, double *A_lefts, double *A_rights, const int *B_num_intervals, double *B_lefts, double *B_rights);
void domain_1d_for_R_test(int *idx, int *m, double *x,
						  const int *num_char_params, const char **char_params,
						  const int *num_int_params, int *int_params,
						  int *num_double_params, double *double_params);
void search_fused_test(const double *arr, const int *length, const double *x, int *res);
void translate_unfuse_test(const double *x, const  int *num_intervals, const double *fused, const double *disp, double *res);
void search_unfused_test(const double *lefts, const double *rights, const int *length, const double *x, int *res);
void translate_fuse_test(const double *x, const int *num_intervals, const double *lefts, const double *rights, const double *disp, double *res);
