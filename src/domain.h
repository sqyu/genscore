//
//  domain.h
//
//
//  Written by Shiqing Yu.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif

#define TOL 1e-10

void poly_domain_1d(int a, int b, double c, int larger, int abs, int nonnegative,
					  int *num_intervals, double **lefts_pt, double **rights_pt, int *errno);
void domain_1d(const int *idx, const int *p, const double *x,
			   const int *num_char_params, const char **char_params,
			   const int *num_int_params, int *int_params,
			   const int *num_double_params, double *double_params,
			   int *num_intervals, double **lefts_pt, double **rights_pt,
			   double **cache, int *errno);
/*void h(const int *p, const double *x, double *hx, double *hpx, const double *h_pow, const double *h_trunc,
	   const int *num_char_params, const char **char_params,
	   const int *num_int_params, int *int_params,
	   int *num_double_params, double *double_params, int *errno);*/

void dist(const int *n, const int *p, const double *x, double *dists, int *dist_ps,
		  const int *num_char_params, const char **char_params,
		  const int *num_int_params, int *int_params,
		  int *num_double_params, double *double_params, int *errno);

void intersection(const int *A_num_intervals, const double *A_lefts, const double *A_rights, const int *B_num_intervals, const double *B_lefts, const double *B_rights, int *res_num_intervals, double **res_lefts, double **res_rights, int *errno);

void merge_sorted_arrays(const int *A_length, const double *A, const int *B_length, const double *B, double **res, int *errno);

void setunion(const int *A_num_intervals, double *A_lefts, double *A_rights, const int *B_num_intervals, double *B_lefts, double *B_rights, int *res_num_intervals, double **res_lefts, double **res_rights, int *errno);


void poly_domain_1d_for_R(int *a, int *b, double *c, int *larger, int *abs, int *nonnegative, int *num_intervals, double *lefts, double *rights, int *print, int *errno);

void h(const int *p, const double *x, double *hx, double *hpx, const double *h_pow, const double *h_trunc,
	   const int *num_char_params, const char **char_params,
	   const int *num_int_params, int *int_params,
	   int *num_double_params, double *double_params, int *errno);

void fuse_endpoints(const int *num_intervals, const double *lefts, const double *rights, double *fused, double *disp, int *errno);

int binarySearch_fused(const double *arr, int l, int r, double x, int *errno);
int naiveSearch_fused(const double *arr, int length, double x, int *errno);
int search_fused(const double *arr, int length, double x, int *errno);
double translate_unfuse(double x, int num_intervals, const double *fused,
						const double *disp, int *errno);
int search_unfused(const double *lefts, const double *rights, int length, double x, int *errno);
double translate_fuse(double x, int num_intervals, const double *lefts,
					  const double *rights, const double *disp, int *errno);
