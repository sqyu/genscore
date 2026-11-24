#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void dist(const int *, const int *, const double *, double *, int *, const int *, const char **, const int *, int *, int *, double *);
extern void elts_ab_c(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_ab_np(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_ab_p(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_exp_c(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_exp_np(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_exp_p(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_gamma_np(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_gamma_p(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_gauss_c(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_gauss_np(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_gauss_p(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_loglog_c(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_loglog_np(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_loglog_p(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_loglog_simplex_c(int *, int *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void elts_loglog_simplex_np(int *, int *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void full(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, double *, int *, int *, double *, int *, double *, int *);
extern void full_asymm(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, double *, int *, int *, double *, int *, double *, int *);
extern void poly_domain_1d_for_R(int *, int *, double *, int *, int *, int *, int *, double *, double *, int *);
extern void profiled(int *, double *, double *, double *, double *, double *, int *, int *, int *, double *, int *, double *, int *, double *, int *);
extern void profiled_asymm(int *, double *, double *, double *, double *, double *, int *, int *, int *, double *, int *, double *, int *, double *, int *);
extern void rab_arms(const int *, const int *, const int *, const int *, const int *, const int *, const int *, const int *, const int *, double *, double *, const double *, const double *, double *, const int *, const char **, const int *, int *, int *, double *, int *);
extern void shunting_yard(int *, const char **, char **, int *);
extern void simplex_centered(int *, int *, double *, double *, double *, double *, double *, double *, int *, int *, int *, double *, int *, double *, int *, double *);
extern void simplex_full(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, double *, int *, int *, double *, int *, double *);

static const R_CMethodDef CEntries[] = {
	{"dist",                   (DL_FUNC) &dist,                   11},
	{"elts_ab_c",              (DL_FUNC) &elts_ab_c,              10},
	{"elts_ab_np",             (DL_FUNC) &elts_ab_np,             14},
	{"elts_ab_p",              (DL_FUNC) &elts_ab_p,              14},
	{"elts_exp_c",             (DL_FUNC) &elts_exp_c,             11},
	{"elts_exp_np",            (DL_FUNC) &elts_exp_np,            12},
	{"elts_exp_p",             (DL_FUNC) &elts_exp_p,             12},
	{"elts_gamma_np",          (DL_FUNC) &elts_gamma_np,          12},
	{"elts_gamma_p",           (DL_FUNC) &elts_gamma_p,           12},
	{"elts_gauss_c",           (DL_FUNC) &elts_gauss_c,           10},
	{"elts_gauss_np",          (DL_FUNC) &elts_gauss_np,          12},
	{"elts_gauss_p",           (DL_FUNC) &elts_gauss_p,           12},
	{"elts_loglog_c",          (DL_FUNC) &elts_loglog_c,          13},
	{"elts_loglog_np",         (DL_FUNC) &elts_loglog_np,         12},
	{"elts_loglog_p",          (DL_FUNC) &elts_loglog_p,          12},
	{"elts_loglog_simplex_c",  (DL_FUNC) &elts_loglog_simplex_c,  20},
	{"elts_loglog_simplex_np", (DL_FUNC) &elts_loglog_simplex_np, 17},
	{"full",                   (DL_FUNC) &full,                   21},
	{"full_asymm",             (DL_FUNC) &full_asymm,             21},
	{"poly_domain_1d_for_R",   (DL_FUNC) &poly_domain_1d_for_R,   10},
	{"profiled",               (DL_FUNC) &profiled,               15},
	{"profiled_asymm",         (DL_FUNC) &profiled_asymm,         15},
	{"rab_arms",			   (DL_FUNC) &rab_arms,				  21},
	{"shunting_yard",          (DL_FUNC) &shunting_yard,           4},
	{"simplex_centered",       (DL_FUNC) &simplex_centered,       16},
	{"simplex_full",           (DL_FUNC) &simplex_full,           25},
	{NULL, NULL, 0}
};

void R_init_genscore(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
