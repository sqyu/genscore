//
//  set_ops.h
//
//
//  Written by Shiqing Yu.
//
//

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#endif


void push_op(char *op_stack, int *op_stack_size, const char new_op, int *errno);
char pop_op(char *op_stack, int *op_stack_size, int *errno);
void push_num(int *num_stack, int *num_stack_size, const int new_num, int *errno);
int pop_num(int *num_stack, int *num_stack_size, int *errno);
void shunting_yard(int *num_eqs, const char **infix_pt, char **postfix_pt, int *errno);
void evaluate_logic(const int *num_eqs, const char *postfix,
					int *num_intervals_list, double **lefts_list, double **rights_list,
					int *res_num_intervals, double **res_lefts, double **res_rights,
					int *errno);
void intersection(const int *A_num_intervals, const double *A_lefts, const double *A_rights, const int *B_num_intervals, const double *B_lefts, const double *B_rights, int *res_num_intervals, double **res_lefts, double **res_rights, int *errno);
void merge_sorted_arrays(const int *A_length, const double *A, const int *B_length, const double *B, double **res, int *errno);
void setunion(const int *A_num_intervals, double *A_lefts, double *A_rights, const int *B_num_intervals, double *B_lefts, double *B_rights, int *res_num_intervals, double **res_lefts, double **res_rights, int *errno);
