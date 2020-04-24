/* 
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
#include "set_ops.h"

#define TOL 1e-6
#define MAXOPSTACK 2999
#define MAXNUMSTACK 999


void push_op(char *op_stack, int *op_stack_size, const char new_op){
	if (*op_stack_size >= MAXOPSTACK)
		error("In parsing notation: operator stack size (%d) over max limit (%d).\n", *op_stack_size, MAXOPSTACK);
	op_stack[(*op_stack_size)++] = new_op;
}

char pop_op(char *op_stack, int *op_stack_size){
	if (*op_stack_size <= 0)
		error("In parsing notation: operator stack size 0, nothing to pop.\n");
	return op_stack[--(*op_stack_size)];
}

void push_num(int *num_stack, int *num_stack_size, const int new_num){
	if (*num_stack_size >= MAXNUMSTACK)
		error("In parsing notation: number stack size (%d) over max limit (%d).\n", *num_stack_size, MAXNUMSTACK);
	num_stack[(*num_stack_size)++] = new_num;
}

int pop_num(int *num_stack, int *num_stack_size){
	if (*num_stack_size <= 0)
		error("In parsing notation: number stack size 0, nothing to pop.\n");
	return num_stack[--(*num_stack_size)];
}

int eq_nums_len(int *num_eqs) {
	// Length of the string with 1 through *num_eqs concatenated
	// Ex: 12 -> 123456789101112 -> returns 15
	if (*num_eqs > MAXNUMSTACK)
		error("Number of equations (%d) exceeded max number (%d) allowed.\n", *num_eqs, MAXNUMSTACK);
	if (*num_eqs < 10) return *num_eqs;
	if (*num_eqs < 100) return (*num_eqs - 9) * 2 + 9;
	return (*num_eqs - 99) * 3 + 189;
}

void shunting_yard(int *num_eqs, const char **infix_pt, char **postfix_pt, int *postfix_length) {
	/*
	 Translates an infix logic notation into a post-fix notation, where the infix string contains equation numbers from 1 through *num_eqs, + for union, * for intersection, spaces and parentheses only.
	 For example, ((1 & 2) | 3) & 4 becomes 1 2 & 3 | 4 &.
	 Note that chained operators with different types are not allowed, e.g. 1 & 2 | 3. Parentheses must be added to avoid ambiguity, i.e. (1 & 2) | 3, 1 & (2 | 3) are allowed and have different meanings.
	 Note that parentheses must match.
	 
	 num_eqs: number of equations
	 infix: infix logic notation.
	 postfix_pt: pointer to the pointer for the post-fix notation output.
	 */
	if (*num_eqs > MAXNUMSTACK)
		error("Number of equations (%d) exceeded max number (%d) allowed.\n", *num_eqs, MAXNUMSTACK);
	const char *infix = *infix_pt;
	int infix_len = strlen(infix), postfix_len = 0;
	char *postfix = *postfix_pt; //(char*)malloc((eq_nums_len(num_eqs) + 2*(*num_eqs)) * sizeof(char)); // Size of postfix cannot exceed infix
	//*postfix_pt = postfix;
	char *op_stack = (char*)malloc(MAXOPSTACK * sizeof(char));
	int pos = 0, op_stack_size = 0, prev_is = 0; // prev_is: 0: start, 1: num, 2: | or &, 3: ( or )
	while (pos < infix_len) {
		if (isspace(infix[pos])) {
			pos++;
			continue;
		} else if (isdigit(infix[pos])) {
			if (prev_is == 1) {
				Rprintf("In parsing notation: Scanned '");
				for (int i = 0; i <= pos; i++)
					Rprintf("%c", infix[i]);
				Rprintf("'.\n");
				error("Numbers cannot be directly followed by another number (e.g. '12 34').\n");
			}
			int this_num = 0; // For sanity check if number larger than num_eqs
			while (pos < infix_len && isdigit(infix[pos])) {
				this_num = this_num * 10 + (int)(infix[pos] - 48); // Finish reading the number
				postfix[postfix_len++] = infix[pos++];
			}
			if (this_num <= 0 || this_num > *num_eqs)
				error("In parsing notation: Equation %d out of range. Equation number must be in [1, %d] since you specified %d equations.\n", this_num, *num_eqs, *num_eqs);
			postfix[postfix_len++] = ' ';
			prev_is = 1;
		} else if (infix[pos] == '|' || infix[pos] == '&') {
			if (prev_is == 2) {
				Rprintf("In parsing notation: Scanned '");
				for (int i = 0; i <= pos; i++)
					Rprintf("%c", infix[i]);
				Rprintf("'.\n");
				error("Operations cannot be directly followed by another operation (e.g. '1 & | 2').\n");
			}
			if (postfix_len < 1) {
				error("In parsing notation: The string cannot start with an operation, and must start with a number instead.\n");
			}
			while (op_stack_size && op_stack[op_stack_size - 1] != '(') { // If last op | or &
				if (infix[pos] != op_stack[op_stack_size - 1])  // If the previous operator is also | or & but different from the current one -- ambiguous notation and parenthese are required
					error("In parsing notation: Ambiguous notation; for chained operations of &/| parenthese required, unless they are of the same time. E.g. '1 & 2 & 3' is okay but '1 & 2 | 3' is not allowed; '(1 & 2) | 3' OR '1 & (2 | 3)' must be used.\n");
				postfix[postfix_len++] = pop_op(op_stack, &op_stack_size);
			}
			push_op(op_stack, &op_stack_size, infix[pos]);
			pos++;
			prev_is = 2;
		} else if (infix[pos] == '(') {
			push_op(op_stack, &op_stack_size, '(');
			pos++;
			prev_is = 3;
		} else if (infix[pos] == ')') {
			while (op_stack_size && op_stack[op_stack_size - 1] != '(')
				postfix[postfix_len++] = pop_op(op_stack, &op_stack_size);
			if (op_stack_size == 0) { // No ( found
				Rprintf("In parsing notation: Mismatched parentheses: extra right parenthesis after '");
				for (int i = 0; i < pos; i++)
					Rprintf("%c", infix[i]);
				Rprintf("'.\n");
				error("Please check your input.\n", pos + 1);
			}
			pos++;
			prev_is = 3;
			op_stack_size--; // Discard the (
		} else
			error("In parsing notation: Invalid character: %c.\n", infix[pos]);
	}
	while (op_stack_size > 0) {
		if (op_stack[op_stack_size - 1] == '(')
			error("In parsing notation: Mismatched parentheses (extra left parenthesis unmatched). Please check your input.\n");
		postfix[postfix_len++] = pop_op(op_stack, &op_stack_size);
	}
	postfix[postfix_len] = '\0';
	free(op_stack);
	*postfix_length = postfix_len;
	return;
}

void evaluate_logic(const int *num_eqs, const char *postfix,
					int *num_intervals_list, double **lefts_list, double **rights_list,
					int *res_num_intervals, double **res_lefts, double **res_rights) {
	/* Assumes that the max number in postfix appeared is not negative and does not exceed num_eqs.
	 This is checked in shunting_yard(). Note that non-negativity cannot be checked here since new
	 intersected/merged domains will be given negative numbers as their index.
	 */
	
	*res_num_intervals = 0;
	*res_lefts = NULL;
	*res_rights = NULL;
	int postfix_len = strlen(postfix), no_operator = 1;
	int pos = 0, num_stack_size = 0, *num_stack = (int*)malloc(MAXOPSTACK * sizeof(int));
	double **intermediate_lefts_list = (double**)malloc(*num_eqs * sizeof(double *)); // Stores pointers to the intermediate domains; number of intermediate results cannot exceed string length
	double **intermediate_rights_list = (double**)malloc(*num_eqs * sizeof(double *));
	int *intermediate_num_intervals = (int*)malloc(*num_eqs * sizeof(int));
	int num_intermediate_results = 0; // Counts number of intermediate domains stored
	while (pos < postfix_len) {
		/*////Rprintf("Now reading %c. Num stack: ", postfix[pos]); //////
		 for (int j = 0; j < num_stack_size; j++)
		 Rprintf("%d ", num_stack[j]);
		 Rprintf("\n");*/
		if (isspace(postfix[pos])) {
			pos++;
			continue;
		} else if (isdigit(postfix[pos])) {
			int this_num = 0;
			while (pos < postfix_len && isdigit(postfix[pos]))
				this_num = this_num * 10 + (int)(postfix[pos++] - 48); // Finish reading the number
			push_num(num_stack, &num_stack_size, this_num);
		} else if (postfix[pos] == '|' || postfix[pos] == '&') {
			no_operator = 0;
			if (num_stack_size < 2) {
				Rprintf("In evaluating notation: There should be at least two numbers before an operator. Got '");
				for (int i = 0; i <= pos; i++)
					Rprintf("%c", postfix[i]);
				Rprintf("'.\n");
				error("Please check your original input.\n");
			}
			int eqs[2] = {pop_num(num_stack, &num_stack_size),
				pop_num(num_stack, &num_stack_size)};
			int num_intervals[2], new_num_intervals;
			double *lefts[2], *rights[2], *new_lefts_pt, *new_rights_pt;
			for (int i = 0; i < 2; i++) {
				if (eqs[i] == 0 || eqs[i] < -num_intermediate_results || eqs[i] > *num_eqs)
					error("In evaluating notation: Equation %d out of range (must be in [%d, -1] or [1, %d]).\n", eqs[i], -num_intermediate_results, *num_eqs);
				if (eqs[i] > 0) {
					num_intervals[i] = num_intervals_list[eqs[i] - 1];
					lefts[i] = lefts_list[eqs[i] - 1];
					rights[i] = rights_list[eqs[i] - 1];
				} else {
					num_intervals[i] = intermediate_num_intervals[-eqs[i] - 1];
					lefts[i] = intermediate_lefts_list[-eqs[i] - 1];
					rights[i] = intermediate_rights_list[-eqs[i] - 1];
				}
			}
			
			if (postfix[pos] == '|') {
				setunion(num_intervals, lefts[0], rights[0],
						 num_intervals + 1, lefts[1], rights[1],
						 &new_num_intervals, &new_lefts_pt, &new_rights_pt);
			} else {
				intersection(num_intervals, lefts[0], rights[0],
							 num_intervals + 1, lefts[1], rights[1],
							 &new_num_intervals, &new_lefts_pt, &new_rights_pt);
			}
			intermediate_lefts_list[num_intermediate_results] = new_lefts_pt;
			intermediate_rights_list[num_intermediate_results] = new_rights_pt;
			intermediate_num_intervals[num_intermediate_results] = new_num_intervals;
			num_intermediate_results++;
			push_num(num_stack, &num_stack_size, -num_intermediate_results);
			pos++;
		} else
			error("In evaluating notation: Invalid character in postfix: %c.\n", postfix[pos]);
	}
	// After consuming the string, there should be exactly one equation number left in the queue
	if (num_stack_size != 1)
		error("In evaluating notation: stack size should be exactly 1 after the whole string is consumed. Got size %d.\n", num_stack_size);
	int res_eq_number = pop_num(num_stack, &num_stack_size);
	// If the last number is an original domain (which happens ONLY when there is no operator, i.e. the whole string is just a number)
	if (res_eq_number == 0 || res_eq_number < -num_intermediate_results || res_eq_number > *num_eqs) {
		if (no_operator)
			error("In evaluating notation: Remaining equation number of range (must be in [1, %d]).\n", *num_eqs);
		else
			error("In evaluating notation: Remaining equation number of range (must be in [%d, -1] or [1, %d]).\n", -num_intermediate_results, *num_eqs);
	}
	if (res_eq_number >= 0) {
		if (no_operator == 0)
			error("In evaluating notation: There is only one number left in the stack, but the original string does contain an operator.\n");
		*res_num_intervals = num_intervals_list[res_eq_number - 1];
		*res_lefts = (double *)malloc(*res_num_intervals * sizeof(double));
		*res_rights = (double *)malloc(*res_num_intervals * sizeof(double));
		for (int j = 0; j < *res_num_intervals; j++) {
			(*res_lefts)[j] = lefts_list[res_eq_number - 1][j];
			(*res_rights)[j] = rights_list[res_eq_number - 1][j];
		}
	} else {
		// Otherwise it is an intermediate result
		*res_num_intervals = intermediate_num_intervals[-res_eq_number - 1];
		*res_lefts = (double *)malloc(*res_num_intervals * sizeof(double));
		*res_rights = (double *)malloc(*res_num_intervals * sizeof(double));
		for (int j = 0; j < *res_num_intervals; j++) {
			(*res_lefts)[j] = intermediate_lefts_list[-res_eq_number - 1][j];
			(*res_rights)[j] = intermediate_rights_list[-res_eq_number - 1][j];
		}
	}
	for (int j = 0; j < num_intermediate_results; j++) {
		//Rprintf("evl(): inter_left[%d]: %p  inter_right[%d]: %p\n", j, (intermediate_lefts_list[j]), j, (intermediate_rights_list[j]));
		free(intermediate_lefts_list[j]); free(intermediate_rights_list[j]);
	}
	free(num_stack); free(intermediate_num_intervals);
	free(intermediate_lefts_list); free(intermediate_rights_list);
}


void intersection(const int *A_num_intervals, const double *A_lefts, const double *A_rights, const int *B_num_intervals, const double *B_lefts, const double *B_rights, int *res_num_intervals, double **res_lefts, double **res_rights){
	/*
	 Calculates the intersections between a list of intervals specified by A_lefts and A_rights with those specified by B_lefts and B_rights. Assumes the intervals within each group are disjoint and sorted.
	 */
	// Max number of intersected intervals cannot exceed total numbers of intervals
	int ai = 0, bi = 0, max_num_intervals = (*A_num_intervals + *B_num_intervals);
	*res_lefts = (double*)malloc(max_num_intervals * sizeof(double));
	*res_rights = (double*)malloc(max_num_intervals * sizeof(double));
	*res_num_intervals = 0;
	while (ai < *A_num_intervals && bi < *B_num_intervals) {
		if (A_rights[ai] > B_lefts[bi] && B_rights[bi] > A_lefts[ai]) {
			(*res_lefts)[*res_num_intervals] = fmax(A_lefts[ai], B_lefts[bi]);
			(*res_rights)[(*res_num_intervals)++] = fmin(A_rights[ai], B_rights[bi]);
		}
		if (A_rights[ai] > B_rights[bi])
			bi++;
		else if (A_rights[ai] < B_rights[bi])
			ai++;
		else {
			ai++; bi++;
		}
	}
}

void merge_sorted_arrays(const int *A_length, const double *A, const int *B_length, const double *B, double **res){
	// Not using length of result as it's always *A_length + *B_length
	*res = (double*)malloc((*A_length + *B_length) * sizeof(double));
	int ai = 0, bi = 0, res_i = 0;
	while (ai < *A_length && bi < *B_length)
		(*res)[res_i++] = A[ai] <= B[bi] ? A[ai++] : B[bi++];
	while (ai < *A_length)
		(*res)[res_i++] = A[ai++];
	while (bi < *B_length)
		(*res)[res_i++] = B[bi++];
}

void setunion(const int *A_num_intervals, double *A_lefts, double *A_rights, const int *B_num_intervals, double *B_lefts, double *B_rights, int *res_num_intervals, double **res_lefts, double **res_rights){
	/*
	 Calculates the intersections between a list of intervals specified by A_lefts and A_rights with those specified by B_lefts and B_rights. Assumes the intervals within each group are disjoint and sorted.
	 */
	// Max number of intersected intervals cannot exceed total numbers of intervals
	int max_num_intervals = (*A_num_intervals + *B_num_intervals);
	*res_lefts = (double*)malloc(max_num_intervals * sizeof(double));
	*res_rights = (double*)malloc(max_num_intervals * sizeof(double));
	
	double *merged_lefts, *merged_rights;
	if (*A_num_intervals == 0) {
		*res_num_intervals = *B_num_intervals;
		for (int i = 0; i < *res_num_intervals; i++) {
			(*res_lefts)[i] = B_lefts[i]; (*res_rights)[i] = B_rights[i];
		}
		return;
	} else if (*B_num_intervals == 0) {
		*res_num_intervals = *A_num_intervals;
		for (int i = 0; i < *res_num_intervals; i++) {
			(*res_lefts)[i] = A_lefts[i]; (*res_rights)[i] = A_rights[i];
		}
		return;
	}
	merge_sorted_arrays(A_num_intervals, A_lefts, B_num_intervals, B_lefts, &merged_lefts);
	merge_sorted_arrays(A_num_intervals, A_rights, B_num_intervals, B_rights, &merged_rights);
	
	(*res_lefts)[0] = merged_lefts[0];
	(*res_rights)[0] = merged_rights[0];
	*res_num_intervals = 1;
	
	for (int i = 1; i < max_num_intervals; i++) {
		if (merged_lefts[i] <= (*res_rights)[*res_num_intervals - 1])
			(*res_rights)[*res_num_intervals - 1] = merged_rights[i];
		else {
			(*res_lefts)[*res_num_intervals] = merged_lefts[i];
			(*res_rights)[*res_num_intervals] = merged_rights[i];
			(*res_num_intervals)++;
		}
	}
	free(merged_lefts); free(merged_rights);
}
