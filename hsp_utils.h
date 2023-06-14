#ifndef HSP_UTILS_H
#define HSP_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include "cliquer.h"

void hsp_log2(int a, int* a1, int* l2);
int hsp_jacobi(int a, int n);
void hsp_paley(int p, int*** Q);
void hsp_allocmat(int*** M, int m, int n);
void hsp_freemat(int*** M, int m, int n);
void hsp_dispmat(int*** M, int m, int n);
int hsp_inArray(int i, int l, int* v);
void hsp_paley_ix(int p, int h, int* ix, int* rx, int*** Q);
int hsp_areOrthogonal(int i, int j, int len, int*** Q);
int hsp_areOrthogonal_ix(int i, int j, int len, int* ix, int*** Q);
int hsp_areOrthogonal_rc(int i, int j, int h, int* row_indices, int* col_indices, int*** Q);
void hsp_set2array(set_t s, int* v);
void hsp_disparray(int** v, int len);
void hsp_residuemat(int la, int lb, int p, int* a, int* b, int*** M);
void hsp_randarray(int a, int b, int len, int* v);

#endif  // HSP_UTILS_H
