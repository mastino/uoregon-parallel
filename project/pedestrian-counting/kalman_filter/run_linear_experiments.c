/*
 * test_linear_algebra.c
 * a small program to test the linear algebra code
 
 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * See LICENSE file for licensing information and boring legal stuff

 * If by some miricale you find this software useful, thanks are accepted in
 * the form of chocolate or introductions to potential employers.

*/

#include "linear_algebra.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

void allocate_mats(double** mat_a, double** mat_b, double** mat_c, double** mat_t, long side_len);
void init_mats(double* mat_a, double* mat_b, double* mat_c, double* mat_t, long side_len);
void free_mats(double* mat_a, double* mat_b, double* mat_c, double* mat_t);

int main(int argc, char **argv) {
  
  double *mat_a, *mat_b, *mat_c, *mat_t;
  int n, rows, cols, rows_a, cols_a, cols_b, side_len=100;
  int iterations = 100;
  int i = 0;

  if(argc >= 2) {
    side_len = atoi(argv[1]);
  }
  if(argc >= 3) {
    iterations = atoi(argv[2]);
  }

  n = rows = cols = rows_a = cols_a = cols_b = side_len;
  allocate_mats(&mat_a, &mat_b, &mat_c, &mat_t, side_len);
  if( (mat_a == NULL) || (mat_b == NULL) || (mat_c == NULL) || (mat_t == NULL) ) {
    printf("ERROR: not enough memory\n");
    exit(1);
  }
  init_mats(mat_a, mat_b, mat_c, mat_t, side_len);

  clock_t start = clock();
  for (i = 0; i < iterations; i++) {
    invert_matrix(mat_a, n, mat_c);
    determinant_matrix(mat_a, n);
    cofactor_matrix(mat_a, n, mat_c);
    add_matrix(mat_a, rows, cols, mat_b, mat_c);
    multiply_matrix_by_scalar(mat_a, rows, cols, 5, mat_c);
    multiply_matrix(mat_a, rows_a, cols_a, mat_b, cols_b, mat_c);
    transpose_matrix(mat_a, rows_a, cols_a, mat_c);
    compute_LUP(mat_a, mat_b, mat_c, mat_t, side_len);
    set_zero(mat_a, rows_a, cols_a);
    set_identity(mat_a, rows_a, cols_a);
  }
  clock_t end = clock();
  
  printf("time = %.4f \nside = %d \niter = %d\n", (float)(end - start) / CLOCKS_PER_SEC, side_len, iterations);
  
  free_mats(mat_a, mat_b, mat_c, mat_t);
  return 0;
  
}

void allocate_mats(double** mat_a, double** mat_b, double** mat_c, double** mat_t, long side_len) {

  long tot_mem = side_len * side_len * sizeof(double);

  *mat_a = (double*)malloc(tot_mem);
  *mat_b = (double*)malloc(tot_mem);
  *mat_c = (double*)malloc(tot_mem);
  *mat_t = (double*)malloc(tot_mem);

}

void init_mats(double* mat_a, double* mat_b, double* mat_c, double* mat_t, long side_len) {

  int tot_elm = side_len*side_len;
  long i = 0;

  for (i = 0; i < tot_elm; i++) {
    mat_a[i] = mat_b[i] = mat_c[i] = mat_t[i] = 123.456; // TODO make this random?
  }

}

void free_mats(double* mat_a, double* mat_b, double* mat_c, double* mat_t){
  free(mat_a);
  free(mat_b);
  free(mat_c);
  free(mat_t);
}