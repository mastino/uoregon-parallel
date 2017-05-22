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

void test_inverse(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_cofactor(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_multiply(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_add(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_transpose(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_determinant(double* mat_a, long side_len);
void test_zero_and_id(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_compute_LUP(double* mat_a, double* mat_b, double* mat_c, long side_len);
void test_compute_LUP_inline(double* mat_a, double* mat_b, double* mat_c, long side_len);

void allocate_mats(double** mat_a, double** mat_b, double** mat_c, long side_len);
void init_mats(double* mat_a, double* mat_b, double* mat_c, long side_len);
void free_mats(double* mat_a, double* mat_b, double* mat_c);

int main(int argc, char **argv) {
  
  double *mat_a, *mat_b, *mat_c;
  long side_len = 100; //TODO get from command line

  if(argc >= 2) {
    side_len = atoi(argv[1]);
  }

  allocate_mats(&mat_a, &mat_b, &mat_c, side_len);
  if( (mat_a == NULL) || (mat_b == NULL) || (mat_c == NULL) ) {
    printf("ERROR: not enough memory\n");
    exit(1);
  }
  init_mats(mat_a, mat_b, mat_c, side_len);

  // test_zero_and_id();
  // test_inverse();
  // test_cofactor();
  test_determinant(mat_a, side_len);
  // test_transpose();
  // test_add();
  // test_multiply();
  // test_compute_LUP();
  // test_compute_LUP_inline();

  free_mats(mat_a, mat_b, mat_c);
  return 0;
}

void test_add(double* mat_a, double* mat_b, double* mat_c, long side_len) {
  int col_A = 3, row_A = 3;
  double A[] = {2,2,2,
                2,2,2,
                2,2,2};
                
  int col_B = 2, row_B = 3;
  double B[] = {2,2,
                2,2,
                2,2};
                
  int col_C = 2, row_C = 3;
  double C[] = {2,2,
                2,2,
                2,2};


  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nB:\n");
  print_matrix(B, row_B, col_B);
  printf("\nC:\n");
  print_matrix(C, row_C, col_C);
  printf("\n");

  add_matrix(A, row_A, col_A, B, C);
  printf("\nC <- A + B:\n");
  print_matrix(C, row_C, col_C);
  printf("\n");
}

void test_multiply(double* mat_a, double* mat_b, double* mat_c, long side_len) {
  int col_A = 3, row_A = 3;
  double A[] = {1,2,3,
                4,5,6,
                7,8,9};
                
  int col_B = 2, row_B = 3;
  double B[] = {1,2,
                3,4,
                5,6};
                
  int col_C = 2, row_C = 3;
  double C[] = {2,2,
                2,2,
                2,2};


  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nB:\n");
  print_matrix(B, row_B, col_B);
  printf("\nC:\n");
  print_matrix(C, row_C, col_C);
  printf("\n");

  multiply_matrix(A, row_A, col_A, B, col_B, C);
  printf("\nC <- A * B:\n");
  print_matrix(C, row_C, col_C);
  printf("\n");
}

void test_transpose(double* mat_a, double* mat_b, double* mat_c, long side_len) {
  int col_A = 3, row_A = 4;
  double A[] = {1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10,11,12};
                                
  int col_C = 4, row_C = 3;
  double C[] = {2,2,2,2,
                2,2,2,2,
                2,2,2,2};


  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nC:\n");
  print_matrix(C, row_C, col_C);
  printf("\n");

  transpose_matrix(A, row_A, col_A, C);
  printf("\nC <- A^T:\n");
  print_matrix(C, row_C, col_C);
  printf("\n");
}

void test_zero_and_id(double* mat_a, double* mat_b, double* mat_c, long side_len) {
  int col_A = 3, row_A = 3;
  double A[] = {1, 2, 3,
                4, 5, 6,
                7, 8, 9};


  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\n");

  set_zero(A, row_A, col_A);
  printf("\nzero:\n");
  print_matrix(A, row_A, col_A);
  printf("\n");
  set_identity(A, row_A, col_A);
  printf("\nidentity:\n");
  print_matrix(A, row_A, col_A);
  printf("\n");
}

void test_determinant(double* mat_a, long side_len) {

  determinant_matrix(mat_a, side_len);

}


void test_cofactor(double* mat_a, double* mat_b, double* mat_c, long side_len) {

  int col_A = 4, row_A = 4;
  double A[] = {3,0,2,-1,
                1,2,0,-2,
                4,0,6,-3,
                5,0,2, 0};

  double expected[] = {12,-50,-30,-44,
                0,10,0,0,
                -4,10,10,8,
                0,20,10,20};

  double result[] = {3,0,2,-1,
                1,2,0,-2,
                4,0,6,-3,
                5,0,2, 0};


  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nExpected cofactor:\n");
  print_matrix(expected, row_A, col_A);
  cofactor_matrix(A, row_A, result);
  printf("\nCalculated cofactor:\n");
  print_matrix(result, row_A, col_A);
  printf("\n");

}


void test_inverse(double* mat_a, double* mat_b, double* mat_c, long side_len) {

  int col_A = 4, row_A = 4;
  double A[] = {458.1233,0,-1,0,
                0,458.1233,0,0,
                0,0,1.63,0,
                0,0,0,1.63};

  double expected[] = {0.00218,0,0.00134,0,
                      0.0,0.00218,0.0,0,
                      0.0,0,0.613,0.0,
                      0.0,0,0.0,0.613};

  double result[] = {3,0,2,-1,
                1,2,0,-2,
                4,0,6,-3,
                5,0,2, 0};


  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nExpected inverse:\n");
  print_matrix(expected, row_A, col_A);
  invert_matrix(A, row_A, result);
  printf("\nCalculated inverse:\n");
  print_matrix(result, row_A, col_A);
  printf("\n");

}


void test_compute_LUP(double* mat_a, double* mat_b, double* mat_c, long side_len) {

  int col_A = 4, row_A = 4, tot=16;
  double A[] = {11,9,24,2,
                1,5,2,6,
                3,17,18,1,
                2,5,7, 1};

  double exp_L[] = {1,0,0,0,
                    0.27273,1,0,0,
                    0.09091,0.2875,1,0,
                    0.18182,0.23125,0.0036,1};


  double exp_U[] = {11,9,24,2,
                    0,14.54545,11.45455,0.45455,
                    0,0,-3.475,5.6875,
                    0,0,0,0.51079};


  double exp_P[] = {1,0,0,0,
                    0,0,1,0,
                    0,1,0,0,
                    0,0,0,1};

  double L[tot], U[tot], P[tot];
  int num_pivots;


  num_pivots = compute_LUP(A, L, U, P, col_A);

  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nExpected L:\n");
  print_matrix(exp_L, row_A, col_A);
  printf("\nCalculated L:\n");
  print_matrix(L, row_A, col_A);
  printf("\n");
  printf("\nExpected U:\n");
  print_matrix(exp_U, row_A, col_A);
  printf("\nCalculated U:\n");
  print_matrix(U, row_A, col_A);
  printf("\n");
  printf("\nExpected P:\n");
  print_matrix(exp_P, row_A, col_A);
  printf("\nCalculated P:\n");
  print_matrix(P, row_A, col_A);
  printf("\n");
  printf("\nExpected odd number of pivots\n");
  printf("Num pivots used: %d\n", num_pivots);
  printf("\n");

}


void test_compute_LUP_inline(double* mat_a, double* mat_b, double* mat_c, long side_len) {

  int col_A = 4, row_A = 4, tot=16;
  double A[] = {11,9,24,2,
                1,5,2,6,
                3,17,18,1,
                2,5,7, 1};

  double exp_L[] = {1,0,0,0,
                    0.27273,1,0,0,
                    0.09091,0.2875,1,0,
                    0.18182,0.23125,0.0036,1};


  double exp_U[] = {11,9,24,2,
                    0,14.54545,11.45455,0.45455,
                    0,0,-3.475,5.6875,
                    0,0,0,0.51079};


  double exp_P[] = {1,0,0,0,
                    0,0,1,0,
                    0,1,0,0,
                    0,0,0,1};

  double L[tot], U[tot], P[tot];
  int num_pivots;


  num_pivots = compute_LUP_inline(A, L, U, P, col_A);

  printf("\nA:\n");
  print_matrix(A, row_A, col_A);
  printf("\nExpected L:\n");
  print_matrix(exp_L, row_A, col_A);
  printf("\nCalculated L:\n");
  print_matrix(L, row_A, col_A);
  printf("\n");
  printf("\nExpected U:\n");
  print_matrix(exp_U, row_A, col_A);
  printf("\nCalculated U:\n");
  print_matrix(U, row_A, col_A);
  printf("\n");
  printf("\nExpected P:\n");
  print_matrix(exp_P, row_A, col_A);
  printf("\nCalculated P:\n");
  print_matrix(P, row_A, col_A);
  printf("\n");
  printf("\nExpected odd number of pivots\n");
  printf("Num pivots used: %d\n", num_pivots);
  printf("\n");

}

void allocate_mats(double** mat_a, double** mat_b, double** mat_c, long side_len) {

  long tot_mem = side_len * side_len * sizeof(double);

  *mat_a = (double*)malloc(tot_mem);
  *mat_b = (double*)malloc(tot_mem);
  *mat_c = (double*)malloc(tot_mem);

}

void init_mats(double* mat_a, double* mat_b, double* mat_c, long side_len) {

  int tot_elm = side_len*side_len;
  long i = 0;

  for (i = 0; i < tot_elm; i++) {
    mat_a[i] = mat_b[i] = mat_c[i] = 123.456; // TODO make this random?
  }

}

void free_mats(double* mat_a, double* mat_b, double* mat_c){
  free(mat_a);
  free(mat_b);
  free(mat_c);
}