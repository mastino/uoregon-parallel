/*
 * kalman_filter.c
 * super basic Kalman filter implementation
 
 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * See LICENSE file for licensing information and boring legal stuff
 * this code is based heavily on a version by Hayk Martirosyan
 *    https://github.com/hmartiro/kalman-cpp

 * If by some miricale you find this software useful, thanks are accepted in
 * the form of chocolate, coffee, or introductions to potential employers.

 * Things that need to be established
    n - num states
    m - num measurements

    A - system dynamics nxn
    C - H matrix - the measurement one, also output? mxn
    Q - process noise covariance nxn
    R - measurement noise covariance mxm
    P - error covariance nxn
    K - kalman gain nxm

    x     - estimated state n x m
    x_hat - the next prediction n x m
    y     - measurements m

    t  - time
    dt - time step
*/

#include "kalman_filter.h"
#include <omp.h>

char allocate_matrices(TYPE** A, TYPE** C, TYPE** Q, TYPE** R, TYPE** P, TYPE** K, int n, int m) {

  *A = (TYPE*) malloc(n * n * sizeof(TYPE)); 
  *C = (TYPE*) malloc(m * n * sizeof(TYPE));
  *Q = (TYPE*) malloc(n * n * sizeof(TYPE));
  *R = (TYPE*) malloc(m * m * sizeof(TYPE));
  *P = (TYPE*) malloc(n * n * sizeof(TYPE));
  *K = (TYPE*) malloc(n * m * sizeof(TYPE));

  return !( (*A == 0) || (*C == 0) || (*Q == 0) || (*R == 0) || (*P == 0) || (*K == 0) );

}

char allocate_vectors(TYPE** x, TYPE** y, TYPE** x_hat, int n, int m) {
  *x     = (TYPE*) malloc(n * sizeof(TYPE));
  *y     = (TYPE*) malloc(m * sizeof(TYPE));
  *x_hat = (TYPE*) malloc(n * sizeof(TYPE));

  set_zero(*x_hat, n, 1);

  return !( (*x == 0) || (*y == 0) || (*x_hat == 0) );
}

char allocate_temp_matrices(TYPE** x_hat_new, TYPE** A_T, TYPE** C_T, TYPE** id,
                            TYPE** temp_1, TYPE** temp_2, TYPE** temp_3, TYPE** temp_4, TYPE** temp_5,
                            int n, int m) {
  char fail = 0;
  int  size = n > m ? n*n : m*m;

  *x_hat_new = (TYPE*) malloc(n * sizeof(TYPE));     // n x 1
  *A_T       = (TYPE*) malloc(n * n * sizeof(TYPE)); // n x n
  *C_T       = (TYPE*) malloc(n * m * sizeof(TYPE)); // m x n
  *id        = (TYPE*) malloc(n * n * sizeof(TYPE)); // n x n identity  
  set_identity(*id, n, n);
  fail = fail || (x_hat_new == 0) || (A_T == 0) || (C_T == 0) || (id == 0);

  size = size * sizeof(TYPE);
  *temp_1     = (TYPE*) malloc(size); // n x n or m x m if bigger
  *temp_2     = (TYPE*) malloc(size); // n x n or m x m if bigger
  *temp_3     = (TYPE*) malloc(size); // n x n or m x m if bigger
  *temp_4     = (TYPE*) malloc(size); // n x n or m x m if bigger
  *temp_5     = (TYPE*) malloc(size); // n x n or m x m if bigger
  fail = fail || (temp_1 == 0) || (temp_2 == 0) || (temp_3 == 0) || (temp_4 == 0) || (temp_5 == 0);

  return !fail;
}

void destroy_matrices(TYPE* A, TYPE* C, TYPE* Q, TYPE* R, TYPE* P, TYPE* K) {
  free(A);
  free(C);
  free(Q);
  free(R);
  free(P);
  free(K);
}

void destroy_vectors(TYPE* x, TYPE* y, TYPE* x_hat) {
  free(x);
  free(y);
  free(x_hat);
}

void destroy_temp_matrices(TYPE* x_hat_new, TYPE* A_T, TYPE* C_T, TYPE* id,
                           TYPE* temp_1, TYPE* temp_2, TYPE* temp_3, TYPE* temp_4, TYPE* temp_5) {
  free(x_hat_new);
  free(A_T);
  free(C_T);
  free(id);
  free(temp_1);
  free(temp_2);
  free(temp_3);
  free(temp_4);
  free(temp_5);
}

//update the filter
//param y is a vector same size as x and x_hat
//post
void update(TYPE* y, TYPE* x_hat, 
            double* t, double dt, int n, int m,
            TYPE* A, TYPE* C, TYPE* Q, TYPE* R, TYPE* P, TYPE* K,
            TYPE* x_hat_new, TYPE* A_T, TYPE* C_T, TYPE* id,
            TYPE* temp_1, TYPE* temp_2, TYPE* temp_3, TYPE* temp_4, TYPE* temp_5) {

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      // 2
      multiply_matrix(A, n, n, x_hat, 1, x_hat_new);
    }
  
    #pragma omp section
    {
      // 3
      multiply_matrix(A, n, n, P, n, temp_1);
      multiply_matrix(temp_1, n, n, A_T, n, temp_2);
      add_matrix(temp_2, n, n, Q, P);
    }

  
    #pragma omp section
    {
      multiply_matrix(P, n, n, C_T, m, temp_3);
    }
        
    #pragma omp section
    {
      // 6
      multiply_matrix(C, m, n, P, n, temp_1);
      multiply_matrix(temp_1, m, n, C_T, m, temp_2);
      add_matrix(temp_2, m, m, R, temp_1);  
      invert_matrix(temp_1, m, temp_2); // (C*P*C_T+R)^-1
    }
    
  } // first section set

  // 7
  multiply_matrix(temp_3, n, m, temp_2, m, K);
  
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      // 8
      multiply_matrix(C, m, n, x_hat_new, 1, temp_5);
      multiply_matrix_by_scalar(temp_5, m, 1, -1, temp_4);
      add_matrix(y, m, 1, temp_4, temp_5);
      multiply_matrix(K, n, m, temp_5, 1, temp_4);
      add_matrix(x_hat_new, n, 1, temp_4, x_hat);
    }

    #pragma omp section
    {
      // 9
      multiply_matrix(K, n, m, C, n, temp_1);
      multiply_matrix_by_scalar(temp_1, n, n, -1, temp_2);
      add_matrix(id, n, n, temp_2, temp_1);
      multiply_matrix(temp_1, n, n, P, n, temp_2);
      copy_mat(temp_2, P, n * n);
    }
  } // second section set

}
