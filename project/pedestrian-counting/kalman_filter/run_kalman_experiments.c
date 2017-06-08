/*
 * test-kalman.c
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

    x     - estimated state n x 1
    x_hat - the next prediction n x 1
    y     - measurements m x 1

    t  - time
    dt - time step
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include "kalman_filter.h"
#include "linear_algebra.h"

#define FILE_NAME "../projectile_motion.csv"
#define FILE_NAME_2 "projectile_motion.csv"

typedef struct Points {
  double* x;
  double* y;
  int size;
} Points;

void test_projectile(int iterations);
Points get_projectile_measurements(FILE *file);

int main(int argc, char* argv[]) {
  int iterations = 100;
  if(argc >= 2) {
   iterations = atoi(argv[1]);
  }

  test_projectile(iterations);
  return 0;
}


void test_projectile(int iterations) {

  FILE* file = fopen(FILE_NAME, "r");

  if( file == NULL ) {
    file = fopen(FILE_NAME_2, "r");
  } 
  if(file == NULL ) {
    printf("ERROR open file\n");
    exit(1);
  }

  Points measurements;
  char success = 0;
  int i = 0;

  int n = 6; // Number of states
  int m = 2; // Number of measurements

  double dt = 0.01; // Time step TODO this should probably come from the file
  double t  = 0;
  double time_tot = 0.0;

  TYPE A_init[] = {1, dt, 0, 0, 0, 0, 0, 1, dt, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, dt, 0, 0, 0, 0, 0, 1, dt, 0, 0, 0, 0, 0, 1};
  TYPE C_init[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
  TYPE Q_init[] = {1e-2, 0, 0, 0, 0, 0, 0, 5.0, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 5.0, 0, 0, 0, 0, 0, 0, 1e-2};
  TYPE R_init[] = {5.0, 0, 0, 5.0}; 
  TYPE P_init[] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  TYPE x_hat_init[] = {0, 0, 0, 0, 0, -9.81};

  TYPE *A, *C, *Q, *R, *P, *K, *x, *y, *x_hat,
       *x_hat_new, *A_T, *C_T, *id,
       *temp_1, *temp_2, *temp_3, *temp_4, *temp_5;

  success = allocate_matrices(&A, &C, &Q, &R, &P, &K, n, m);
  success = success && allocate_vectors(&x, &y, &x_hat, n, m);
  success = success && allocate_temp_matrices(&x_hat_new, &A_T, &C_T, &id,
                                              &temp_1, &temp_2, &temp_3, &temp_4, &temp_5, n, m);

  if( !success ) {
    printf("ERROR allocating matrices\n");
    exit(1);
  }

  measurements = get_projectile_measurements(file);

  for (int j = 0; j < iterations; j++) {
    
    x_hat_init[0] = measurements.x[0];
    x_hat_init[3] = measurements.y[0];

    copy_mat(A_init, A, n * n);
    copy_mat(C_init, C, n * m);
    copy_mat(Q_init, Q, n * n);
    copy_mat(R_init, R, m * m);
    copy_mat(P_init, P, n * n);
    copy_mat(x_hat_init, x_hat, n);

    transpose_matrix(A, n, n, A_T);
    transpose_matrix(C, m, n, C_T);

    clock_t start = clock();
    for(i = 0; i < measurements.size; i++) {
      y[0] = measurements.x[i];
      y[1] = measurements.y[i];

      update(y, x_hat, &t, dt, n, m, A,  C,  Q,  R,  P,  K,
             x_hat_new, A_T, C_T, id, temp_1, temp_2, temp_3, temp_4, temp_5);

      t += dt;
    }
    clock_t end = clock();
    time_tot += (float)(end - start) / CLOCKS_PER_SEC;
  }
  
  printf("time = %.4f \niter = %d\n", time_tot, iterations);

  destroy_matrices(A, C, Q, R, P, K);
  destroy_vectors(x, y, x_hat);
  destroy_temp_matrices(x_hat_new, A_T, C_T, id,
                        temp_1, temp_2, temp_3, temp_4, temp_5);

  free(measurements.x);
  free(measurements.y);
}

Points get_projectile_measurements(FILE *file) {

  int n, i, j;
  char* tok;

  char line[1024];

  fgets(line, 1024, file); // header
  fgets(line, 1024, file); 
  n = atoi(line);

  Points data_in;
  data_in.size = n;
  data_in.x = malloc(n * sizeof(TYPE));
  data_in.y = malloc(n * sizeof(TYPE));

  i = 0;
  while (fgets(line, 1024, file)) { //t

    tok = strtok(line, ",");
    for (j = 0; j < 13; j++) {
      if (tok == NULL)
        break;

      if(j == 2) {
        data_in.x[i] = atof(tok);

      } else if(j == 8) {
        data_in.y[i] = atoi(tok);
        break;
      }
      tok = strtok(NULL, ",");
    }

    i++;
  }

  return data_in;
}



