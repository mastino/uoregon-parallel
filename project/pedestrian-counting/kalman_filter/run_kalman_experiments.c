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
#include "kalman_filter.h"
#include "linear_algebra.h"

#define FILE_NAME "projectile_motion.csv"

typedef struct Points {
  double* x;
  double* y;
  int size;
} Points;

void test_original();
void test_projectile();
Points get_projectile_measurements(FILE *file);

int main(int argc, char* argv[]) {
  // test_original();
  test_projectile();
  return 0;
}


void test_projectile() {

  FILE* file = fopen(FILE_NAME, "r");

  if( file == NULL ) {
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

  // system dynamics matrix A (nxn)
  // 1  dt 0  0  0  0
  // 0  1  dt 0  0  0
  // 0  0  1  0  0  0
  // 0  0  0  1  dt 0
  // 0  0  0  0  1  dt
  // 0  0  0  0  0  1
  TYPE A_init[] = {1, dt, 0, 0, 0, 0, 0, 1, dt, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, dt, 0, 0, 0, 0, 0, 1, dt, 0, 0, 0, 0, 0, 1};

  // measurement matrix H
  // 1  0  0  0  0  0
  // 0  0  0  1  0  0
  TYPE C_init[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};

  // Reasonable covariance matrices
  // process noise covariance Q
  // 1e-2  0     0     0     0      0
  // 0     5.0   0     0     0      0
  // 0     0     1e-2  0     0      0
  // 0     0     0     1e-2  0      0
  // 0     0     0     0     5.0    0
  // 0     0     0     0     0      1e-2
  TYPE Q_init[] = {1e-2, 0, 0, 0, 0, 0, 0, 5.0, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 5.0, 0, 0, 0, 0, 0, 0, 1e-2};

  // measurement noise covariance R
  // 5 0
  // 0 5
  TYPE R_init[] = {5.0, 0, 0, 5.0};

  // error covariance P
  // 1     0     0     0     0      0
  // 0     1     0     0     0      0
  // 0     0     1     0     0      0
  // 0     0     0     1     0      0
  // 0     0     0     0     1      0
  // 0     0     0     0     0      1   
  TYPE P_init[] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};


  TYPE x_hat_init[] = {0, 0, 0, 0, 0, -9.81};


  TYPE *A, *C, *Q, *R, *P, *K, *x, *y, *x_hat,
       *x_hat_new, *A_T, *C_T, *id,
       *temp_1, *temp_2, *temp_3, *temp_4;

  success = allocate_matrices(&A, &C, &Q, &R, &P, &K, n, m);
  success = success && allocate_vectors(&x, &y, &x_hat, n, m);
  success = success && allocate_temp_matrices(&x_hat_new, &A_T, &C_T, &id,
                                              &temp_1, &temp_2, &temp_3, &temp_4, n, m);

  if( !success ) {
    printf("ERROR allocating matrices\n");
    exit(1);
  }

  // get data
  measurements = get_projectile_measurements(file);

  printf("n: %d\n", measurements.size);
  x_hat_init[0] = measurements.x[0];
  x_hat_init[3] = measurements.y[0];

  copy_mat(A_init, A, n * n);
  copy_mat(C_init, C, n * m);
  copy_mat(Q_init, Q, n * n);
  copy_mat(R_init, R, m * m);
  copy_mat(P_init, P, n * n);
  copy_mat(x_hat_init, x_hat, n);

  printf("\nA:\n");
  print_matrix(A, n, n);
  printf("\nC:\n");
  print_matrix(C, m, n);
  printf("\nQ:\n");
  print_matrix(Q, n, n);
  printf("\nR:\n");
  print_matrix(R, m, m);
  printf("\nP:\n");
  print_matrix(P, n, n);

  printf("t     = %f\n", t);
  printf("x_hat = ");
  print_matrix(x_hat, 1, n);
  printf("\n");

  for(i = 0; i < measurements.size; i++) {
    y[0] = measurements.x[i];
    y[1] = measurements.y[i];

    update(y, x_hat, &t, dt, n, m, A,  C,  Q,  R,  P,  K,
           x_hat_new, A_T, C_T, id, temp_1, temp_2, temp_3, temp_4);
    t += dt;

    printf("t     = %f\n", t);
    printf("y = ");
    print_matrix(y, 1, m);
    printf("x_hat = ");
    print_matrix(x_hat, 1, n);
    printf("\n");
  }

  destroy_matrices(A, C, Q, R, P, K);
  destroy_vectors(x, y, x_hat);
  destroy_temp_matrices(x_hat_new, A_T, C_T, id,
                        temp_1, temp_2, temp_3, temp_4);

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





void test_original() {

  int n = 3; // Number of states
  int m = 1; // Number of measurements
  int i;     // iterator

  double dt = 1.0/30; // Time step
  double t  = 0.0; 

  // Discrete LTI projectile motion, measuring position only
  TYPE A_init[] = {1, dt, 0, 0, 1, dt, 0, 0, 1};
  TYPE C_init[] = {1, 0, 0};

  // Reasonable covariance matrices
  TYPE Q_init[] = {.05, .05, .0, .05, .05, .0, .0, .0, .0};
  TYPE R_init[] = {5};
  TYPE P_init[] = {.1, .1, .1, .1, 10000, 10, .1, 10, 100};

  TYPE *A, *C, *Q, *R, *P, *K, *x, *y, *x_hat,
       *x_hat_new, *A_T, *C_T, *id,
       *temp_1, *temp_2, *temp_3, *temp_4;

  // List of noisy position measurements (y)
  TYPE measurements[] = {
      1.04202710058, 1.10726790452, 1.2913511148, 1.48485250951, 1.72825901034,
      1.74216489744, 2.11672039768, 2.14529225112, 2.16029641405, 2.21269371128,
      2.57709350237, 2.6682215744, 2.51641839428, 2.76034056782, 2.88131780617,
      2.88373786518, 2.9448468727, 2.82866600131, 3.0006601946, 3.12920591669,
      2.858361783, 2.83808170354, 2.68975330958, 2.66533185589, 2.81613499531,
      2.81003612051, 2.88321849354, 2.69789264832, 2.4342229249, 2.23464791825,
      2.30278776224, 2.02069770395, 1.94393985809, 1.82498398739, 1.52526230354,
      1.86967808173, 1.18073207847, 1.10729605087, 0.916168349913, 0.678547664519,
      0.562381751596, 0.355468474885, -0.155607486619, -0.287198661013, -0.602973173813
  };

  int num_measurements = 45;

  TYPE x_hat_init[] = {measurements[0], 0, -9.81};
  char success = 0;

  success = allocate_matrices(&A, &C, &Q, &R, &P, &K, n, m);
  success = success && allocate_vectors(&x, &y, &x_hat, n, m);
  success = success && allocate_temp_matrices(&x_hat_new, &A_T, &C_T, &id,
                                              &temp_1, &temp_2, &temp_3, &temp_4, n, m);

  if( !success ) {
    printf("ERROR allocating matrices\n");
    exit(1);
  }

  copy_mat(A_init, A, n * n);
  copy_mat(C_init, C, n * m);
  copy_mat(Q_init, Q, n * n);
  copy_mat(R_init, R, m * m);
  copy_mat(P_init, P, n * n);
  copy_mat(x_hat_init, x_hat, n);


  printf("\nA:\n");
  print_matrix(A, n, n);
  printf("\nC:\n");
  print_matrix(C, m, n);
  printf("\nQ:\n");
  print_matrix(Q, n, n);
  printf("\nR:\n");
  print_matrix(R, m, m);
  printf("\nP:\n");
  print_matrix(P, n, n);


  printf("t     = %f\n", t);
  printf("x_hat = ");
  print_matrix(x_hat, 1, n);
  printf("\n");

  for(i = 0; i < num_measurements; i++) {
    y[0] = measurements[i];

    update(y, x_hat, &t, dt, n, m, A,  C,  Q,  R,  P,  K,
           x_hat_new, A_T, C_T, id, temp_1, temp_2, temp_3, temp_4);
    t += dt;

    printf("t     = %f\n", t);
    printf("x_hat = ");
    print_matrix(x_hat, 1, n);
    printf("\n");
  }

  destroy_matrices(A, C, Q, R, P, K);
  destroy_vectors(x, y, x_hat);
  destroy_temp_matrices(x_hat_new, A_T, C_T, id,
                        temp_1, temp_2, temp_3, temp_4);

}
