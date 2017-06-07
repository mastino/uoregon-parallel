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
                            TYPE** temp_1, TYPE** temp_2, TYPE** temp_3, TYPE** temp_4, int n, int m) {
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
  fail = fail || (temp_1 == 0) || (temp_2 == 0) || (temp_3 == 0) || (temp_4 == 0);

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
                           TYPE* temp_1, TYPE* temp_2, TYPE* temp_3, TYPE* temp_4) {
  free(x_hat_new);
  free(A_T);
  free(C_T);
  free(id);
  free(temp_1);
  free(temp_2);
  free(temp_3);
  free(temp_4);
}

//update the filter
//param y is a vector same size as x and x_hat
//post
void update(TYPE* y, TYPE* x_hat, 
            double* t, double dt, int n, int m,
            TYPE* A, TYPE* C, TYPE* Q, TYPE* R, TYPE* P, TYPE* K,
            TYPE* x_hat_new, TYPE* A_T, TYPE* C_T, TYPE* id,
            TYPE* temp_1, TYPE* temp_2, TYPE* temp_3, TYPE* temp_4) {

  predict(x_hat, n, m, A, Q, P, 
          x_hat_new, A_T, temp_1, temp_2);
  correct(y, x_hat, n, m, C, R, P, K,
          x_hat_new, C_T, id, temp_1, temp_2, temp_3, temp_4);

}

//predict the state and update P
//param
//post
void predict(TYPE* x_hat, 
            int n, int m,
            TYPE* A, TYPE* Q, TYPE* P,
            TYPE* x_hat_new, TYPE* A_T,
            TYPE* temp_1, TYPE* temp_2) {

  transpose_matrix(A, n, n, A_T); // do this separately since they never? change

  //x_hat_new = A * x_hat
  multiply_matrix(A, n, n, x_hat, 1, x_hat_new);

  //P = A*P*A_T + Q;
  multiply_matrix(A, n, n, P, n, temp_1);
  multiply_matrix(temp_1, n, n, A_T, n, temp_2);
  add_matrix(temp_2, n, n, Q, P);

}

//correct the filter based on measurement
//param 
//post
void correct(TYPE* y, TYPE* x_hat, 
            int n, int m,
            TYPE* C, TYPE* R, TYPE* P, TYPE* K,
            TYPE* x_hat_new, TYPE* C_T, TYPE* id,
            TYPE* temp_1, TYPE* temp_2, TYPE* temp_3, TYPE* temp_4) {


  transpose_matrix(C, m, n, C_T); // do this separately since they never? change  

  // K = P*C_T*(C*P*C_T+R)^-1
  multiply_matrix(C, m, n, P, n, temp_1);
  multiply_matrix(temp_1, m, n, C_T, m, temp_2);
  add_matrix(temp_2, m, m, R, temp_1);  
  invert_matrix(temp_1, m, temp_2); // (C*P*C_T+R)^-1
  multiply_matrix(P, n, n, C_T, m, temp_1); // P*C_T
  multiply_matrix(temp_1, n, m, temp_2, m, K);

  // x_hat = x_hat_new + K * (y - C*x_hat_new);
  multiply_matrix(C, m, n, x_hat_new, 1, temp_3);
  multiply_matrix_by_scalar(temp_3, m, 1, -1, temp_4);
  add_matrix(y, m, 1, temp_4, temp_3);
  multiply_matrix(K, n, m, temp_3, 1, temp_4);
  add_matrix(x_hat_new, n, 1, temp_4, x_hat);

  // P = (I - K*C)*P;
  multiply_matrix(K, n, m, C, n, temp_1);
  multiply_matrix_by_scalar(temp_1, n, n, -1, temp_2);
  add_matrix(id, n, n, temp_2, temp_1);
  multiply_matrix(temp_1, n, n, P, n, temp_2);
  copy_mat(temp_2, P, n * n);
}

//predict the state and update P
//param
//post
void predict_inline(TYPE* x_hat, 
            int n, int m,
            TYPE* A, TYPE* Q, TYPE* P,
            TYPE* x_hat_new, TYPE* A_T,
            TYPE* temp_1, TYPE* temp_2) {
  
  int i, j, k;
  int a_row, c_ind, c_row;

  // transpose_matrix(A, n, n, A_T); // do this separately since they never? change
  for (i = 0; i < n; i++) {
    a_row = n * i;
    for (j = 0; j < n; j++) {
      A_T[n * j + i] = A[a_row + j];
    }
  }

  // x_hat_new = A * x_hat
  // multiply_matrix(A, n, n, x_hat, 1, x_hat_new);
  for (i = 0; i < n; i++) {
    a_row = n * i;
    c_row = 1 * i;
    for (j = 0; j < 1; j++) {
      c_ind = j + c_row;
      x_hat_new[c_ind] = 0;
      for (k = 0; k < n; k++) {
        x_hat_new[c_ind] += A[a_row + k] * x_hat[1 * k + j];
      }
    } 
  }


  // P = A*P*A_T + Q;
  // multiply_matrix(A, n, n, P, n, temp_1);
  for (i = 0; i < n; i++) {
    a_row = n * i;
    c_row = n * i;
    for (j = 0; j < n; j++) {
      c_ind = j + c_row;
      temp_1[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_1[c_ind] += A[a_row + k] * P[n * k + j];
      }
    } 
  }
  // multiply_matrix(temp_1, n, n, A_T, n, temp_2);
  for (i = 0; i < n; i++) {
    a_row = n * i;
    c_row = n * i;
    for (j = 0; j < n; j++) {
      c_ind = j + c_row;
      temp_2[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_2[c_ind] += temp_1[a_row + k] * A_T[n * k + j];
      }
    } 
  }
  // add_matrix(temp_2, n, n, Q, P);
  for (i = 0; i < n; i++) {
    c_row = n * i;
    for (j = 0; j < n; j++) {
      c_ind = c_row + j;
      P[c_ind] = temp_2[c_ind] + Q[c_ind];
    }
  }

}



//correct the filter based on measurement
//param 
//post
//note writing this function made me want to cry
//     please send chocolate
void correct_inline(TYPE* y, TYPE* x_hat, 
            int n, int m,
            TYPE* C, TYPE* R, TYPE* P, TYPE* K,
            TYPE* x_hat_new, TYPE* C_T, TYPE* id,
            TYPE* temp_1, TYPE* temp_2, TYPE* temp_3, TYPE* temp_4) {

  // stuff that is generally useful (particularly add and mult)
  int i, j, k;
  int c_ind, a_row, c_row;

  // stuff for inverting
  TYPE cofactor[n*n];
  TYPE adjoint[n*n];
  TYPE det;

  // stuff for determinant
  int num_pivots;
  int size_a = m*m;
  TYPE L[size_a];
  TYPE U[size_a];
  TYPE P2[size_a];

  // stuff for LUP
  int ind_max, curr_row, next_row;
  int i2, j2, k2;
  TYPE max_a, abs_a, coeff;
  TYPE temp_row[m];

  // for cofactor
  TYPE det2 = 0;
  int r2, c2, row, rr;
  int n_b = m-1;
  int size_b = (m-1) * (m-1);
  int sign = 1;
  TYPE mat_b[size_b];

  // stuff for the inner determinant and LUP
  int i3, j3, k3;
  TYPE L_small[size_b];
  TYPE U_small[size_b];
  TYPE P_small[size_b];
  int i4, j4, k4;
  TYPE temp_row_small[n_b];
  int i5, j5, row5; // I relate much better to that 127 Hours movie having done this

  // stuff for other
  int row2, ind;


  // transpose_matrix(C, m, n, C_T); // do this separately since they never? change  
  for (i = 0; i < m; i++) {
    a_row = n * i;
    for (j = 0; j < n; j++) {
      C_T[m * j + i] = C[a_row + j];
    }
  }

  /************************** K = P*C_T*(C*P*C_T+R)^-1 ********************/

  // multiply_matrix(C, m, n, P, n, temp_1);
  for (i = 0; i < m; i++) {
    a_row = n * i;
    c_row = n * i;
    for (j = 0; j < n; j++) {
      c_ind = j + c_row;
      temp_1[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_1[c_ind] += C[a_row + k] * P[n * k + j];
      }
    } 
  }

  // multiply_matrix(temp_1, m, n, C_T, m, temp_2); 
  for (i = 0; i < m; i++) {
    a_row = n * i;
    c_row = m * i;
    for (j = 0; j < m; j++) {
      c_ind = j + c_row;
      temp_2[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_2[c_ind] += temp_1[a_row + k] * C_T[m * k + j];
      }
    } 
  }

  // add_matrix(temp_2, m, m, R, temp_1); 
  for (i = 0; i < m; i++) {
    c_row = m * i;
    for (j = 0; j < m; j++) {
      c_ind = c_row + j;
      temp_1[c_ind] = temp_2[c_ind] + R[c_ind];
    }
  }

  // invert_matrix(temp_1, m, temp_2); // (C*P*C_T+R)^-1
  if (m == 1) {
    temp_2[0] = 1 / temp_1[0];
  } else {
    
    // determinant: det = determinant_matrix(temp_1, m);
      det = 1.0;
        // LUP: num_pivots = compute_LUP(temp_1, L, U, P, m);
        num_pivots = 0;

        // set_identity(P2, m, m);
        for (i2 = 0; i2 < m; i2++) {
          a_row = m * i2;
          for (j2 = 0; j2 < m; j2++) {
            P2[a_row + j2] = (double)(i2 == j2);
          }
        }
        // set_identity(L, m, m);
        for (i2 = 0; i2 < m; i2++) {
          a_row = m * i2;
          for (j2 = 0; j2 < m; j2++) {
            L[a_row + j2] = (double)(i2 == j2);
          }
        }
        // copy_mat(temp_1, U, size_a);
        for (i2 = 0; i2 < size_a; i2++)
          U[i2] = temp_1[i2];

        for(i = 0; i < m; i++) {
          curr_row = i * m;
          // max_a = get_abs(U[curr_row + i]);
          max_a = (((U[curr_row + i] < 0) * -2) + 1) * U[curr_row + i];
          ind_max = i;

          for (j = i+1; j < m; j++) {
            // abs_a = get_abs(U[j * m + i]);
            abs_a = (((U[j * m + i] < 0) * -2) + 1) * U[j * m + i];
            if (abs_a > max_a) {
              max_a = abs_a;
              ind_max = j;
            }
          }

          if (ind_max != i) {
            num_pivots++;
            ind_max *= m;

            // copy_mat(&P2[curr_row], temp_row, m);
            for (i2 = 0; i2 < size_a; i2++)
              temp_row[i2] = (&P2[curr_row])[i2];
            // copy_mat(&P2[ind_max], &P2[curr_row], m);
            for (i2 = 0; i2 < size_a; i2++)
              (&P2[curr_row])[i2] = (&P2[ind_max])[i2];
            // copy_mat(temp_row, &P[ind_max], m);
            for (i2 = 0; i2 < size_a; i2++)
              (&P[ind_max])[i2] = temp_row[i2];

            // copy_mat(&U[curr_row], temp_row, m);
            for (i2 = 0; i2 < size_a; i2++)
              temp_row[i2] = (&U[curr_row])[i2];
            // copy_mat(&U[ind_max], &U[curr_row], m);
            for (i2 = 0; i2 < size_a; i2++)
              (&U[ind_max])[i2] = (&U[curr_row])[i2];
            // copy_mat(temp_row, &U[ind_max], m);
            for (i2 = 0; i2 < size_a; i2++)
              (&U[ind_max])[i2] = temp_row[i2];
          }

          for(j = i+1; j < m; j++) {
            next_row = j * m;
            coeff = (U[next_row+i]/U[curr_row+i]);
            L[next_row+i] = coeff;
            for (k = i; k < m; k++) {
              U[next_row + k] -= coeff * U[curr_row + k];
            }
          }

        } //end main for
        // end LUP
      det = (num_pivots%2) == 1 ? -1.0 : 1.0;
      for (i = 0; i < m; i++) 
        det *= U[i*m+i];
    // end determinant

    det = 1 / det;

    // cofactor_matrix(temp_1, m, cofactor);
    for (i2 = 0; i2 < m; i2++) {
      row = m * i2;
      for (j2 = 0; j2 < m; j2++) {

          k2 = 0;
          for (r2 = 0; r2 < m; r2++) {
            if(r2 != i2){
              rr = m * r2;
              for (c2 = 0; c2 < m; c2++) {
                if(c2 != j2) mat_b[k2++] = temp_1[rr + c2];
              }
            }
          }

        // det2 = determinant_matrix(mat_b, n_b);

        det2 = 1.0;
        
        // num_pivots = compute_LUP(mat_b, L_small, U_small, P_small, n_b);
          num_pivots = 0;
          // set_identity(P_small, n_b, n_b);
          for (i5 = 0; i5 < n_b; i5++) {
            row5 = n_b * i5;
            for (j5 = 0; j5 < n_b; j5++) {
              P_small[row5 + j5] = (double)(i5 == j5);
            }
          }
          // set_identity(L_small, n_b, n_b);
          for (i5 = 0; i5 < n_b; i5++) {
            row5 = n_b * i5;
            for (j5 = 0; j5 < n_b; j5++) {
              P_small[row5 + j5] = (double)(i5 == j5);
            }
          }
          // copy_mat(mat_b, U_small, size_b);
          for (i5 = 0; i5 < size_b; i5++)
            U_small[i5] = mat_b[i5];

          for(i4 = 0; i4 < n_b; i4++) {
            curr_row = i4 * n_b;
            // max_a = get_abs(U[curr_row + i4]);
            max_a = (((U[curr_row + i4] < 0) * -2) + 1) * U[curr_row + i4];
            ind_max = i4;

            for (j4 = i4+1; j4 < n_b; j4++) {
              // abs_a = get_abs(U[j4 * n_b + i4]);
              abs_a = (((U[j4 * n_b + i4] < 0) * -2) + 1) * U[j4 * n_b + i4];
              if (abs_a > max_a) {
                max_a = abs_a;
                ind_max = j4;
              }
            }

            if (ind_max != i4) {
              num_pivots++;
              ind_max *= n_b;

              // copy_mat(&P_small[curr_row], temp_row, n_b);
              for (i5 = 0; i5 < n_b; i5++)
                temp_row[i5] = (&P_small[curr_row])[i5];
              // copy_mat(&P_small[ind_max], &P_small[curr_row], n_b);
              for (i5 = 0; i5 < n_b; i5++)
                (&P_small[curr_row])[i5] = (&P_small[ind_max])[i5];
              // copy_mat(temp_row, &P_small[ind_max], n_b);
              for (i5 = 0; i5 < n_b; i5++)
                (&P_small[ind_max])[i5] = temp_row[i5];

              // copy_mat(&U_small[curr_row], temp_row, n_b);
              for (i5 = 0; i5 < n_b; i5++)
                temp_row[i5] = (&U_small[curr_row])[i5];
              // copy_mat(&U_small[ind_max], &U_small[curr_row], n_b);
              for (i5 = 0; i5 < n_b; i5++)
                (&U_small[curr_row])[i5] = (&U_small[ind_max])[i5];
              // copy_mat(temp_row, &U_small[ind_max], n_b);
              for (i5 = 0; i5 < n_b; i5++)
                (&U_small[ind_max])[i5] = temp_row[i5];
            }

            for(j4 = i4+1; j4 < n_b; j4++) {
              next_row = j4 * n_b;
              coeff = (U[next_row+i4]/U[curr_row+i4]);
              L_small[next_row+i4] = coeff;
              for (k4 = i4; k4 < n_b; k4++) {
                U_small[next_row + k4] -= coeff * U_small[curr_row + k4];
              }
            }

          } //end main for
        // end inner LUP

        det2 = (num_pivots%2) == 1 ? -1.0 : 1.0;

        for (i3 = 0; i3 < n_b; i3++) {
          det2 *= U_small[i3*n_b+i3];
        }
        // end det2

        cofactor[row + j2] = sign * det2;
        sign = sign * -1;
      }
      sign = sign * -1;
    }
    // end cofactor

    // transpose_matrix(cofactor, m, m, adjoint);
    for (i5 = 0; i5 < m; i5++) {
      row5 = m * i5;
      for (j5 = 0; j5 < m; j5++) {
        adjoint[m * j5 + i5] = cofactor[row5 + j5];
      }
    }

    // multiply_matrix_by_scalar(adjoint, m, m, det, temp_2);  
    for (i2 = 0; i2 < m; i2++) {
      row = m * i2;
      for (j2 = 0; j2 < m; j2++) {
        ind = row + j2;
        temp_2[ind] = adjoint[ind] * det;
      }
    }


  }
  // end invert

  // multiply_matrix(P, n, n, C_T, m, temp_1); // P*C_T
  for (i = 0; i < n; i++) {
    a_row = n * i;
    c_row = m * i;
    for (j = 0; j < m; j++) {
      c_ind = j + c_row;
      temp_1[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_1[c_ind] += P[a_row + k] * C_T[m * k + j];
      }
    } 
  }

  // multiply_matrix(temp_1, n, m, temp_2, m, K);
  for (i = 0; i < n; i++) {
    a_row = m * i;
    c_row = m * i;
    for (j = 0; j < m; j++) {
      c_ind = j + c_row;
      K[c_ind] = 0;
      for (k = 0; k < m; k++) {
        K[c_ind] += temp_1[a_row + k] * temp_2[m * k + j];
      }
    } 
  }


  /************************** x_hat = x_hat_new + K * (y - C*x_hat_new); **************************/
  // multiply_matrix(C, m, n, x_hat_new, 1, temp_3);
  for (i = 0; i < m; i++) {
    a_row = n * i;
    c_row = 1 * i;
    for (j = 0; j < 1; j++) {
      c_ind = j + c_row;
      temp_3[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_3[c_ind] += C[a_row + k] * x_hat_new[1 * k + j];
      }
    } 
  }

  // multiply_matrix_by_scalar(temp_3, m, 1, -1, temp_4);
  for (i = 0; i < m; i++) {
    row = 1 * i;
    for (j = 0; j < 1; j++) {
      ind = row + j;
      temp_4[ind] = temp_3[ind] * -1;
    }
  }


  // add_matrix(y, m, 1, temp_4, temp_3);
  for (i = 0; i < m; i++) {
    row = 1 * i;
    for (j = 0; j < 1; j++) {
      ind = row + j;
      temp_3[ind] = y[ind] + temp_4[ind];
    }
  }

  // multiply_matrix(K, n, m, temp_3, 1, temp_4);
  for (i = 0; i < n; i++) {
    a_row = m * i;
    c_row = 1 * i;
    for (j = 0; j < 1; j++) {
      c_ind = j + c_row;
      temp_4[c_ind] = 0;
      for (k = 0; k < m; k++) {
        temp_4[c_ind] += K[a_row + k] * temp_3[1 * k + j];
      }
    } 
  }

  // add_matrix(x_hat_new, n, 1, temp_4, x_hat);
  for (i = 0; i < n; i++) {
    row = 1 * i;
    for (j = 0; j < 1; j++) {
      ind = row + j;
      x_hat[ind] = x_hat_new[ind] + temp_4[ind];
    }
  }


  /************************** P = (I - K*C)*P; **************************/
  // multiply_matrix(K, n, m, C, n, temp_1);
  for (i = 0; i < n; i++) {
    a_row = m * i;
    c_row = n * i;
    for (j = 0; j < n; j++) {
      c_ind = j + c_row;
      temp_1[c_ind] = 0;
      for (k = 0; k < m; k++) {
        temp_1[c_ind] += K[a_row + k] * C[n * k + j];
      }
    } 
  }


  // multiply_matrix_by_scalar(temp_1, n, n, -1, temp_2);
  for (i = 0; i < n; i++) {
    row = n * i;
    for (j = 0; j < n; j++) {
      ind = row + j;
      temp_2[ind] = temp_1[ind] * -1;
    }
  }

  // add_matrix(id, n, n, temp_2, temp_1);
  for (i = 0; i < n; i++) {
    row = n * i;
    for (j = 0; j < n; j++) {
      ind = row + j;
      temp_1[ind] = id[ind] + temp_2[ind];
    }
  }

  // multiply_matrix(temp_1, n, n, P, n, temp_2);
  for (i = 0; i < n; i++) {
    a_row = n * i;
    c_row = n * i;
    for (j = 0; j < n; j++) {
      c_ind = j + c_row;
      temp_2[c_ind] = 0;
      for (k = 0; k < n; k++) {
        temp_2[c_ind] += temp_1[a_row + k] * P[n * k + j];
      }
    } 
  }

  // copy_mat(temp_2, P, n * n);
  size_a = n*n;
  for (i = 0; i < size_a; i++)
    P[i] = temp_2[i];
}
