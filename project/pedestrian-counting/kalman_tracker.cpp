/* kalman_tracker.cpp
 * a class that uses kalman filters to track the objects
 
 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * If by some miracle you find this software useful, thanks are accepted in
 * the form of chocolate or introductions to potential employers.

 */



#ifndef KALMAN_TRACKER_H
#define KALMAN_TRACKER_H

#include "kalman_tracker.h"

using namespace std;

Kalman_Tracker::Kalman_Tracker() {
  allocate_matrices(&A, &C, &Q, &R, &P, &K, n, m);
  allocate_vectors(&x, &y, &x_hat, n, m);
  allocate_temp_matrices(&x_hat_new, &A_T, &C_T, &id,
                         &temp_1, &temp_2, &temp_3, &temp_4, n, m);
}


Kalman_Tracker::~Kalman_Tracker() {
  destroy_matrices(A, C, Q, R, P, K);
  destroy_vectors(x, y, x_hat);
  destroy_temp_matrices(x_hat_new, A_T, C_T, id,
                        temp_1, temp_2, temp_3, temp_4);
}

void Kalman_Tracker::set_matrix(TYPE init_mat, int mat_name) {

  int num_elms = 0;
  TYPE *mat_ptr = NULL;

  switch(mat_name) {
    case KT_A:
      num_elms = n*n;
      mat_ptr = A;
      break;
    case KT_C:
      num_elms = n*m;
      mat_ptr = C;
      break;
    case KT_Q:
      num_elms = n*n;
      mat_ptr = Q;
      break;
    case KT_R:
      num_elms = m*m;
      mat_ptr = R;
      break;
    case KT_P:
      num_elms = n*n;
      mat_ptr = P;
      break;
    case KT_K:
      num_elms = n*m;
      mat_ptr = K;
      break;
    case KT_x:
      num_elms = n;
      mat_ptr = x;
      break;
    case KT_y:
      num_elms = n;
      mat_ptr = y;
      break;
    case KT_x_hat:
      num_elms = n;
      mat_ptr = x_hat;
      break;
  }

  copy_matrix(init_mat, mat_ptr, num_elms);

}

bool Kalman_Tracker::get_mat(int mat_name, TYPE* &mat_ptr, int &r, int &c){

  bool out = true;

  switch(mat_name) {
    case KT_A:
      r = n;
      c = n;
      mat_ptr = A;
      break;
    case KT_C:
      r = m;
      c = n;
      mat_ptr = C;
      break;
    case KT_Q:
      r = n;
      c = n;
      mat_ptr = Q;
      break;
    case KT_R:
      r = m;
      c = m;
      mat_ptr = R;
      break;
    case KT_P:
      r = n;
      c = n;
      mat_ptr = P;
      break;
    case KT_K:
      r = n;
      c = m;
      mat_ptr = K;
      break;
    case KT_x:
      r = n;
      c = 1;
      mat_ptr = x;
      break;
    case KT_y:
      r = m;
      c = 1;
      mat_ptr = y;
      break;
    case KT_x_hat:
      r = n;
      c = 1;
      mat_ptr = x_hat;
      break;
    default:
      out = false;
  }

  return out;

}

void Kalman_Tracker::update(){
    update(y, x_hat, &t, dt, n, m, A,  C,  Q,  R,  P,  K,
           x_hat_new, A_T, C_T, id, temp_1, temp_2, temp_3, temp_4);
}

void Kalman_Tracker::predict(){
  predict(x_hat, n, m, A, Q, P, 
          x_hat_new, A_T, temp_1, temp_2);
}

void Kalman_Tracker::correct(){
  correct(y, x_hat, n, m, C, R, P, K,
          x_hat_new, C_T, id, temp_1, temp_2, temp_3, temp_4);
}


/* class variables
  TYPE *A, *C, *Q, *R, *P, *K, *x, *y, *x_hat,
       *x_hat_new, *A_T, *C_T, *id,
       *temp_1, *temp_2, *temp_3, *temp_4;

  int n; // Number of states
  int m; // Number of measurements

  double dt; // Time step
  double t; 
*/
