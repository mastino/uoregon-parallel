/* kalman_tracker.h
 * a class that uses kalman filters to track the objects
 
 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * If by some miracle you find this software useful, thanks are accepted in
 * the form of chocolate or introductions to potential employers.

 */



#ifndef KALMAN_TRACKER_H
#define KALMAN_TRACKER_H

#include <opencv/cv.hpp>
#include "opencv2/core.hpp"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "object.h"
#include "useful_functions.h"
#include "../kalman-filters/basic-c/kalman_filter.h"
#include "../kalman-filters/basic-c/linear_algebra.h"

using namespace std;

// for matrix setting
#define KT_A 1
#define KT_C 2
#define KT_Q 3
#define KT_R 4
#define KT_P 5
#define KT_K 6

#define KT_x 7
#define KT_y 8
#define KT_x_hat 9


class Kalman_Tracker {

public:
  Kalman_Tracker();
  ~Kalman_Tracker();

  //@set the values in the matrix
  //@pre allocate memory for the matrix using the constructor
  //@post deep copy of init_mat into the selected mat
  //@param init_mat pointer to 1D array with data to copy
  //       mat_name the name of the mat using the above macros
  void set_matrix(TYPE init_mat, int mat_name);

  //@get the values in the matrix
  //@post shallow copy of selected mat to mat_ptr (i.e. don't change it!)
  //@param mat_name the name of the mat using the above macros (input)
  //       mat_ptr - ptr for matrix - should not be pointing anywhere yet (output)
  //       r, c - num rows and columns in the selected mat (outputs)
  bool Kalman_Tracker::get_mat(int mat_name, TYPE* &mat_ptr, int &r, int &c);
  
  void update();
  void predict();
  void correct();

private:
  // arrays
  TYPE *A, 
       *C, 
       *Q, 
       *R, 
       *P, 
       *K, 
       *x, 
       *y, 
       *x_hat,
       *x_hat_new, 
       *A_T, 
       *C_T, 
       *id,
       *temp_1, 
       *temp_2, 
       *temp_3, 
       *temp_4;

  int n; // Number of states
  int m; // Number of measurements

  double dt; // Time step
  double t; 

};

#endif