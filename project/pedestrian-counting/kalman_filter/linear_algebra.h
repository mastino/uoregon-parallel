/*
 * linear_algebra.h
 * simple header for matrix function code
 
 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * credit to www.mathwords.com for definitions and help with linear algebra

 * See LICENSE file for licensing information and boring legal stuff

 * If by some miricale you find this software useful, thanks are accepted in
 * the form of chocolate, coffee, or introductions to potential employers.

*/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <stdio.h>


#define TYPE double // in retrospect this will probably always be double

void print_matrix(TYPE* mat_a, int rows_a, int cols_a);
 
//get inverse matrix
//pre matrix a has been created and filled with data
//     matrices must be square and the same size
//params mat_a - matrix to invert
//        n -  the matrix side size
//post mat_c has the inverse matrix of mat_a
//TODO clean up and use a better algo
void invert_matrix(TYPE* mat_a, int n, TYPE* mat_c);

//get matrix determinant recursively
//pre matrix a has been created and filled with data
//     matrix must be square
//params mat_a - matrix to invert
//        n -  the matrix side size
//returms determinant of matrix a
//TODO use better algorithm or make this not recursive
TYPE determinant_matrix_recur(TYPE* mat_a, int n);

//get matrix determinant iteratively (LU factorization)
//pre matrix a has been created and filled with data
//     matrix must be square
//params mat_a - matrix to invert
//        n -  the matrix side size
//returms determinant of matrix a
//TODO make this dolittle LU factoization
TYPE determinant_matrix(TYPE* mat_a, int n);

//get cofactor matrix
//pre matrix a has been created and filled with data
//     matrices must be square and the smae size
//params mat_a - matrix to cofactor
//        n -  the matrix side size
//post mat_c has the cofactor matrix of mat_a
//TODO clean up
void cofactor_matrix(TYPE* mat_a, int n, TYPE* mat_c);

//add matrices together
//pre all matrices are initialized, c shouldn't have any important data in it
//     all matrices should be the same dimensions
//post mat_c has the result of multipling mat_a and mat_b
void add_matrix(TYPE* mat_a, int rows, int cols, TYPE* mat_b, TYPE* mat_c);

//multiply matrices together
//pre all matrices are initialized, c shouldn't have any important data in it
//     all matrices should be the same dimensions
//post mat_c has the result of multipling mat_a and scalar
void multiply_matrix_by_scalar(TYPE* mat_a, int rows, int cols, TYPE scalar, TYPE* mat_c);

//multiply matrices together
//pre all matrices are initialized, c shouldn't have any important data in it
//     rows in b == cols in a
//     c is initialized to the same size as b
//post mat_c has the result of multipling mat_a and mat_b
void multiply_matrix(TYPE* mat_a, int rows_a, int cols_a, TYPE* mat_b, int cols_b, TYPE* mat_c);

//transpose a matrix
//pre all matrices are initialized, c shouldn't have any important data in it
//     rows in c == cols in a
//post mat_c has the transpose of mat_a
void transpose_matrix(TYPE* mat_a, int rows_a, int cols_a, TYPE* mat_c);

//compute LUP factorization of the matrix
//pre all matrices are initialized, L, U, P shouldn't have anything in them
//    all are nxn
//output L is lower, U is upper, P is the P matrix which i don't understand
//returns number of row swaps or -1 if failure
int compute_LUP(TYPE* mat_a, TYPE* L, TYPE* U, TYPE* P, int n);
int compute_LUP_inline(TYPE* mat_a, TYPE* L, TYPE* U, TYPE* P, int n);

//set a matrix to zero
//pre matrix_a has been allocated to rows_a x cols_a
//post mat_a is all zeros
void set_zero(TYPE* mat_a, int rows_a, int cols_a);

//set a matrix to the identity
//pre matrix_a has been allocated to rows_a x cols_a
//post mat_a has ones in the diagonal and zeros elsewhere
void set_identity(TYPE* mat_a, int rows_a, int cols_a);

//deep copy of a to c
void copy_mat(TYPE* mat_a, TYPE* mat_c, int total_elms);

//returns abs(a)
TYPE get_abs(TYPE a);



#endif
