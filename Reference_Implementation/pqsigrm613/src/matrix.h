#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define MATRIX_NULL 0 
#define ELEMBLOCKSIZE 8

#define INV_SUCCESS 1
#define INV_FAIL 0

#define getElement(A, i, j) 		(!!((A)->elem[(i) * (A)->rwdcnt + (j) / ELEMBLOCKSIZE] & (0x80 >> ((j) % ELEMBLOCKSIZE))))
#define flipElement(A, i, j) 	((A)->elem[(i) * (A)->rwdcnt + (j) / ELEMBLOCKSIZE] ^= (0x80 >> ((j) % ELEMBLOCKSIZE)))
#define setElement(A, i, j, val) 	((getElement((A), (i), (j)) == (val))? 0 : flipElement((A), (i),(j)))
#define initZero(R) 			memset((R)->elem,0,(R)->alloc_size)

typedef struct {
   int nrows;//number of rows.
   int ncols;//number of columns.
   int words_in_row;//number of words in a row
   int alloc_size;//number of allocated bytes
   unsigned char *elem;//row index.
} matrix;

matrix* new_matrix(int nrows, int ncols) ;
void delete_matrix(matrix * mtx) ;

matrix* rref(matrix* mtx);
matrix* transpose(matrix *dest, matrix *src);
int inverse(matrix *mtx, matrix *mtxInv);
int is_nonsingular(matrix *mtx);

void get_pivot(matrix* mtx, uint16_t *lead, uint16_t *lead_diff);

matrix* copy_matrix(matrix* dest, matrix* src);

int mat_mat_prod(matrix * mtx1, matrix * mtx2, matrix * prod); 
void vec_mat_prod(matrix *dest, matrix* m, matrix *vec);
int mat_mat_add(matrix *m1, matrix *m2, matrix *res);

int export_matrix(unsigned char* dest, matrix* mtx);
matrix* import_matrix(matrix* dest_mtx, const unsigned char* src);

void dual(matrix* G, matrix* H_sys, uint16_t *lead, uint16_t *lead_diff);
void row_interchange(matrix* mtx, int row_idx1, int row_idx2);
void partial_replace(matrix* dest, const int r1, const int c1,const int r2, const int c2, matrix* src, const int r3, const int c3);
#endif
