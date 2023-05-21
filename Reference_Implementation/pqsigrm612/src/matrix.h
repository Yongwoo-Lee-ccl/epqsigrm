#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define MATRIX_NULL 0 

#define INV_SUCCESS 1
#define INV_FAIL 0

#define get_element(A, i, j)  ((A->elem)[(i)][(j)])
#define set_element(A, i, j, val) ((A->elem)[(i)][(j)] = (val))

typedef struct {
    uint16_t nrows;//number of rows.
    uint16_t ncols;//number of columns.
    uint8_t** elem;//row index.
} matrix;

matrix* new_matrix(uint32_t nrows, uint32_t ncols) ;
void init_zero(matrix *self);
void delete_matrix(matrix *self) ;
void randomize(matrix *self, uint8_t* randstr);

matrix* copy_matrix(matrix* self, matrix* src);
void export_matrix(matrix* self, uint8_t* dest);
void import_matrix(matrix* self, const uint8_t* src);

matrix* rref(matrix* self);
matrix* transpose(matrix *self, matrix* dest);
int inverse(matrix *self, matrix *dest);
int is_nonsingular(matrix *self);

void get_pivot(matrix* self, uint16_t *lead, uint16_t *lead_diff);

void mat_mat_prod(matrix* self, matrix* mat1, matrix* mat2); 
void vec_mat_prod(matrix* self, matrix* mat,  matrix *vec);
void mat_mat_add(matrix* self, matrix* mat1, matrix *mat2);

void dual(matrix* self, matrix* dual_sys);
void row_interchange(matrix* mtx, uint32_t row_idx1, uint32_t row_idx2);
void partial_replace(matrix* self, const uint32_t r1, const uint32_t c1,
        const uint32_t r2, const uint32_t c2, 
        matrix* src, const int r3, const int c3);

void codeword(matrix* self, uint8_t* seed, matrix* dest);

#endif
