#include "matrix.h"

matrix* new_matrix (int nrows, int ncols)
{
  matrix* A;

  A = (matrix*) malloc (sizeof (matrix));
  A->ncols = ncols;
  A->nrows = nrows;  
  A->words_in_row = (1 + (ncols - 1) / ELEMBLOCKSIZE);
  A->alloc_size = nrows * A->words_in_row * sizeof (unsigned char);
  A->elem = (unsigned char *)malloc(A->alloc_size);
  init_zero(A);
  return A;
}

void partial_replace(matrix* dest, const int r1, const int c1,const int r2, const int c2, matrix* src, const int r3, const int c3){
	for(int i = 0; i < r2 - r1; i++)
		for(int j = 0; j < c2 - c1; j++)
			set_element(dest, r1 + i, c1+j, get_element(src, r3 + i, c3 + j));
}

void delete_matrix(matrix* A)
{
  free(A->elem);
  free(A);
}

int mat_mat_prod(matrix * mtx1, matrix * mtx2, matrix * prod) {
	int row, col, k;
	int val;
	if(mtx1->ncols != mtx2->nrows) 	return -1;
	
	for (row = 0; row < mtx1->nrows; row++)
		for (col = 0; col < mtx2->ncols; col++) {
			val = 0;
			for (k = 0; k < mtx1->ncols; k++)
				val ^= get_element(mtx1, row, k) & get_element(mtx2, k, col);
			set_element(prod, row, col, val);
		}
	return 0;
}

//assume vector is transposed
void vec_mat_prod(matrix *dest, matrix* m, matrix *vec){
	unsigned char bit = 0;
	unsigned char offset;
	int row, col;
	for(row = 0; row < m->nrows; row++){
		bit = 0;
		for(col=0; col < m->words_in_row - 1; col++)
			bit ^= m->elem[row*m->words_in_row + col] & vec->elem[col];
	
		offset = 0xff << (ELEMBLOCKSIZE*m->words_in_row - m->ncols);
		bit ^= (m->elem[row*m->words_in_row + col] & vec->elem[col])&offset;

		bit ^= (bit >> 4);
		bit ^= (bit >> 2);
		bit ^= (bit >> 1);
		bit &= (unsigned char)1;
		
		set_element(dest, 0, row, bit);
	}
}

void row_interchange(matrix* A, int row_idx1, int row_idx2){
	int col_idx;
	unsigned char temp;
	for(col_idx=0; col_idx<A->words_in_row; ++col_idx){
		temp 	 								= A->elem[row_idx1 * A->words_in_row + col_idx];
		A->elem[row_idx1 * A->words_in_row + col_idx] = A->elem[row_idx2 * A->words_in_row + col_idx];
		A->elem[row_idx2 * A->words_in_row + col_idx] = temp;
	}
}

void row_add(matrix* A, int dest_row_idx, int adding_row_idx){
	int col_idx;
	for(col_idx=0; col_idx<A->words_in_row; ++col_idx){
		A->elem[dest_row_idx * A->words_in_row + col_idx] 
			^= A->elem[adding_row_idx * A->words_in_row + col_idx];
	}
}

matrix * rref(matrix* A)
{
	// Considering column is longer than row
	int succ_row_idx=0;
	int col_idx, row_idx=0;
	int i;
	for (col_idx = 0; col_idx < (A->ncols); ++col_idx) {
		
		// finding first row s.t. i th elem of the row is 1
		for(; row_idx < A->nrows; ++row_idx)
			if(get_element(A, row_idx, col_idx) == 1) 
				break;
		// When reaches the last row,
		// increase column index and search again
		if (row_idx == A->nrows){ 
			row_idx=succ_row_idx;
			continue;
		}
		// if row_idx is not succ_row_idx, 
		// interchange between:
		// <succ_row_idx> th row <-> <row_idx> th row
		if(row_idx != succ_row_idx){
			row_interchange(A, succ_row_idx, row_idx);
		}
				
		// By adding <succ_row_idx> th row in the other nrows 
		// s.t. A(i, <succ_row_idx>) == 1,
		// making previous columns as element row.
		for(i=0; i<A->nrows; ++i){
			if(i == succ_row_idx) continue;

			if(get_element(A, i, col_idx) == 1){
				row_add(A, i, succ_row_idx);
			}
		}
		row_idx = ++succ_row_idx;
	}
	//Gaussian elimination is finished. So return A.
	return A;
}

int inverse(matrix *mtx, matrix *mtxInv){
	if(mtx->nrows != mtx->ncols) 			return INV_FAIL;
	if(mtxInv->nrows != mtxInv->ncols) 	return INV_FAIL;
	if(mtx->nrows != mtxInv->nrows) 		return INV_FAIL;

	matrix* temp = new_matrix(mtx->nrows, mtx->ncols);
	copy_matrix(temp, mtx);

	int r, c;
	for(r = 0; r< mtxInv->alloc_size;++r){
		mtxInv->elem[r] = 0;
	}
	for ( r = 0; r <  mtxInv->nrows; ++r)
	{
		set_element(mtxInv, r, r, 1);
	}

	for (c = 0; c < temp->ncols; ++c)
	{
		if(get_element(temp, c, c) == 0)
		{	
			for (r = c+1; r < mtx->nrows; ++r)
			{
				if(get_element(temp, r, c) != 0){
					row_interchange(temp, r, c);
					row_interchange(mtxInv, r, c);
					break;
				}
			}
			if(r >= temp->nrows) 		return INV_FAIL;
		}
		

		for(r = 0; r < temp->nrows; r++){
			if(r == c) continue;
			if(get_element(temp, r, c) != 0){
				row_add(temp, r, c);
				row_add(mtxInv, r, c);
			}
		}
	}
	delete_matrix(temp);
	return INV_SUCCESS;
}

int is_nonsingular(matrix *mtx){

	matrix* temp = new_matrix(mtx->nrows, mtx->ncols);
	copy_matrix(temp, mtx);

	int r, c;
	unsigned char bit_one =	0x80;

	for (c = 0; c < temp->ncols; ++c)
	{
		if(get_element(temp, c, c) == 0)
		{	
			for (r = c+1; r < mtx->nrows; ++r)
			{
				if(get_element(temp, r, c) != 0){
					row_interchange(temp, r, c);
					break;
				}
			}
			if(r >= temp->nrows) 		return INV_FAIL;
		}
		

		for(r = 0; r < temp->nrows; r++){
			if(r == c) continue;
			if(get_element(temp, r, c) != 0){
				row_add(temp, r, c);
			}
		}
	}
	delete_matrix(temp);
	return INV_SUCCESS;
}


matrix* copy_matrix(matrix* dest, matrix* src){
	if(dest->nrows != src->nrows || dest->ncols!=src->ncols) return MATRIX_NULL;
	
	memcpy(dest->elem, src->elem, dest->alloc_size);
	return dest;
}

matrix* transpose(matrix *dest, matrix *src){
	if((dest->nrows != src->ncols) || (dest->ncols != src->nrows))
		return MATRIX_NULL;
	int row, col;
	for(row=0; row < dest->nrows; ++row)
		for(col=0; col < dest->ncols; ++col)
			set_element(dest, row, col, get_element(src, col, row));
	return dest;
}

// Exports a matrix into unsigned char destination.
int export_matrix(unsigned char* dest, matrix* src_mtx){
	memcpy(dest, src_mtx->elem, src_mtx->alloc_size);
	return src_mtx->alloc_size;
}

matrix* import_matrix(matrix* dest_mtx, const unsigned char* src){
	memcpy(dest_mtx->elem, src, dest_mtx->alloc_size);

	return dest_mtx;
}

int mat_mat_add(matrix *m1, matrix *m2, matrix *res){
	if((m1->nrows != m2->nrows) || (m1->ncols != m2->ncols))
		return -1;
	
	for(int i =0; i< res->alloc_size; i++)
		res->elem[i] = (m1->elem[i])^(m2->elem[i]);

	return 0;
}

void get_pivot(matrix* mtx, uint16_t *lead, uint16_t *lead_diff){
	int row=0, col=0;
	int lead_idx=0, diff_idx=0;
	while((col < mtx->ncols) && (row < mtx->nrows) && (lead_idx < mtx->nrows) && (diff_idx < (mtx->ncols - mtx->nrows))){
		if(get_element(mtx, row, col) == 1){
			lead[lead_idx++] = col;
			row++;
		}
		else{
			lead_diff[diff_idx++] =col;
		}
		col++;
	}

	while(col < mtx->ncols){
		lead_diff[diff_idx++]=col++;
	}
}

void dual(matrix* G, matrix* H_sys, uint16_t *lead, uint16_t *lead_diff){
	int row, col, flg = 0; 
	rref(G);
	if(lead == 0 || lead_diff == 0){
		lead = (uint16_t*)malloc(sizeof(uint16_t)*G->nrows);
		lead_diff = (uint16_t*)malloc(sizeof(uint16_t)*(G->ncols - G->nrows));	
		flg = 1;
	}
	get_pivot(G, lead, lead_diff);
	// Fill not-identity part (P')
	for ( row = 0; row < H_sys->nrows; row++) 
		for ( col = 0; col < G->nrows; col++) 
			set_element(H_sys, row, lead[col], get_element(G, col, lead_diff[row]));

	for ( row = 0; row < H_sys->nrows; row++) 
			set_element(H_sys, row, lead_diff[row], 1);
	
	if(flg){
		free(lead);
		free(lead_diff);
	}
}
