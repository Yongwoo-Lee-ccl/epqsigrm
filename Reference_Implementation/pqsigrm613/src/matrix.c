#include "matrix.h"

matrix* new_matrix (uint32_t nrows, uint32_t ncols)
{
    matrix* mat;
    mat = (matrix*) malloc (sizeof (matrix));
    mat->nrows = nrows; 
    mat->ncols = ncols;
    // Pad zeros up to colsize * 64 to fit with opt. implementation
    mat->colsize  = (ncols + 63)/64;
    mat->elem = (uint8_t**)malloc(nrows * sizeof(uint8_t*));
    for (uint32_t i = 0; i < nrows; i++)
    {
        mat->elem[i] = (uint8_t*)calloc(mat->colsize * 64, sizeof(uint8_t));
    }
    init_zero(mat);
    return mat;
}


void init_zero(matrix *self){
    for (uint32_t i = 0; i < self->nrows; i++)
    {
        for (uint32_t j = 0; j < self->ncols; j++)
        {
            set_element(self, i, j, 0);
        }
    }
}

// generate a random matrix using a random string.
// random string contains (nrows * ncols)-bit, while ncols is multiple of 8 
void randomize(matrix *self, uint8_t* randstr){
    uint32_t bit_index = 0;
    for (uint32_t i = 0; i < self->nrows; i++){
        for (uint32_t j = 0; j < self->ncols; j++){
            uint32_t byte_index = bit_index >> 3; // byte_index = bit_index / 8;
            uint32_t bit_offset = bit_index & 7; // bit_offset = bit_index % 8;
            uint8_t rand_bit = (randstr[byte_index] >> bit_offset) & 1;
            set_element(self, i, j, rand_bit);
            bit_index++;
        }
    }
}

void delete_matrix(matrix* self)
{
    for (uint32_t i = 0; i < self->nrows; i++) {
        free(self->elem[i]);
    }
    free(self->elem);
    free(self);
}

void copy_matrix(matrix* self, matrix* src){
    for (uint32_t i = 0; i < src->nrows; i++)
    {
        memcpy(self->elem[i], src->elem[i], src->ncols);
    }
}

void export_matrix(matrix* self, uint8_t* dest){
    uint32_t byte_index = 0;
    for (uint32_t i = 0; i < self->nrows; i++){
        for (uint32_t j = 0; j < self->colsize * 64; j += 8){
            uint8_t byte = 0;
            for (uint32_t k = 0; (k < 8) && (j+k < self->ncols); k++)
            {
                byte |= (self->elem[i][j+k] << k);
            }
            dest[byte_index++] = byte; 
        }
    }
}

void import_matrix(matrix* self, const uint8_t* src){
    uint32_t byte_index = 0;
    for (uint32_t i = 0; i < self->nrows; i++){
        for (uint32_t j = 0; j < self->colsize * 64; j+=8){
            uint8_t byte = src[byte_index++];
            for (uint32_t k = 0; (k < 8) && (j+k < self->ncols); k++)
            {
                self->elem[i][j+k] = ( byte >> k) & 1;
            }
        }
    }
}

void row_addition_internal(matrix* self, const uint32_t r1, const uint32_t r2){
    for (uint32_t j = 0; j < self->ncols; j++) {
        uint8_t bit = get_element(self, r1, j) ^ get_element(self, r2, j);
        set_element(self, r1, j, bit);
    }
}

// make a matrix into rref form, inplace
void rref(matrix* self)
{
    // Assume column is longer than row
    uint32_t succ_row_idx = 0;
    uint32_t col_idx, row_idx = 0;
    for (col_idx = 0; col_idx < (self->ncols); ++col_idx) {
        
        // finding first row s.t. i th elem of the row is 1
        for(; row_idx < self->nrows; ++row_idx)
            if(get_element(self, row_idx, col_idx) == 1) 
                break;
        // When reaches the last row, increase column index and search again
        if (row_idx == self->nrows){ 
            row_idx = succ_row_idx;
            continue;
        }
        // if row_idx is not succ_row_idx, 
        // interchange between:
        // <succ_row_idx> th row <-> <row_idx> th row
        if(row_idx != succ_row_idx){
            row_interchange(self, succ_row_idx, row_idx);
        }
                
        // By adding <succ_row_idx> th row in the other nrows 
        // s.t. self(i, <succ_row_idx>) == 1,
        // making previous columns as element row.
        for(uint32_t i = 0; i < self->nrows; ++i){
            if(i == succ_row_idx) continue;

            if(get_element(self, i, col_idx) == 1){
                row_addition_internal(self, i, succ_row_idx);
            }
        }
        row_idx = ++succ_row_idx;
    }
}

// Input should be in rref form.
void get_pivot(matrix* self, uint16_t* lead, uint16_t* lead_diff){
    uint16_t row=0, col=0;
    uint16_t lead_idx=0, diff_idx=0;

    while((col < self->ncols) 
            && (row < self->nrows) 
            && (lead_idx < self->nrows) 
            && (diff_idx < (self->ncols - self->nrows))){

        if(get_element(self, row, col) == (uint8_t)1){
            lead[lead_idx++] = col++;
            row++;
        }
        else{
            lead_diff[diff_idx++] = col++;
        }
    }

    while(col < self->ncols){
        if(lead_idx < self->nrows) {
            lead[lead_idx++] = col++;
        }
        else{
            lead_diff[diff_idx++] = col++;
        }
    }
}

// assume vector is transposed
// self is also transposed
void vec_mat_prod(matrix* self, matrix* mat, matrix* vec){

    for(uint32_t i = 0; i < mat->nrows; i++) {
        uint8_t bit = 0;
        for (uint32_t j = 0; j < mat->ncols; j++) {
            bit ^= get_element(mat, i, j) & get_element(vec, 0, j);
        }
        set_element(self, 0, i, bit);
    }
}

void vec_vec_add(matrix* self, matrix* vec){
    for (uint32_t j = 0; j < self->ncols; j++) {
        uint8_t bit = get_element(self, 0, j) ^ get_element(vec, 0, j);
        set_element(self, 0, j, bit);
    }
}

void dual(matrix* self, matrix* dual_sys){
    uint16_t lead[self->nrows];
    uint16_t lead_diff[self->ncols - self->nrows];    

    init_zero(dual_sys);

    rref(self);
    get_pivot(self, lead, lead_diff);

    // Fill not-identity part (P')
    for (uint32_t row = 0; row < dual_sys->nrows; row++) {
        for (uint32_t col = 0; col < self->nrows; col++) {
            set_element(dual_sys, row, lead[col], get_element(self, col, lead_diff[row]));
        }
    }
    
    for (uint32_t row = 0; row < dual_sys->nrows; row++) {
        set_element(dual_sys, row, lead_diff[row], 1);    
    }
}

void row_interchange(matrix* self, uint32_t row1, uint32_t row2) {
    uint8_t* temp = self->elem[row1];
    self->elem[row1] = self->elem[row2];
    self->elem[row2] = temp;
}

void partial_replace(matrix* self, const uint32_t r1, const uint32_t r2,
        const uint32_t c1, const uint32_t c2, 
        matrix* src, const int r3, const int c3){
    for(uint32_t i = 0; i < r2 - r1; i++) {
        for(uint32_t j = 0; j < c2 - c1; j++) {
            set_element(self, r1 + i, c1+j, get_element(src, r3 + i, c3 + j));
        }
    }
}

void codeword(matrix* self, uint8_t* seed, matrix* dest){
    for (uint32_t i = 0; i < self->nrows; i++) {
        uint32_t byte_index = i >> 3; // byte_index = bit_index / 8;
        uint32_t bit_offset = i & 7; // bit_offset = bit_index % 8;
        uint8_t rand_bit = (seed[byte_index] >> bit_offset) & 1;
        if (rand_bit == 1) {
            for (uint32_t j = 0; j < self->ncols; j++)
            {
                dest->elem[0][j] ^= self->elem[i][j];
            }
        }
    }
}

uint8_t is_zero(matrix* self){
    for (uint32_t i = 0; i < self->nrows; i++) {
        for (size_t j = 0; j < self->ncols; j++) {
            if(get_element(self, i, j) != 0) {
                return 0;
            }
        }
    }
    return 1;
}

void col_permute(matrix* self, const int r1, const int r2
	, const int c1, const int c2, uint16_t* Q) {	
	matrix* copy = new_matrix(r2 - r1, c2 - c1);
	for (uint32_t r = 0; r < r2 - r1; r++) {
		for (uint32_t c = 0; c < c2 - c1; c++) {
			uint8_t bit = get_element(self, r1 + r, c1 + c);
			set_element(copy, r, c, bit);
		}
	}
	
	for(uint32_t c = 0; c < c2 - c1; c++) {
		for(uint32_t r = 0; r < r2 - r1; r++) {
            uint8_t bit =  get_element(copy, r, Q[c]);
			set_element(self, r1 + r, c1 + c, bit);
		}
	}

	delete_matrix(copy);
}
