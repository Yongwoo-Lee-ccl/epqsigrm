#include "matrix.h"

matrix* new_matrix (uint32_t nrows, uint32_t ncols)
{
    matrix* mat;
    mat = (matrix*) malloc (sizeof (matrix));
    mat->nrows = nrows; 
    mat->ncols = ncols;
    mat->elem = (uint8_t**)malloc(nrows * sizeof(uint8_t*));
    for (uint32_t i = 0; i < nrows; i++)
    {
        mat->elem[i] = (uint8_t*)malloc(ncols * sizeof(uint8_t));
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
    for (uint32_t i = 0; i < self->nrows; i++)
    {
        free(self->elem[i]);
    }
    free(self->elem);
    free(self);
}

matrix* copy_matrix(matrix* self, matrix* src){
    assert(self->nrows >= src->nrows);
    assert(self->ncols == src->ncols);
    
    for (uint32_t i = 0; i < src->nrows; i++)
    {
        memcpy(self->elem[i], src->elem[i], src->ncols);
    }
    
    return self;
}

void export_matrix(matrix* self, uint8_t* dest){
    assert(self->ncols % 8 == 0);

    uint32_t byte_index = 0;
    for (uint32_t i = 0; i < self->nrows; i++){
        for (uint32_t j = 0; j < self->ncols; j+=8){
            dest[byte_index] = (self->elem[i][j+7] << 7) |
                                (self->elem[i][j+6] << 6) |
                                (self->elem[i][j+5] << 5) |
                                (self->elem[i][j+4] << 4) |
                                (self->elem[i][j+3] << 3) |
                                (self->elem[i][j+2] << 2) |
                                (self->elem[i][j+1] << 1) |
                                self->elem[i][j];
            byte_index++;
        }
    }
}

void import_matrix(matrix* self, const uint8_t* src){
    assert(self->ncols % 8 == 0);
    
    uint32_t byte_index = 0;
    for (uint32_t i = 0; i < self->nrows; i++){
        for (uint32_t j = 0; j < self->ncols; j+=8){
            self->elem[i][j+7] = (src[byte_index] >> 7) & 1;
            self->elem[i][j+6] = (src[byte_index] >> 6) & 1;
            self->elem[i][j+5] = (src[byte_index] >> 5) & 1;
            self->elem[i][j+4] = (src[byte_index] >> 4) & 1;
            self->elem[i][j+3] = (src[byte_index] >> 3) & 1;
            self->elem[i][j+2] = (src[byte_index] >> 2) & 1;
            self->elem[i][j+1] = (src[byte_index] >> 1) & 1;
            self->elem[i][j] = src[byte_index] & 1;
            byte_index++;
        }
    }
}

void row_addition_internal(matrix* self, const uint32_t r1, const uint32_t r2){
    for (uint32_t j = 0; j < self->ncols; j++)
    {
        uint8_t bit = get_element(self, r1, j) ^ get_element(self, r2, j);
        set_element(self, r1, j, bit);
    }
}

// make a matrix into rref form, inplace
matrix* rref(matrix* self, matrix* syndrome)
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
            if(syndrome != NULL){
                row_interchange(syndrome, succ_row_idx, row_idx);
            }
        }
                
        // By adding <succ_row_idx> th row in the other nrows 
        // s.t. self(i, <succ_row_idx>) == 1,
        // making previous columns as element row.
        for(uint32_t i = 0; i < self->nrows; ++i){
            if(i == succ_row_idx) continue;

            if(get_element(self, i, col_idx) == 1){
                row_addition_internal(self, i, succ_row_idx);
                if(syndrome != NULL){
                    row_interchange(syndrome, i, succ_row_idx);
                }
            }
        }
        row_idx = ++succ_row_idx;
    }
    //Gaussian elimination is finished. So return self.
    return self;
}

matrix* transpose(matrix *self, matrix* dest){
    for(uint32_t row = 0; row < dest->nrows; ++row){
        for(uint32_t col = 0; col < dest->ncols; ++col){
            set_element(dest, row, col, get_element(self, col, row));
        }
    }
        
    return dest;
}

int inverse(matrix *self, matrix *dest){
    if(self->nrows != self->ncols)  return INV_FAIL;
    if(dest->nrows != dest->ncols)  return INV_FAIL;
    if(self->nrows != dest->nrows)  return INV_FAIL;

    matrix* temp = new_matrix(self->nrows, self->ncols);
    copy_matrix(temp, self);

    uint32_t r, c;
    init_zero(dest);

    for (r = 0; r <  dest->nrows; ++r)
    {
        set_element(dest, r, r, 1);
    }

    for (c = 0; c < temp->ncols; ++c)
    {
        if(get_element(temp, c, c) == 0)
        {    
            for (r = c+1; r < self->nrows; ++r)
            {
                if(get_element(temp, r, c) != 0){
                    row_interchange(temp, r, c);
                    row_interchange(dest, r, c);
                    break;
                }
            }
            if(r >= temp->nrows)         return INV_FAIL;
        }
        
        for(r = 0; r < temp->nrows; r++){
            if(r == c) continue;
            if(get_element(temp, r, c) != 0){
                row_addition_internal(temp, r, c);
                row_addition_internal(dest, r, c);
            }
        }
    }
    
    // fprintf(stderr, "delete temp\n");
    delete_matrix(temp);
    return INV_SUCCESS;
}

int is_nonsingular(matrix *self){

    matrix* temp = new_matrix(self->nrows, self->ncols);
    copy_matrix(temp, self);

    uint32_t r, c;

    for (c = 0; c < temp->ncols; ++c)
    {
        if(get_element(temp, c, c) == 0)
        {    
            for (r = c+1; r < self->nrows; ++r)
            {
                if(get_element(temp, r, c) != 0){
                    row_interchange(temp, r, c);
                    break;
                }
            }
            if(r >= temp->nrows)         
                return INV_FAIL;
        }

        for(r = 0; r < temp->nrows; r++){
            if(r == c) continue;
            if(get_element(temp, r, c) != 0){
                row_addition_internal(temp, r, c);
            }
        }
    }

    delete_matrix(temp);
    return INV_SUCCESS;
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

void mat_mat_prod(matrix* self, matrix* mtx1, matrix* mtx2) {
    assert(mtx1->ncols == mtx2->nrows);
    
    for (uint32_t i = 0; i < mtx1->nrows; i++){ 
        for (uint32_t j = 0; j < mtx2->ncols; j++) {
            uint8_t val = 0;
            for (uint32_t k = 0; k < mtx1->ncols; k++)
                val ^= get_element(mtx1, i, k) & get_element(mtx2, k, j);
            set_element(self, i, j, val);
        }
    }
}

// assume vector is transposed
// self is also transposed
void vec_mat_prod(matrix* self, matrix* mat, matrix* vec){
    assert(mat->ncols == vec->ncols);
    assert(self->ncols == mat->nrows);

    for(uint32_t i = 0; i < mat->nrows; i++){
        uint8_t bit = 0;
        for (uint32_t j = 0; j < mat->ncols; j++)
        {
            bit ^= get_element(mat, i, j) & get_element(vec, 0, j);
        }
        
        set_element(self, 0, i, bit);
    }
}

void mat_mat_add(matrix* self, matrix *mat1, matrix *mat2){
    assert((mat1->nrows == mat2->nrows) && (mat1->ncols == mat2->ncols));
    
    for (uint32_t i = 0; i < self->nrows; i++)
    {
        for (uint32_t j = 0; j < self->ncols; j++)
        {
            uint8_t bit = get_element(mat1, i, j) ^ get_element(mat2, i, j);
            set_element(self, i, j, bit);
        }
        
        
    }
}

void dual(matrix* self, matrix* dual_sys){

    uint16_t lead[self->nrows];
    uint16_t lead_diff[self->ncols - self->nrows];    

    init_zero(dual_sys);

    rref(self, NULL);
    get_pivot(self, lead, lead_diff);

    // Fill not-identity part (P')
    for (uint32_t row = 0; row < dual_sys->nrows; row++) 
        for (uint32_t col = 0; col < self->nrows; col++) 
            set_element(dual_sys, row, lead[col], get_element(self, col, lead_diff[row]));
    
    for (uint32_t row = 0; row < dual_sys->nrows; row++) 
            set_element(dual_sys, row, lead_diff[row], 1);    
}

void row_interchange(matrix* self, uint32_t row1, uint32_t row2) {
    uint8_t* temp = self->elem[row1];
    self->elem[row1] = self->elem[row2];
    self->elem[row2] = temp;
}

void partial_replace(matrix* self, const uint32_t r1, const uint32_t r2,
        const uint32_t c1, const uint32_t c2, 
        matrix* src, const int r3, const int c3){
    for(uint32_t i = 0; i < r2 - r1; i++)
        for(uint32_t j = 0; j < c2 - c1; j++)
            set_element(self, r1 + i, c1+j, get_element(src, r3 + i, c3 + j));
}

void codeword(matrix* self, uint8_t* seed, matrix* dest){
    for (uint32_t i = 0; i < self->nrows; i++)
    {
        uint32_t byte_index = i >> 3; // byte_index = bit_index / 8;
        uint32_t bit_offset = i & 7; // bit_offset = bit_index % 8;
        uint8_t rand_bit = (seed[byte_index] >> bit_offset) & 1;
        if (rand_bit == 1)
        {
            for (uint32_t j = 0; j < self->ncols; j++)
            {
                dest->elem[0][j] ^= self->elem[i][j];
            }
        }
    }
}

uint8_t is_zero(matrix* self){
    for (uint32_t i = 0; i < self->nrows; i++)
    {
        for (size_t j = 0; j < self->ncols; j++)
        {
            if(get_element(self, i, j) != 0){
                return 0;
            }
        }
    }
    return 1;
}

// Function to perform Gaussian elimination
uint16_t rank(const matrix* self) {
    matrix* copy = new_matrix(self->nrows, self->ncols);
    for (uint16_t i = 0; i < self->nrows; ++i) {
        for (uint16_t j = 0; j < self->ncols; ++j) {
            set_element(copy, i, j, get_element(self, i, j));
        }
    }

    uint16_t rank = 0;  // Initialize rank as 0

    for (uint16_t r = 0; r < copy->nrows; ++r) {
        uint16_t lead = 0;  // Current leading column

        while (lead < copy->ncols) {
            uint16_t i = r;

            while (i < copy->nrows && get_element(copy, i, lead) == 0) {
                ++i;
            }

            if (i < copy->nrows) {
                row_interchange(copy, r, i);

                for (uint16_t j = r + 1; j < copy->nrows; ++j) {
                    if (get_element(copy, j, lead) != 0) {
                        for (uint16_t k = lead; k < copy->ncols; ++k) {
                            set_element(copy, j, k, get_element(copy, j, k) ^ get_element(copy, r, k));
                        }
                    }
                }

                ++rank;
                break;
            }

            ++lead;
        }
    }

    delete_matrix(copy);

    return rank;
}

uint32_t size_in_byte(const matrix* self){
    // Assume ncols is a multiple of 8
    // otherwise, do a ceiling
    return self->nrows * (self->ncols + 7)/8;
}