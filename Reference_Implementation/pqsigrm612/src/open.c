#include "api.h"
#include "common.h"

char* convertToHexString2(const unsigned char* array, size_t length) {
    char* hexString = (char*) malloc(length * 2 + 1);  // Allocate memory for the hex string

    if (hexString == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    for (size_t i = 0; i < length; ++i) {
        sprintf(hexString + (i * 2), "%02X", array[i]);  // Convert each byte to a 2-digit hexadecimal number
    }

    hexString[length * 2] = '\0';  // Null-terminate the string

    return hexString;
}


void print_partial_open(matrix* mat, uint32_t r1, uint32_t r2, uint32_t c1, uint32_t c2){
	for(uint32_t i = r1; i < r2; i++){
		for(uint32_t j = c1; j < c2; j++){
			printf("%u", get_element(mat, i, j));
		}
		printf("\n");
	}
}

void import_pk(const unsigned char *pk, matrix *H_pub){
    import_matrix(H_pub, pk);
}

void print_matrix_open(matrix* mtx){
    uint32_t row = mtx->nrows;
    uint32_t col = mtx->ncols;
    for (size_t i = 0; i < row; i++)
    {
        for (size_t j = 0; j < col; j++)
        {
            printf("%d",get_element(mtx, i,j));
        }printf("\n");
    }
    
}

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk){
    matrix *sign = new_matrix(1, CODE_N);
    matrix *Hpub = new_matrix(CODE_N - CODE_K - 1, CODE_N);

    matrix *syndrome_by_hash = new_matrix(1, CODE_N - CODE_K - 1);
    matrix *syndrome_by_e     = new_matrix(1, CODE_N - CODE_K - 1);

    uint64_t sign_i;
    uint64_t mlen_rx;
    unsigned char* m_rx;
    
    // import signed msg 
    // sign is (mlen, M, e, sign_i)
    mlen_rx = *(uint64_t*)sm;
    m_rx = (unsigned char*)malloc(mlen_rx);
    memcpy(m_rx, sm + sizeof(uint64_t), mlen_rx);
    import_matrix(sign, sm + sizeof(uint64_t) + mlen_rx);

    sign_i = *(uint64_t*)(sm + sizeof(uint64_t) + mlen_rx + sign->nrows * sign->ncols/8);    


    
    if(hamming_weight(sign) > WEIGHT_PUB) {
        fprintf(stderr, "larger weight\n");
        // return VERIF_REJECT;
    }
    
    //import public key
    import_pk(pk, Hpub);

    uint8_t randstr[syndrome_by_hash->ncols/8 + 1];
    hash_message(randstr, m_rx, mlen_rx, sign_i);
    randomize(syndrome_by_hash, randstr);


    
    matrix* test_e = new_matrix(1, Hpub->ncols);
    matrix* test_s = new_matrix(1, Hpub->nrows);

    for (uint32_t i = 0; i < syndrome_by_hash->ncols; i++)
    {
        uint8_t bit = get_element(syndrome_by_hash, 0, i);
        set_element(test_e, 0, i, bit);
    }

    vec_mat_prod(test_s, Hpub, test_e);


    mat_mat_add(test_s, syndrome_by_hash, test_s);
    printf("is error generated correctly: %d\n", is_zero(test_s));

    vec_mat_prod(syndrome_by_e, Hpub, sign);

    printf("is a codeword: %d\n", is_zero(syndrome_by_e));
    // int count_diff = 0;
    for(uint32_t i=0; i < CODE_N-CODE_K - 1; ++i){
        if(get_element(syndrome_by_hash, 0, i) != get_element(syndrome_by_e, 0, i)){
            fprintf(stderr, "different hash at %d\n", i);
            // count_diff += 1;
            return VERIF_REJECT;
        }
    }
    memcpy(m, m_rx, mlen_rx);
    *mlen = mlen_rx;

    delete_matrix(sign);
    delete_matrix(Hpub);

    delete_matrix(syndrome_by_hash);
    delete_matrix(syndrome_by_e);
    free(m_rx);

    return 0;
}
