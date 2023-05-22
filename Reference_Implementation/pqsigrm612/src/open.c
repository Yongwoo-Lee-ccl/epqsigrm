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

	matrix *H_pub = new_matrix(CODE_N - CODE_K - 1, CODE_N);

	matrix *syndrome_by_hash = new_matrix(1, CODE_N - CODE_K - 1);
	matrix *syndrome_by_e	 = new_matrix(1, CODE_N - CODE_K - 1);

	uint64_t sign_i;
	uint64_t mlen_rx;
	unsigned char* m_rx;
	
	// import signed msg 
	// sign is (mlen, M, e, sign_i)
	mlen_rx = *(unsigned long long*)sm;
	m_rx = (unsigned char*)malloc(mlen_rx);
	memcpy(m_rx, sm + sizeof(uint64_t), mlen_rx);
	import_matrix(sign, sm + sizeof(unsigned long long) + mlen_rx);
	sign_i = *(unsigned long long*)(sm + sizeof(unsigned long long) + mlen_rx + sign->nrows * sign->ncols/8);	
	
	if(hamming_weight(sign) > WEIGHT_PUB) {
		fprintf(stderr, "larger weight\n");
		return VERIF_REJECT;
	}
	
	uint8_t randstr[syndrome_by_hash->ncols/8 + 1];
	printf("open!!!!\n%s\nmlen: %llu\nsign_i: %lu\n", convertToHexString2(m_rx, mlen_rx), mlen_rx, sign_i);
	hash_message(randstr, m_rx, mlen_rx, sign_i);
	randomize(syndrome_by_hash, randstr);
	
	//import public key
	import_pk(pk, H_pub);

	vec_mat_prod(syndrome_by_e, H_pub, sign);

	for(uint32_t i=0; i < CODE_N-CODE_K - 1; ++i){
		if(get_element(syndrome_by_hash, 0, i) != get_element(syndrome_by_e, 0, i)){
			fprintf(stderr, "different hash at %d\n", i);
			// printf("hashed value: \n");
			// print_matrix_open(syndrome_by_hash);
			// printf("H*e value: \n");
			// print_matrix_open(syndrome_by_e);
			// return VERIF_REJECT;
		}
		return VERIF_REJECT;
	}
	memcpy(m, m_rx, mlen_rx);
	*mlen = mlen_rx;

	delete_matrix(sign);
	delete_matrix(H_pub);

	delete_matrix(syndrome_by_hash);
	delete_matrix(syndrome_by_e);
	free(m_rx);

	return 0;
}
