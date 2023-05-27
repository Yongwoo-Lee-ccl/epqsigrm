// common.c for pqsigRM
#include "common.h"

void SHAKE256(unsigned char *digest, const unsigned char *message, size_t message_len, size_t hashlen)
{
    EVP_MD_CTX *mdctx;
    mdctx = EVP_MD_CTX_create();
    EVP_DigestInit_ex(mdctx, EVP_shake256(), NULL);
    EVP_DigestUpdate(mdctx, message, message_len);
    EVP_DigestFinalXOF(mdctx, digest, hashlen);
    EVP_MD_CTX_destroy(mdctx);
}

unsigned char* hash_message(unsigned char *s, const unsigned char *m, 
	unsigned long long mlen, unsigned long long sign_i){
	// Hash the given message
	// syndrome s = h(h(M)|i)
	unsigned char* buffer = (unsigned char*)malloc(mlen+sizeof(uint64_t));
	memcpy(buffer, m, mlen);
	memcpy(buffer + mlen, (unsigned char*)(&sign_i), sizeof(unsigned long long));
	SHAKE256(s, buffer, mlen+sizeof(uint64_t), (1 + (CODE_N - CODE_K - 1)/64)*64/8);

	free(buffer);
	return s;
}

int hamming_weight(matrix* error){

    uint64_t wgt = 0;
    unsigned char mask = 1;
    for (int c = 0; c < error->ncols; ++c){
        wgt += error->elem[0][c] & mask;
    }

    return wgt;
}

void swap16(uint16_t *Q, const int i, const int j){
	uint16_t temp;
	temp = Q[i];
	Q[i] = Q[j];
	Q[j] = temp;
}

void permutation_gen(uint16_t *Q, uint32_t len){
	uint32_t i; 
	for(i=0; i<len; i++){
		Q[i] = i;
	}
	for(i=0; i<len; i++){
		swap16(Q, i, random16(len));
	}
}

int static compare(const void* first, const void* second){
	return (*(uint16_t*)first > *(uint16_t*)second)?1:-1;
}

void partial_permutation_gen(uint16_t* Q){
	permutation_gen(Q, CODE_N/4);
	uint16_t partial_elem[PARM_P];
	uint16_t partial_perm[PARM_P];

	memcpy(partial_perm, Q, sizeof(uint16_t)*PARM_P);
	memcpy(partial_elem, Q, sizeof(uint16_t)*PARM_P);

	qsort(partial_elem, PARM_P, sizeof(uint16_t), compare);
	qsort(Q, CODE_N/4, sizeof(uint16_t), compare);

	for (uint32_t i = 0; i < PARM_P; ++i)
	{
		Q[partial_elem[i]] = partial_perm[i];
	}
}

uint16_t random16(uint16_t n){
	uint16_t r;
	randombytes((unsigned char*)&r, 2);
	return r%n;
}

void col_permute(matrix* m, const int r1, const int r2
	, const int c1, const int c2, uint16_t* Q)
{	
	matrix* copy = new_matrix(r2 - r1, c2 - c1);
	for (uint32_t r = 0; r < r2 - r1; r++)
	{
		for (uint32_t c = 0; c < c2 - c1; c++)
		{
			uint8_t bit = get_element(m, r1 + r, c1 + c);
			set_element(copy, r, c, bit);
		}
	}
	
	for(uint32_t c = 0; c < c2 - c1; c++){
		for(uint32_t r = 0; r < r2 - r1; r++){
            uint8_t bit =  get_element(copy, r, Q[c]);
			set_element(m, r1 + r, c1 + c, bit);
		}
	}

	delete_matrix(copy);
}



