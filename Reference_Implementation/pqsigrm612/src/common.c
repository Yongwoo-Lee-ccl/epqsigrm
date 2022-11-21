// common.c for pqsigRM
#include "common.h"

unsigned char* hash_message(unsigned char *s, const unsigned char *m, 
	unsigned long long mlen, unsigned long long sign_i){
	// Hash the given message
	// syndrome s = h(h(M)|i) | (h(h(M)|i)) | ...
	unsigned char* stemp = (unsigned char*)malloc(HASHSIZEBYTES*4);

	SHA512(m, mlen, stemp);
	*(unsigned long long*)(stemp+HASHSIZEBYTES) = sign_i;// concatenate i i.e. h(M)|i
	
	SHA512(stemp, HASHSIZEBYTES+sizeof(unsigned long long), stemp); //h(h(M)|i)
	SHA512(stemp                , HASHSIZEBYTES, stemp+HASHSIZEBYTES);//(h(h(M)|i))
	SHA512(stemp+HASHSIZEBYTES  , HASHSIZEBYTES, stemp+HASHSIZEBYTES*2);
	SHA512(stemp+HASHSIZEBYTES*2, HASHSIZEBYTES, stemp+HASHSIZEBYTES*3);
	memcpy(s, stemp, 1+(CODE_N-CODE_K-1)/8);
	
	free(stemp);
	return s;
}

int hamming_weight(matrix* error){
	assert(error->ncols % ELEMBLOCKSIZE == 0);

    uint64_t wgt = 0;

    for (int i = 0; i < error->words_in_row; ++i){
        wgt += (1 & (error->elem[i]))
            + (1 & (error->elem[i] >> 1))
            + (1 & (error->elem[i] >> 2))
            + (1 & (error->elem[i] >> 3))
            + (1 & (error->elem[i] >> 4))
            + (1 & (error->elem[i] >> 5))
            + (1 & (error->elem[i] >> 6))
            + (1 & (error->elem[i] >> 7))
            + (1 & (error->elem[i] >> 8))
            + (1 & (error->elem[i] >> 9))
            + (1 & (error->elem[i] >> 10))
            + (1 & (error->elem[i] >> 11))
            + (1 & (error->elem[i] >> 12))
            + (1 & (error->elem[i] >> 13))
            + (1 & (error->elem[i] >> 14))
            + (1 & (error->elem[i] >> 15))
            + (1 & (error->elem[i] >> 16))
            + (1 & (error->elem[i] >> 17))
            + (1 & (error->elem[i] >> 18))
            + (1 & (error->elem[i] >> 19))
            + (1 & (error->elem[i] >> 20))
            + (1 & (error->elem[i] >> 21))
            + (1 & (error->elem[i] >> 22))
            + (1 & (error->elem[i] >> 23))
            + (1 & (error->elem[i] >> 24))
            + (1 & (error->elem[i] >> 25))
            + (1 & (error->elem[i] >> 26))
            + (1 & (error->elem[i] >> 27))
            + (1 & (error->elem[i] >> 28))
            + (1 & (error->elem[i] >> 29))
            + (1 & (error->elem[i] >> 30))
            + (1 & (error->elem[i] >> 31))
            + (1 & (error->elem[i] >> 32))
            + (1 & (error->elem[i] >> 33))
            + (1 & (error->elem[i] >> 34))
            + (1 & (error->elem[i] >> 35))
            + (1 & (error->elem[i] >> 36))
            + (1 & (error->elem[i] >> 37))
            + (1 & (error->elem[i] >> 38))
            + (1 & (error->elem[i] >> 39))
            + (1 & (error->elem[i] >> 40))
            + (1 & (error->elem[i] >> 41))
            + (1 & (error->elem[i] >> 42))
            + (1 & (error->elem[i] >> 43))
            + (1 & (error->elem[i] >> 44))
            + (1 & (error->elem[i] >> 45))
            + (1 & (error->elem[i] >> 46))
            + (1 & (error->elem[i] >> 47))
            + (1 & (error->elem[i] >> 48))
            + (1 & (error->elem[i] >> 49))
            + (1 & (error->elem[i] >> 50))
            + (1 & (error->elem[i] >> 51))
            + (1 & (error->elem[i] >> 52))
            + (1 & (error->elem[i] >> 53))
            + (1 & (error->elem[i] >> 54))
            + (1 & (error->elem[i] >> 55))
            + (1 & (error->elem[i] >> 56))
            + (1 & (error->elem[i] >> 57))
            + (1 & (error->elem[i] >> 58))
            + (1 & (error->elem[i] >> 59))
            + (1 & (error->elem[i] >> 60))
            + (1 & (error->elem[i] >> 61))
            + (1 & (error->elem[i] >> 62))
            + (1 & (error->elem[i] >> 63));
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

void col_permute(matrix* m, const int rf, const int rr, const int cf, 
	const int cr, uint16_t* Q)
{
	matrix* mcpy = new_matrix(m->nrows, m->ncols); 
	memcpy(mcpy->elem, m->elem, m->alloc_size);

	for(uint32_t c = cf; c < cr; c++){
		for(uint32_t r = rf; r < rr; r++){
			set_element(m, r, c, get_element(mcpy, r, cf + Q[c-cf]));
		}
	}

	delete_matrix(mcpy);
}



