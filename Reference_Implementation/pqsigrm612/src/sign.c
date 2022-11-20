#include "api.h"
#include "common.h"
#include "nearest_vector.h"

int wgt(float *yc, float *yr)
{
	int i, w=0;
	for(i=0; i<CODE_N; i++)
		if(yc[i] != yr[i]) w++;
	return w;
}

void print_matrix_sign(matrix* mtx){
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

void import_sk(const unsigned char *sk, matrix *Sinv
		, uint16_t **Q, uint16_t **part_perm1, uint16_t **part_perm2
		, uint16_t **s_lead)
{
	import_matrix(Sinv, sk);
	*Q 			= (uint16_t*)(sk+Sinv->alloc_size);
	*part_perm1 = (uint16_t*)(sk+Sinv->alloc_size+sizeof(uint16_t)*CODE_N);
	*part_perm2 = (uint16_t*)(sk+Sinv->alloc_size+sizeof(uint16_t)*CODE_N
					+sizeof(uint16_t)*CODE_N/4);
	*s_lead 	= (uint16_t*)(sk+Sinv->alloc_size+sizeof(uint16_t)*CODE_N
					+sizeof(uint16_t)*CODE_N/2);
}

/*
* yr and yc are equal at the end.
*/
void y_init(float *yc, float *yr, matrix* syndrome, uint16_t *Q){
		for(uint32_t i=0; i < CODE_N-CODE_K; i++) {
			yc[i] = (get_element(syndrome, 0, i) == 0)? 1.:-1.;
		}
		for(uint32_t i=CODE_N-CODE_K; i < CODE_N; i++) {
			yc[i] = 1.;
		}

		// yr first, yc next
		for(uint32_t i =0; i < CODE_N; i++) {
			yr[Q[i]] = yc[i];
		}
		for(uint32_t i =0; i < CODE_N; i++) {
			yc[i] = yr[i];
		}
}

int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
	const unsigned char *m, unsigned long long mlen,
	const unsigned char *sk){

	// read secret key(bit stream) into appropriate type.
	matrix* Sinv = new_matrix(CODE_N - CODE_K, CODE_N - CODE_K);
	uint16_t *Q, *part_perm1, *part_perm2, *s_lead;

	import_sk(sk, Sinv, &Q, &part_perm1, &part_perm2, &s_lead);
	
	// assert(0);
	// Do signing, decode until the a error vector wt <= w is achieved
	
	uint64_t sign_i;

	// unsigned char sign[CODE_N];
	matrix *synd_mtx= new_matrix(1, CODE_N - CODE_K);

	float yc[CODE_N], yr[CODE_N];
	
	init_decoding(CODE_N);
	while(1){
		//random number
		randombytes((unsigned char*)&sign_i, sizeof(uint64_t));
		// Find syndrome
		hash_message((unsigned char*)synd_mtx->elem, m, mlen, sign_i);
		y_init(yc, yr, synd_mtx, Q);
		// printf("s_lead:\n");
		// for (size_t i = 0; i < CODE_N - CODE_K; i++)
		// {
		// 	printf("%4d ", s_lead[i]);
		// }
		// printf("\nyc\n");
		// for (size_t i = 0; i < CODE_N; i++)
		// {
		// 	printf("%d", (yc[i] == -1)?1:0);
		// }printf("\nsynd\n");
		// for (size_t i = 0; i < CODE_N-CODE_K; i++)
		// {
		// 	printf("%d", get_element(scrambled_synd_mtx, 0, i));
		// }
		
		// decode and find e
		// In the recursive decoding procedure,
		// Y is 1 when the received codeword is 0, o.w, -1
		recursive_decoding_mod(yc, RM_R, RM_M, 0, CODE_N, part_perm1, part_perm2);
		
		// Check Hamming weight of e'
		if(wgt(yr, yc) <= WEIGHT_PUB) break;
	}

	// printf("syndrome\n");
	// print_matrix_sign(synd_mtx);
	// compute Qinv*e'
	matrix *sign = new_matrix(1, CODE_N);
	for(uint32_t i=0; i < CODE_N; i++){
		set_element(sign, 0, i, (yr[Q[i]] != yc[Q[i]]));
	}
	// printf("sign:\n");
	// print_matrix_sign(sign);
	// export message
	// sing is (mlen, M, e, sign_i)
	// M includes its length, i.e., mlen
	*(unsigned long long*)sm = mlen;
	memcpy(sm+sizeof(unsigned long long), m, mlen);

	memcpy(sm+sizeof(unsigned long long)+mlen, sign->elem, sign->alloc_size);
	*(unsigned long long*)(sm + sizeof(unsigned long long) + mlen + sign->alloc_size) 
		= sign_i;

	*smlen = sizeof(unsigned long long) + mlen + sign->alloc_size + sizeof(unsigned long long);
	
	delete_matrix(Sinv);
	delete_matrix(synd_mtx);

	return 0;	
}
