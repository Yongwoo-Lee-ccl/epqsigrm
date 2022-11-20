#include "api.h"
#include "common.h"

void export_sk(unsigned char *sk,uint16_t *Q, uint16_t *part_perm1, uint16_t* part_perm2){
	//export private in order: Q, part_perm1, pert_perm2
	memcpy		(sk, 
					Q, sizeof(uint16_t)*CODE_N);
	memcpy		(sk+sizeof(uint16_t)*CODE_N, 
					part_perm1, sizeof(uint16_t)*CODE_N/4);
	memcpy		(sk+sizeof(uint16_t)*CODE_N+sizeof(uint16_t)*CODE_N/4, 
					part_perm2, sizeof(uint16_t)*CODE_N/4);
}

void export_pk(unsigned char *pk, matrix *H_pub){
	export_matrix(pk, H_pub);
}

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk){
	matrix* G_M = new_matrix(CODE_K, CODE_N);

	matrix* H_M = new_matrix(CODE_N - CODE_K, CODE_N);
	matrix* H_pub = new_matrix(CODE_N - CODE_K, CODE_N);

	uint16_t Q[CODE_N];
	
	uint16_t part_perm1[(CODE_N/4)];
	uint16_t part_perm2[(CODE_N/4)];

	uint16_t s_lead[CODE_N - CODE_K];
	uint16_t s_diff[CODE_K];

	// generate secret parital permutations
	partial_permutation_gen(part_perm1);
	partial_permutation_gen(part_perm2);
	
	// Generate a partially permute generator matrix G_M
	// fprintf(stderr, "gen mod start\n");
	rm_gen_mod(G_M, part_perm1, part_perm2);
	// fprintf(stderr, "gen mod\n");

	// Parity check matrix of the modified RM code
	dual(G_M, H_M, 0, 0);
	

	// Generate a Scrambling matrix and its inverse. 
	permutation_gen(Q, CODE_N);
	
	matrix* Hcpy = new_matrix(H_M->nrows, H_M->ncols); 
	memcpy(Hcpy->elem, H_M->elem, H_M->alloc_size);
	rref(Hcpy); 
	get_pivot(Hcpy, s_lead, s_diff);

	col_permute(Hcpy, 0, CODE_N-CODE_K, 0, CODE_N, Q);
	rref(Hcpy);

	uint16_t pivot[CODE_N - CODE_K];
	uint16_t d_pivot[CODE_K];
	get_pivot(Hcpy, pivot, d_pivot);

	for (uint32_t i = 0; i < CODE_N - CODE_K; i++)
	{
		if(pivot[i] != i){
			uint16_t tmp = Q[i];
			Q[i] = Q[pivot[i]];
			Q[pivot[i]] = tmp;
		}		
	}
	
	col_permute(H_M, 0, CODE_N-CODE_K, 0, CODE_N, Q);
	rref(H_M);

	for (uint32_t i = 0; i < CODE_N - CODE_K; i++)
	{
		if(get_element(H_M, i, i) != 1){
			printf("not identity!, %d, %d, %d\n", i, i, get_element(H_M, i, i));
		}		
	}

	copy_matrix(H_pub, H_M);
	// fprintf(stderr, "Hpubgen\n");
	
	// fprintf(stderr, "mat mul\n");	

	export_sk(sk, Q, part_perm1, part_perm2);

	// printf("slead_cpy:\n");
	// for (size_t i = 0; i < CODE_N - CODE_K; i++)
	// {
	// 	printf("%4d ", slead_cpy[i]);
	// }printf("\n");
	// printf("sk: %p, pk: %p\n", sk, pk);
	
	export_pk(pk, H_pub);

	delete_matrix(G_M);
	delete_matrix(H_M);
	delete_matrix(H_pub); 

	return 0;
}
