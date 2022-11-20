#include "api.h"
#include "common.h"

void export_sk(unsigned char *sk, matrix *Sinv, uint16_t *Q, 
	uint16_t *part_perm1, uint16_t* part_perm2, uint16_t *s_lead){
	//export private in order: Sinv, Q, part_perm1, pert_perm2, pivot of H_M
	export_matrix(sk        			, Sinv);
	
	memcpy		(sk+Sinv->alloc_size, 
									Q, sizeof(uint16_t)*CODE_N);
	memcpy		(sk+Sinv->alloc_size+sizeof(uint16_t)*CODE_N, 
								   part_perm1, sizeof(uint16_t)*CODE_N/4);

	memcpy		(sk+Sinv->alloc_size+sizeof(uint16_t)*CODE_N+sizeof(uint16_t)*CODE_N/4, 
								   part_perm2, sizeof(uint16_t)*CODE_N/4);
	memcpy		(sk+Sinv->alloc_size+sizeof(uint16_t)*CODE_N+(sizeof(uint16_t)*CODE_N/4)*2, 
								   s_lead, sizeof(uint16_t)*(CODE_N - CODE_K));
}

void export_pk(unsigned char *pk, matrix *H_pub){
	export_matrix(pk, H_pub);
}

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk){
	matrix* G_M = new_matrix(CODE_K, CODE_N);

	matrix* H_M = new_matrix(CODE_N - CODE_K, CODE_N);
	matrix* H_pub = new_matrix(CODE_N - CODE_K, CODE_N);

	matrix* S = new_matrix(CODE_N - CODE_K, CODE_N - CODE_K);
	matrix* Sinv = new_matrix(CODE_N - CODE_K, CODE_N - CODE_K);

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
	rref(H_M); 
	get_pivot(H_M, s_lead, s_diff);
	// printf("s_lead:\n");
	// for (size_t i = 0; i < CODE_N - CODE_K; i++)
	// {
	// 	printf("%4d ", s_lead[i]);
	// }printf("\n");

	// fprintf(stderr, "pivot\n");

	// Generate a Scrambling matrix and its inverse. 
	do{
		randombytes((unsigned char*)S->elem, S->alloc_size);
	}while(is_nonsingular(S) != INV_SUCCESS);
	inverse(S, Sinv);
	// fprintf(stderr, "s inv\n");
	permutation_gen(Q, CODE_N);
	// fprintf(stderr, "Qgen\n");
	col_permute(H_M, 0, CODE_N-CODE_K, 0, CODE_N, Q);
	copy_matrix(H_pub, H_M);
	// fprintf(stderr, "Hpubgen\n");
	
	// fprintf(stderr, "mat mul\n");	

	export_sk(sk, Sinv, Q, part_perm1, part_perm2, s_lead);

	// printf("slead_cpy:\n");
	// for (size_t i = 0; i < CODE_N - CODE_K; i++)
	// {
	// 	printf("%4d ", slead_cpy[i]);
	// }printf("\n");
	// printf("sk: %p, pk: %p\n", sk, pk);
	
	mat_mat_prod(S, H_M, H_pub);
	export_pk(pk, H_pub);

	delete_matrix(G_M);
	delete_matrix(H_M);
	delete_matrix(H_pub); 
	delete_matrix(S);
	delete_matrix(Sinv);

	return 0;
}
