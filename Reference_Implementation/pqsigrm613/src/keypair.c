#include "api.h"
#include "common.h"

void export_sk(unsigned char *sk,uint16_t *Q, uint16_t *part_perm1, uint16_t* part_perm2, matrix* Hrep){
	//export private in order: Q, part_perm1, pert_perm2
	memcpy		(sk, 
					Q, sizeof(uint16_t)*CODE_N);
	memcpy		(sk+sizeof(uint16_t)*CODE_N, 
					part_perm1, sizeof(uint16_t)*CODE_N/4);
	memcpy		(sk+sizeof(uint16_t)*CODE_N+sizeof(uint16_t)*CODE_N/4, 
					part_perm2, sizeof(uint16_t)*CODE_N/4);
	export_matrix(Hrep, sk + sizeof(uint16_t)*CODE_N + (sizeof(uint16_t)*CODE_N/4)*2);
}

void export_pk(unsigned char *pk, matrix *Hpub){
	export_matrix(Hpub, pk);
}

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk){
	// nrows of Gm is set to K + 1 for public key
	matrix* Gm = new_matrix(CODE_K, CODE_N);
	matrix* Hm = new_matrix(CODE_N - CODE_K, CODE_N);

	matrix* Gpub = new_matrix(CODE_K + 1, CODE_N);
	matrix* Hpub = new_matrix(CODE_N - CODE_K - 1, CODE_N);

	uint16_t Q[CODE_N];
	
	uint16_t part_perm1[(CODE_N/4)];
	uint16_t part_perm2[(CODE_N/4)];

	uint16_t s_lead[CODE_N - CODE_K];
	uint16_t s_diff[CODE_K];

	// generate secret parital permutations
	partial_permutation_gen(part_perm1);
	partial_permutation_gen(part_perm2);
	
	// Generate a partially permute generator matrix Gm
	rm_gen(Gm, RM_R, RM_M, 0, CODE_K, 0, CODE_N);

	// replace of RM(r,r)
	matrix* Grep = new_matrix((1<<RM_R) - K_REP, (1<<RM_R));
	matrix* Hrep = new_matrix(K_REP, 1<<RM_R);

	uint8_t is_odd = 0;
	uint8_t randstr[Grep->ncols * Grep->nrows/8];

	// check if Hrep has a odd row.
	while(1){
		randombytes(randstr, Grep->ncols * Grep->nrows/8);
		randomize(Grep, randstr);
		
		dual(Grep, Hrep);
		for (uint32_t i = 0; i < Hrep->nrows; i++)
		{	
			uint8_t parity = 0;
			for (uint32_t j = 0; j < Hrep->ncols; j++)
			{
				parity ^= get_element(Hrep, i, j);
			}
			if (parity == 1){
				is_odd = 1;
				break;
			}
		}
		if(is_odd) break;
	}
	rref(Hrep);
	dual(Hrep, Grep);

	// TODO: two random rows 

	// replace the code (starting from second row)
	for (uint32_t i = 0; i < CODE_N; i += Grep->ncols)
	{
		partial_replace(Gm, K_REP, K_REP + Grep->nrows, i, i + Grep->ncols, Grep, 0, 0); 
	}
	
	// Partial permutation
	for (uint32_t i = 0; i < 4; ++i)
	{
		col_permute(Gm, 0, rm_dim[RM_R][RM_M -2], 
			i*(CODE_N/4),(i+1)*(CODE_N/4), part_perm1);
	}
	col_permute(Gm, CODE_K - rm_dim[RM_R-2][RM_M-2], CODE_K, 
		3*CODE_N/4, CODE_N, part_perm2);
	
	// Parity check matrix of the modified RM code
	dual(Gm, Hm);

	// pick a random codeword from the dual code
	matrix* code_from_dual = new_matrix(1, Hm->ncols);
	matrix* random_syndrome = new_matrix(1, Hm->nrows);
	uint8_t seed[(Hm->nrows + 7)/8];
	while(1){
		randombytes(seed, (Hm->nrows + 7)/8);
		codeword(Hm, seed, code_from_dual);
		
		vec_mat_prod(random_syndrome, Hm, code_from_dual);
		if(! is_zero(random_syndrome)){
			break;
		}
	}
	delete_matrix(random_syndrome);

	copy_matrix(Gpub, Gm);
	// partial_replace(Gpub, CODE_K, CODE_K + 1, 0, CODE_N, code_from_dual, 0, 0);
	permutation_gen(Q, CODE_N);

	// Generate the dual code of Gm and the public key
	dual(Gpub, Hpub);
	matrix* Hcpy = new_matrix(Hpub->nrows, Hpub->ncols); 
	copy_matrix(Hcpy, Hpub);

	col_permute(Hcpy, 0, Hcpy->nrows, 0, Hcpy->ncols, Q);
	rref(Hcpy);

	uint16_t pivot[Hcpy->nrows];
	uint16_t d_pivot[Hcpy->ncols - Hcpy->nrows];
	get_pivot(Hcpy, pivot, d_pivot);

	for (uint32_t i = 0; i < Hcpy->nrows; i++)
	{
		if(pivot[i] != i){
			uint16_t tmp = Q[i];
			Q[i] = Q[pivot[i]];
			Q[pivot[i]] = tmp;
		}		
	}
	
	col_permute(Hpub, 0, Hpub->nrows, 0, Hpub->ncols, Q);

	rref(Hpub);
	// matrix* Sinv = new_matrix(Hpub->nrows, Hpub->nrows);
	// inverse(S, Sinv);
	//printf("rank: %d / %d\n", rank(Hpub), Hpub->nrows);

	// export_sk(sk, Q, part_perm1, part_perm2, Hrep, Sinv);
	export_sk(sk, Q, part_perm1, part_perm2, Hrep);
	export_pk(pk, Hpub);

	//print_matrix_keypair(Hpub, 1000, 1064, 2000, 2064);

	delete_matrix(Gm);
	delete_matrix(Hm);
	delete_matrix(Gpub);
	delete_matrix(Hpub);
	delete_matrix(Grep);
	delete_matrix(Hrep);
	delete_matrix(code_from_dual);
	delete_matrix(Hcpy);

	return 0;
}
