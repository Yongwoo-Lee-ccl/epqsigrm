#include "api.h"
#include "common.h"

void import_signed_msg(matrix *errorMtx, unsigned long long *sign_i, const unsigned char *sm){
	import_matrix(errorMtx, sm);
	*sign_i = *((unsigned long long*)(sm+ERRORSIZEBYTES));
}


void import_pk(const unsigned char *pk, matrix *H_pub){
	import_matrix(H_pub, pk);
}

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk){
	matrix *errorMtx = new_matrix(1, CODE_N);

	matrix *H_pub = new_matrix(CODE_N-CODE_K, CODE_N);

	matrix *syndrome_by_hash = new_matrix(1, CODE_N - CODE_K);
	matrix *syndrome_by_e	 = new_matrix(1, CODE_N - CODE_K);

	unsigned long long sign_i;
	unsigned long long mlen_rx;
	unsigned char* m_rx;
	
	int i;

	memcpy(&mlen_rx, sm, sizeof(unsigned long long));
	m_rx = (unsigned char*)malloc(mlen_rx);

	memcpy(m_rx, sm + sizeof(unsigned long long), mlen_rx);

	import_signed_msg(errorMtx, &sign_i, sm + sizeof(unsigned long long) + mlen_rx);
	
	if(hamming_weight(errorMtx) > WEIGHT_PUB) 
		return VERIF_REJECT;
	

	hash_message(syndrome_by_hash->elem, m_rx, mlen_rx, sign_i);
	
	//import public key
	import_pk(pk, H_pub);
	
	vec_mat_prod(syndrome_by_e, H_pub, errorMtx);

	for(i=0; i<CODE_N-CODE_K; ++i)
		if(get_element(syndrome_by_hash, 0, i) != get_element(syndrome_by_e, 0, i))
			return VERIF_REJECT;

	memcpy(m, m_rx, mlen_rx);
	*mlen = mlen_rx;

	delete_matrix(errorMtx);
	delete_matrix(H_pub);

	delete_matrix(syndrome_by_hash);
	delete_matrix(syndrome_by_e);
	free(m_rx);

	return 0;
}
