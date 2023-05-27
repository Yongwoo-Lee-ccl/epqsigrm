#include "api.h"
#include "common.h"

void import_pk(const unsigned char *pk, matrix *H_pub){
    import_matrix(H_pub, pk);
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
        return VERIF_REJECT;
    }
    
    //import public key
    import_pk(pk, Hpub);

    uint8_t randstr[syndrome_by_hash->ncols/8 + 1];
    hash_message(randstr, m_rx, mlen_rx, sign_i);
    randomize(syndrome_by_hash, randstr);

    vec_mat_prod(syndrome_by_e, Hpub, sign);

    for(uint32_t i=0; i < CODE_N-CODE_K - 1; ++i){
        if(get_element(syndrome_by_hash, 0, i) != get_element(syndrome_by_e, 0, i)){
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
