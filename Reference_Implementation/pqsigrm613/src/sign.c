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

void import_sk(const unsigned char *sk, uint16_t **Q, uint16_t **part_perm1, uint16_t **part_perm2, matrix* Hrep)
{
    *Q             = (uint16_t*)(sk);
    *part_perm1 = (uint16_t*)(sk+sizeof(uint16_t)*CODE_N);
    *part_perm2 = (uint16_t*)(sk+sizeof(uint16_t)*CODE_N + sizeof(uint16_t)*CODE_N/4);
    import_matrix(Hrep, sk+sizeof(uint16_t)*CODE_N + (sizeof(uint16_t)*CODE_N/4)*2);
}

// yr and yc are equal at the end.
void y_init(float *yc, float *yr, matrix* syndrome, uint16_t *Q){
    for(uint32_t i=0; i < syndrome->ncols; i++) {
        uint8_t bit = get_element(syndrome, 0, i);
        yc[i] = (bit == 0)? 1.:-1.;
    }
    for(uint32_t i = syndrome->ncols; i < CODE_N; i++) {
        yc[i] = 1.;
    }

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
    uint16_t *Q, *part_perm1, *part_perm2, *s_lead;
    matrix* Hrep = new_matrix(K_REP, (1<<RM_R));
    // matrix* Sinv = new_matrix(CODE_N - CODE_K - 1, CODE_N - CODE_K - 1);

    // import_sk(sk, &Q, &part_perm1, &part_perm2, Hrep, Sinv);
    import_sk(sk, &Q, &part_perm1, &part_perm2, Hrep);
    
    // Do signing, decode until the a error vector wt <= w is achieved
    uint64_t sign_i;
    matrix *challenge= new_matrix(1, CODE_N - CODE_K - 1);
    matrix *syndrome = new_matrix(1, CODE_N - CODE_K - 1);

    float yc[CODE_N];
    float yr[CODE_N];
    float mempool[CODE_N];

    init_decoding(mempool);
    uint32_t iter = 0;
    uint8_t randstr[challenge->ncols/8];

    while(1){
        // random number
        randombytes((unsigned char*)&sign_i, sizeof(uint64_t));
        hash_message(randstr, m, mlen, sign_i);
        randomize(challenge, randstr);
        // Find syndrome Sinv * challenge
        // vec_mat_prod(syndrome, Sinv, challenge);
        syndrome = challenge;
        // printf("syndrome:\n");
        // for (uint32_t i = 0; i < syndrome->ncols; i++)
        // {
        //     printf("%d", get_element(syndrome, 0, i));
        // }printf("\n");

        y_init(yc, yr, syndrome, Q);
        
        // decode and find e
        // In the recursive decoding procedure,
        // Y is 1 when the received codeword is 0, o.w, -1
        recursive_decoding_mod(yc, RM_R, RM_M, 0, CODE_N, part_perm1, part_perm2, Hrep);
        
        // Check Hamming weight of e'
        if(wgt(yr, yc) <= WEIGHT_PUB) break;
    }

    // compute Qinv*e'
    matrix *sign = new_matrix(1, CODE_N);
    for(uint32_t i=0; i < CODE_N; i++){
        set_element(sign, 0, i, (uint8_t)(yr[Q[i]] != yc[Q[i]]));
        // set_element(sign, 0, i, (yc[Q[i]] >= 0)?0 : 1);
    }

    {// Decoding: verification
        // generate RM code
        matrix* Gm = new_matrix(CODE_K, CODE_N);
        rm_gen(Gm, RM_R, RM_M, 0, CODE_K, 0, CODE_N);

        // partial replacement
        matrix* Grep = new_matrix((1<<RM_R) - K_REP, (1<<RM_R));
        dual(Hrep, Grep);
        //print_matrix_sign(Hrep, 0, Hrep->nrows, 0, Hrep->ncols);
        //printf("\n^Hrep\n");
        //print_matrix_sign(Grep, 0, Grep->nrows, 0, Grep->ncols);
        //printf("\n^grep\n");
        for (uint32_t i = 0; i < CODE_N; i += Grep->ncols)
        {
            partial_replace(Gm, K_REP, K_REP + Grep->nrows, i, i + Grep->ncols, Grep, 0, 0); 
        }

        // partial permutation
        for (uint32_t i = 0; i < 4; ++i)
        {
            col_permute(Gm, 0, rm_dim[RM_R][RM_M -2], 
                i*(CODE_N/4),(i+1)*(CODE_N/4), part_perm1);
        }
        col_permute(Gm, CODE_K - rm_dim[RM_R-2][RM_M-2], CODE_K, 3*CODE_N/4, CODE_N, part_perm2);
        /*
        for (size_t i = 0; i < CODE_N/4 ; i++)
        {
            printf("%d ", part_perm1[i]);
        }printf("\n^ partperm1\n");
        for (size_t i = 0; i < CODE_N/4 ; i++)
        {
            printf("%d ", part_perm2[i]);
        }printf("\n^ partperm2\n");
        */

        // find dual 
        // don't add a codeword from dual
        matrix* Hm = new_matrix(CODE_N-CODE_K, CODE_N);
        dual(Gm, Hm);
        //print_matrix_sign(Hm, 1000, 1064, 2000, 2064);


        // a codeword
        matrix* codeword = new_matrix(1, CODE_N);
        for(size_t i = 0; i < codeword->ncols; i++){
            set_element(codeword, 0, i, (yc[Q[i]] >= 0)?0 : 1);
        }
        col_permute(Hm, 0, Hm->nrows, 0, Hm->ncols, Q);
        rref(Hm, NULL);
        //print_matrix_sign(Hm, 1000, 1064, 2000, 2064);

        matrix* syndrome_test = new_matrix(1, Hm->nrows);
        vec_mat_prod(syndrome_test, Hm, codeword);

        
        // mat_mat_add(syndrome_test, syndrome, syndrome_challenge);
        uint8_t res = is_zero(syndrome_test);
        //printf("is decoding ok?: %d\n", res);
    }// Decoding: verification


    // export message
    // sign is (mlen, M, e, sign_i)
    // M includes its length, i.e., mlen
    *(unsigned long long*)sm = mlen;
    memcpy(sm+sizeof(unsigned long long), m, mlen);
    export_matrix(sign, sm+sizeof(unsigned long long)+mlen);

    *(unsigned long long*)(sm + sizeof(unsigned long long) + mlen + sign->nrows * sign->ncols/8) 
        = sign_i;

    *smlen = sizeof(unsigned long long) + mlen + sign->nrows * sign->ncols/8 + sizeof(unsigned long long);
    
    delete_matrix(Hrep);
    delete_matrix(challenge);
    delete_matrix(sign);

// 추가시작
//    delete_matrix(syndrome); 
//    delete_matrix(Gm);          //왜 오류??
//    delete_matrix(Grep);        //왜 오류??
//    delete_matrix(Hm);          //왜 오류??
//    delete_matrix(codeword);
//    delete_matrix(syndrome_test);//왜 오류??
// 추가끝

    return 0;    
}
