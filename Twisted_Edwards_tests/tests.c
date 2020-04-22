#include <stdio.h>
#include "ecurve.h"
#include <assert.h>
#include "tests.h"
#include "ak_curves.h"
#include "esign.h"
#include <stdlib.h>
#include <string.h>
#include <ak_parameters.h>
#include <time.h>
#include <math.h>
int cmp_epoint(ak_epoint ep1, ak_epoint ep2, ak_ecurve ec);
void test_add() {
	struct ecurve gost = gost_3410_2012_512_paramSetC;
	ak_ecurve ec = &gost;
	struct epoint ep1, ep2, ep3, ep4;
	ak_epoint epoint1 = &ep1, epoint2 = &ep2, epoint3 = &ep3, epoint4 = &ep4;

	ak_epoint_set_epoint(epoint1, &ec->point, ec);
	ak_epoint_set_as_unit(epoint2, ec);
	ak_epoint_set_as_unit(epoint3, ec);
	ak_epoint_set_epoint(epoint4, epoint1, ec);

	ak_epoint_add(epoint1, epoint2, ec);
	ak_epoint_reduce(epoint1, ec);

	assert(ak_epoint_is_ok(epoint1, ec) && "test add : P + 0 : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint1, epoint4, ec) && "test add : P + 0 : Something is wrong!\n");

	ak_epoint_add(epoint2, epoint2, ec);
	ak_epoint_reduce(epoint2, ec);

	assert(ak_epoint_is_ok(epoint2, ec) && "test add : 0 + 0 : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint2, epoint3, ec) && "test add : 0 + 0 : Something is wrong!\n");

}

void test_double() {
	struct ecurve gost = gost_3410_2012_512_paramSetC;
	ak_ecurve ec = &gost;
	struct epoint ep1, ep2;
	ak_epoint epoint1 = &ep1, epoint2 = &ep2;

	ak_epoint_set_as_unit(epoint1, ec);
	ak_epoint_set_as_unit(epoint2, ec);

	ak_epoint_double(epoint1, ec);
	ak_epoint_reduce(epoint1, ec);

	assert(ak_epoint_is_ok(epoint1, ec) && "test double : 2*0 = 0 : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint1, epoint2, ec) && "test double : 2*0 = 0 : Something is wrong!\n");

}
void test_pow(void (*f)( ak_epoint , ak_epoint, ak_uint64 *, size_t,
                                       ak_ecurve, ak_uint64, ak_uint64 ),ak_uint64 l, ak_uint64 w) {
	struct ecurve gost = gost_3410_2012_512_paramSetC;
	ak_ecurve ec = &gost;
	struct epoint ep1, ep2, ep3, ep4;
	ak_epoint epoint1 = &ep1, epoint2 = &ep2, epoint3 = &ep3, epoint4 = &ep4;

	ak_epoint_set_epoint(epoint1, &ec->point, ec);
	ak_epoint_set_as_unit(epoint2, ec);
	ak_epoint_set_epoint(epoint3, epoint1, ec);
	ak_epoint_set_epoint(epoint4, epoint1, ec);

    //ak_epoint_pow(epoint1, epoint1, ec->q, ec->size, ec);
    f(epoint1, epoint1, ec->q, ec->size, ec,l,w);
	ak_epoint_reduce(epoint1, ec);

	assert(ak_epoint_is_ok(epoint1, ec) && "test pow : qP=0 : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint1, epoint2, ec) && "test pow : qP=0 : Something is wrong!");

	ak_mpzn512 one = ak_mpzn512_one, two, sum;
	ak_mpzn_add_montgomery(one, one, ec->q, ec->p, ec->size);
    //ak_epoint_pow(epoint3, epoint4, one, ec->size, ec);
    f(epoint3, epoint4, one, ec->size, ec,l,w);
	ak_epoint_reduce(epoint3, ec);
	assert(ak_epoint_is_ok(epoint3, ec) && "test pow : (q+1)P=P : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint3, epoint4, ec) && "test pow : (q+1)P=P : Something is wrong!");

	struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 25; i++) {
		ak_epoint_set_epoint(epoint3, &ec->point, ec);
		ak_epoint_set_epoint(epoint4, &ec->point, ec);
		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        //ak_epoint_pow(epoint3, epoint3, one, ec->size, ec);
        f(epoint3, epoint3, one, ec->size, ec,l,w);
		ak_epoint_double(epoint3, ec);
		ak_epoint_reduce(epoint3, ec);
		ak_epoint_double(epoint4, ec);
        //ak_epoint_pow(epoint4, epoint4, one, ec->size, ec);
        f(epoint4, epoint4, one, ec->size, ec,l,w);
		ak_epoint_reduce(epoint4, ec);
		assert(ak_epoint_is_ok(epoint3, ec) && "test pow : 2((r1)P)=r1((2)P) : Point 1 is not on a curve!\n");
		assert(ak_epoint_is_ok(epoint4, ec) && "test pow : 2((r1)P)=r1((2)P) : Point 2 is not on a curve!\n");
		assert(!cmp_epoint(epoint3, epoint4, ec) && "test pow : 2((r1)P)=r1((2)P) : Something is wrong!");

		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
		ak_mpzn_set_random_modulo(two, ec->q, ec->size, &r);
		ak_epoint_set_epoint(epoint3, &ec->point, ec);
		ak_epoint_set_epoint(epoint4, &ec->point, ec);
		ak_mpzn_add_montgomery(sum, one, two, ec->p, ec->size);


        //ak_epoint_pow(epoint4, epoint4, sum, ec->size, ec);
        f(epoint4, epoint4, sum, ec->size, ec,l,w);
		ak_epoint_reduce(epoint4, ec);
        //ak_epoint_pow(epoint1, epoint3, one, ec->size, ec);
        //ak_epoint_pow(epoint2, epoint3, two, ec->size, ec);
        f(epoint1, epoint3, one, ec->size, ec,l,w);
        f(epoint2, epoint3, two, ec->size, ec,l,w);
        ak_epoint_add(epoint1, epoint2, ec);
		ak_epoint_reduce(epoint1, ec);

		assert(ak_epoint_is_ok(epoint1, ec) && "test pow : (r1+r2)P=r1P+r2P : Point 1 is not on a curve!\n");
		assert(ak_epoint_is_ok(epoint4, ec) && "test pow : (r1+r2)P=r1P+r2P : Point 2 is not on a curve!\n");
		assert(!cmp_epoint(epoint1, epoint4, ec) && "test pow : (r1+r2)P=r1P+r2P : Something is wrong!");
	}
	ak_random_context_destroy(&r);
}
void test_form_changing() {
	struct ecurve e_gost = gost_3410_2012_512_paramSetC; ak_ecurve ec = &e_gost;
	struct wcurve w_gost = id_tc26_gost_3410_2012_512_paramSetC; ak_wcurve wc = &w_gost;
    struct epoint ep1, ep2;
    ak_epoint epoint1 = &ep1, epoint2 = &ep2;
    struct wpoint wp1, wp2;
    ak_wpoint wpoint1 = &wp1, wpoint2 = &wp2;

	ak_epoint_set_epoint(epoint1, &ec->point, ec);
    ak_epoint_set_as_unit(epoint2, ec);

	ak_mpzn512 one=ak_mpzn512_one;
	
	struct random r;
    ak_random_context_create_lcg(&r);
	for (int i = 0; i <= 50; i++) {
		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
		ak_epoint_set_epoint(epoint1, &ec->point, ec);
		ak_wpoint_set_wpoint(wpoint1, &wc->point, wc);

		ak_epoint_pow(epoint1, epoint1, one, ec->size, ec);
		ak_wpoint_pow(wpoint1, wpoint1, one, wc->size, wc);
		
		ak_epoint_to_wpoint(epoint1, ec, wpoint2, wc);
		ak_wpoint_to_epoint(wpoint1, wc, epoint2, ec);
		ak_epoint_reduce(epoint2, ec);
		ak_wpoint_reduce(wpoint2, wc);
		ak_epoint_reduce(epoint1, ec);
		ak_wpoint_reduce(wpoint1, wc);
		assert(!ak_mpzn_cmp(wpoint1->x, wpoint2->x, ec->size) && "test form : wpoint->X : Something is wrong!");
		assert(!ak_mpzn_cmp(wpoint1->y, wpoint2->y, ec->size) && "test form : wpoint->Y : Something is wrong!");
		assert(!ak_mpzn_cmp(epoint1->X, epoint2->X, ec->size) && "test form : epoint->X : Something is wrong!");
		assert(!ak_mpzn_cmp(epoint1->Y, epoint2->Y, ec->size) && "test form : epoint->Y : Something is wrong!");
	}
	ak_random_context_destroy(&r);
}

int ak_signkey_test_edwards_() {

	ak_uint8 key512[64] = {
	 0xd4, 0x8d, 0xa1, 0x1f, 0x82, 0x67, 0x29, 0xc6, 0xdf, 0xaa, 0x18, 0xfd, 0x7b, 0x6b, 0x63, 0xa2,
	 0x14, 0x27, 0x7e, 0x82, 0xd2, 0xda, 0x22, 0x33, 0x56, 0xa0, 0x00, 0x22, 0x3b, 0x12, 0xe8, 0x72,
	 0x20, 0x10, 0x8b, 0x50, 0x8e, 0x50, 0xe7, 0x0e, 0x70, 0x69, 0x46, 0x51, 0xe8, 0xa0, 0x91, 0x30,
	 0xc9, 0xd7, 0x56, 0x77, 0xd4, 0x36, 0x09, 0xa4, 0x1b, 0x24, 0xae, 0xad, 0x8a, 0x04, 0xa6, 0x0b };

	ak_uint64 e512[ak_mpzn512_size] =
	{ 0xC6777D2972075B8CLL, 0x407ADEDB1D560C4FLL, 0x4339976C647C5D5ALL, 0x7184EE536593F441,
	  0xA71D147035B0C591LL, 0x1B09B6F9C170C533LL, 0x5C4F4A7C4D8DAB53LL, 0x3754F3CFACC9E061 };
	ak_uint64 k512[ak_mpzn512_size] =
	{ 0xA3AF71BB1AE679F1LL, 0x212273A6D14CF70ELL, 0x4434006011842286LL, 0x86748ED7A44B3E79LL,
	  0xD455986E364F3658LL, 0x946312120B39D019LL, 0xCC570456C6801496LL, 0x0359E7F4B1410FEALL };

#ifndef LIBAKRYPT_LITTLE_ENDIAN
	int i = 0;
#endif
	ak_uint64 e[64];
	char *str = NULL;
	struct signkey sk;
	ak_uint8 sign[128];
	struct verifykey pk;
	int error = ak_error_ok, audit = ak_log_get_level();
	struct ecurve ec = gost_3410_2012_512_paramSetC;
	if (audit >= ak_log_maximum)
		ak_error_message(ak_error_ok, __func__, "testing digital signatures started");
	/* 2. Âòîðîé ïðèìåð èç ïðèëîæåíèÿ À ÃÎÑÒ Ð 34.10-2012. */
	if ((error = ak_signkey_context_create_streebog512(&sk,
		(ak_wcurve)&id_tc26_gost_3410_2012_512_paramSetC)) != ak_error_ok) {
		ak_error_message(error, __func__,
			"incorrect creation of 512 bits secret key for GOST R 34.10-2012");
		return ak_false;
	}

	if ((error = ak_signkey_context_set_key(&sk, key512, sizeof(key512), ak_true)) != ak_error_ok) {
		ak_error_message(error, __func__, "incorrect assigning a constant key value");
		ak_signkey_context_destroy(&sk);
		return ak_false;
	}

	memset(sign, 0, 128);
	ak_signkey_context_sign_const_values_edwards(&sk, k512, e512, sign,&ec);

	if ((error = ak_verifykey_context_create_from_signkey(&pk, &sk)) != ak_error_ok) {
		ak_error_message(error, __func__, "incorrect creation of digital signature public key");
		ak_signkey_context_destroy(&sk);
		return ak_false;
	}
	ak_signkey_context_destroy(&sk);

	/* ïîñêîëüêó ïðîâåðêà âûïîëíÿåòñÿ äëÿ ðåçóëüòàòà ôóíêöèè õåøèðîâàíèÿ, òî ìû
	   ïðåîáðàçóåì ïîñëåäîâàòåëüíîñòü 64õ áèòíûõ èíòîâ â ïîñëåäîâàòåëüíîñòü áàéò */
#ifdef LIBAKRYPT_LITTLE_ENDIAN
	memcpy(e, e512, 64);
#else
	for (i = 0; i < ak_mpzn512_size; i++) ((ak_uint64 *)e)[i] = bswap_64(e512[i]);
#endif

	if (ak_verifykey_context_verify_hash_edwards(&pk, e, 64, sign,&ec)) {
		if (audit >= ak_log_maximum)
			ak_error_message(ak_error_ok, __func__,
				"digital signature verification process from GOST R 34.10-2012 (for 512 bit curve) is Ok");
	}
	else {
		ak_error_message(ak_error_not_equal_data, __func__,
			"digital signature verification process from GOST R 34.10-2012 (for 512 bit curve) is wrong");
		ak_log_set_message((str = ak_ptr_to_hexstr(sign, 64, ak_true))); free(str);
		ak_verifykey_context_destroy(&pk);
		return ak_false;
	}
	ak_verifykey_context_destroy(&pk);
	return ak_true;

}
void ak_signkey_test_edwards(){
    assert(ak_signkey_test_edwards_());
}

int cmp_epoint(ak_epoint ep1, ak_epoint ep2, ak_ecurve ec) {
	if (ak_mpzn_cmp(ep1->X, ep2->X, ec->size)) return 1;
	if (ak_mpzn_cmp(ep1->Y, ep2->Y, ec->size)) return 2;
	if (ak_mpzn_cmp(ep1->Z, ep2->Z, ec->size)) return 3;
	return 0;
}

void test_triple(){

    struct ecurve gost = gost_3410_2012_512_paramSetC;
    ak_ecurve ec = &gost;
    struct epoint ep1, ep2;
    ak_epoint epoint1 = &ep1, epoint2 = &ep2;
    ak_mpzn512 one = ak_mpzn512_one;

    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 50; i++) {
        ak_epoint_set(epoint1, ec);
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_epoint_pow(epoint1, epoint1, one, ec->size, ec);
        ak_epoint_set_epoint(epoint2,epoint1,ec);
        ak_epoint_triple(epoint1, ec);
        ak_mpzn_set_ui(one,ak_mpzn512_size,3);
        ak_epoint_pow(epoint2, epoint2, one, ec->size, ec);
        ak_epoint_reduce(epoint1, ec);
        ak_epoint_reduce(epoint2, ec);
        assert(ak_epoint_is_ok(epoint1, ec) && "test trip : Point 1 is not on a curve!\n");
        assert(ak_epoint_is_ok(epoint2, ec) && "test trip : Point 2 is not on a curve!\n");
        assert(!cmp_epoint(epoint1, epoint2, ec) && "test trip : Something is wrong!");

    }
    ak_random_context_destroy(&r);
}

void test_quintuple(){

    struct ecurve gost = gost_3410_2012_512_paramSetC;
    ak_ecurve ec = &gost;
    struct epoint ep1, ep2;
    ak_epoint epoint1 = &ep1, epoint2 = &ep2;
    ak_mpzn512 one = ak_mpzn512_one;

    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 50; i++) {
        ak_epoint_set(epoint1, ec);
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_epoint_pow(epoint1, epoint1, one, ec->size, ec);
        ak_epoint_set_epoint(epoint2,epoint1,ec);
        ak_epoint_quintuple(epoint1, ec);
        ak_mpzn_set_ui(one,ak_mpzn512_size,5);
        ak_epoint_pow(epoint2, epoint2, one, ec->size, ec);
        ak_epoint_reduce(epoint1, ec);
        ak_epoint_reduce(epoint2, ec);
        assert(ak_epoint_is_ok(epoint1, ec) && "test trip : Point 1 is not on a curve!\n");
        assert(ak_epoint_is_ok(epoint2, ec) && "test trip : Point 2 is not on a curve!\n");
        assert(!cmp_epoint(epoint1, epoint2, ec) && "test trip : Something is wrong!");

    }
    ak_random_context_destroy(&r);
}
void M_pow_mont( ak_epoint wq, ak_epoint ep, ak_uint64 *k, size_t size, ak_ecurve ec, ak_uint64 l, ak_uint64 w){
    ak_epoint_pow(wq,ep,k,size,ec);
}
void test_pow_montgomery(){
    test_pow(M_pow_mont,0,0);
}
void M_pow_bin( ak_epoint wq, ak_epoint ep, ak_uint64 *k, size_t size, ak_ecurve ec, ak_uint64 l, ak_uint64 w){
    ak_epoint_pow_binary(wq,ep,k,size,ec);
}
void test_pow_bin(){
    test_pow(M_pow_bin,0,0);
}
void M_pow_NAF( ak_epoint wq, ak_epoint ep, ak_uint64 *k, size_t size, ak_ecurve ec, ak_uint64 l, ak_uint64 w){
    ak_epoint_pow_NAF(wq,ep,k,size,ec);
}

void test_pow_NAF(){
    test_pow(M_pow_NAF,0,0);
}
void M_pow_NAF_2_w( ak_epoint wq, ak_epoint ep, ak_uint64 *k, size_t size, ak_ecurve ec, ak_uint64 l, ak_uint64 w){
    ak_epoint_pow_NAF_powof2(wq,ep,k,size,w,ec);
}
void test_pow_NAF_2_w(){
    for(int i=2;i<5;i++)test_pow(M_pow_NAF_2_w,0,i);
}


void M_pow_NAF_l_w( ak_epoint wq, ak_epoint ep, ak_uint64 *k, size_t size, ak_ecurve ec, ak_uint64 l, ak_uint64 w){
    ak_epoint_pow_NAF_powofL(wq,ep,k,size,l,w,ec);
}
void test_pow_NAF_l_w(){
    for(int i=2;i<5;i++) test_pow(M_pow_NAF_l_w,3,i);
    for(int i=2;i<5;i++) test_pow(M_pow_NAF_l_w,5,i);
}


void test_pow_NAF_mbw(){
    struct ecurve gost = gost_3410_2012_512_paramSetC;
    ak_ecurve ec = &gost;
    struct epoint ep1, ep2;
    ak_epoint epoint1 = &ep1, epoint2 = &ep2;
    const int ss=2;
    ak_uint32 l[ss]={2,3};
    ak_uint32 w[ss]={1,2};
    ak_epoint_set_epoint(epoint1, &ec->point, ec);

    ak_mpzn512 one;
    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 25; i++) {
        ak_epoint_set_epoint(epoint1, &ec->point, ec);
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_epoint_pow(epoint2, epoint1, one, ec->size, ec);
        ak_epoint_pow_NAF_mbw(epoint1, epoint1, one, ec->size,l,w,ss,ec);
        ak_epoint_reduce(epoint1, ec);
        ak_epoint_reduce(epoint2, ec);
        assert(!cmp_epoint(epoint1, epoint2, ec) && "MBW pow : Something is wrong!");

    }
    ak_random_context_destroy(&r);
}

void get_time(){

    struct ecurve gost = gost_3410_2012_512_paramSetC;
    ak_ecurve ec = &gost;
    struct epoint ep1; ak_epoint epoint1 = &ep1;
    ak_epoint_set_epoint(epoint1, &ec->point, ec);
    ak_mpzn512 one = ak_mpzn512_one;
    struct random r;
    ak_random_context_create_lcg(&r);
    ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);

    printf("\nЛесенка Монтгомери");
    ak_epoint_set_epoint(epoint1, &ec->point, ec);
    for (int i = 0; i <= 10; i++) {
        clock_t cl = clock( );
        ak_epoint_pow(epoint1, epoint1, one, ec->size, ec);
        printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
    }


    printf("\nБинарный");
    ak_epoint_set_epoint(epoint1, &ec->point, ec);
    for (int i = 0; i <= 10; i++) {
        clock_t cl = clock( );
        ak_epoint_pow_binary(epoint1, epoint1, one, ec->size, ec);
        printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
    }


    printf("\nNAF 2^2 (отдельный алгоритм) с дополнительными вычислениями");
    ak_epoint_set_epoint(epoint1, &ec->point, ec);
    for (int i = 0; i <= 10; i++) {
        clock_t cl = clock( );
        ak_epoint_pow_NAF(epoint1, epoint1, one, ec->size, ec);
        printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
    }

    for(int j=2;j<7;j++){
        printf("\nNAF 2^%d с дополнительными вычислениями",j);
        ak_epoint_set_epoint(epoint1, &ec->point, ec);
        for (int i = 0; i <= 10; i++) {
            clock_t cl = clock( );
            ak_epoint_pow_NAF_powof2(epoint1, epoint1, one, ec->size, j, ec);
            printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
        }
    }

    for(int base=2;base<5;base++){
        for(int j=2;j<7;j++){
            printf("\nNAF L=%d^%d с дополнительными вычислениями",base,j);
            ak_epoint_set_epoint(epoint1, &ec->point, ec);
            for (int i = 0; i <= 10; i++) {
                clock_t cl = clock( );
                ak_epoint_pow_NAF_powofL(epoint1, epoint1, one, ec->size, base,j, ec);
                printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
            }
        }
        base=(base==3)?4:base;
    }

    for(int j=2;j<6;j++){
        printf("\nNAF с основаниями: 2,3 и одинаковым окном %d с дополнительными вычислениями",j);
        ak_epoint_set_epoint(epoint1, &ec->point, ec);
        ak_uint32 bases[2]={2,3};
        ak_uint32 wind[2]={j,j};
        for (int i = 0; i <= 10; i++) {
            clock_t cl = clock( );
            ak_epoint_pow_NAF_mbw(epoint1, epoint1, one, ec->size, bases,wind, 2, ec);
            printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
        }
    }
    for(int j=2;j<4;j++){
        printf("\nNAF с основаниями: 2,3,5 и одинаковым окном %d с дополнительными вычислениями",j);
        ak_epoint_set_epoint(epoint1, &ec->point, ec);
        ak_uint32 bases[3]={2,3,5};
        ak_uint32 wind[3]={j,j,j};
        for (int i = 0; i <= 10; i++) {
            clock_t cl = clock( );
            ak_epoint_pow_NAF_mbw(epoint1, epoint1, one, ec->size, bases,wind, 3, ec);
            printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
        }
    }
    ak_mpznmax k;
    ak_mpzn_set(k,one,ec->size);
    size_t size=ec->size;
    printf("\nNAF 2^2 (отдельный алгоритм) все значения предвычислены");
    for (int i = 0; i <= 10; i++) {
        ak_epoint_set_epoint(epoint1, &ec->point, ec);
        ak_int32 kNAF[512];
        ak_int32 i=ak_n_to_NAF2(k,kNAF,ec->size);
        struct epoint Q, R,R_;
        ak_epoint_set_as_unit( &Q, ec );
        ak_epoint_set_epoint( &R, epoint1, ec );
        ak_epoint_set_epoint( &R_, epoint1, ec );
        ak_invert_ep(&R_,&R,ec);
        clock_t cl = clock();
        for (i=i-1; i>=0; i--) {
            ak_epoint_double(&Q,ec);
            ak_int32 nj=kNAF[i];
            if (nj>0) {
                ak_epoint_add(&Q,&R,ec);
            }
            else if(nj<0){
                ak_epoint_add(&Q,&R_,ec);
            }
        }
        ak_epoint_set_epoint( epoint1, &Q, ec );
        printf(";%f",(double)((clock()-cl)) / (double)CLOCKS_PER_SEC);
   }
   for(int w=2;w<7;w++){
       printf("\nNAF 2^%d все значения предвычислены",w);

       for (int i = 0; i <= 10; i++) {
           ak_epoint_set_epoint(epoint1, &ec->point, ec);
           ak_mpznmax pows;
           ak_int32 kNAF[512];
           ak_int32 i=ak_n_to_NAF_powof2(k,w,kNAF,ec->size);
           struct epoint Q, R_[32],inv;
           ak_epoint_set_as_unit(&Q, ec );
           ak_epoint_set_epoint(R_, epoint1, ec );
           ak_epoint_set_epoint(R_+1, epoint1, ec );
           ak_epoint_set_epoint(R_+2, epoint1, ec );
           ak_epoint_triple(R_+1, ec);
           if(w>3){
               ak_epoint_quintuple(R_+2, ec);
               for (ak_int32 j=3;j<w;j++){
                   ak_mpzn_set_ui(pows,size,j*2+1);
                   ak_epoint_pow_binary(R_+j,R_,pows,size,ec);
               }
           }
           clock_t cl = clock( );
           for (i=i-1; i>=0; i--) {
               ak_epoint_double(&Q,ec);
               ak_int32 nj=kNAF[i];
               if (nj>0) {
                   ak_epoint_add(&Q,R_+((nj-1)>>1),ec);
               }
               else if(nj<0){
                   ak_invert_ep(&inv,R_+((-nj-1)>>1),ec);
                   ak_epoint_add(&Q,&inv,ec);
               }
           }
           ak_epoint_set_epoint( epoint1, &Q, ec );
           printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);

       }

   }

   for(int base=2;base<5;base++){
       for(int w=2;w<7;w++){
           printf("\nNAF L=%d^%d все значения предвычислены",base,w);
           ak_epoint_set_epoint(epoint1, &ec->point, ec);
           for (int i = 0; i <= 10; i++) {
               int l =base;
               ak_mpznmax pows;
               ak_int32 kNAF[512],LL_=pow(l,w-1);
               ak_int32 len=LL_*(l-1);

               ak_int32 i=ak_n_to_NAF_L_w(k,kNAF,l,w,size,ec);
               struct epoint Q, R_[512],inv;
               ak_epoint_set_as_unit(&Q, ec );
               ak_epoint_set_epoint(R_, epoint1, ec );
               ak_int32 cr=1;
               for (ak_int32 j=2;j<len;j++){
                   if(j%l){ ak_mpzn_set_ui(pows,size,j);
                       ak_epoint_pow_binary(R_+cr,R_,pows,size,ec);
                       cr++;
                   }
               }


               clock_t cl = clock( );
               for (i=i-1; i>=0; i--) {
                   if(l==2) ak_epoint_double(&Q,ec);
                   else{
                       if(l==3) ak_epoint_triple(&Q,ec);
                       else if(l==5) ak_epoint_quintuple(&Q,ec);
                       else{
                           ak_mpzn_set_ui(pows,size,l);
                           ak_epoint_pow_binary(&Q,&Q,pows,size,ec);
                       }
                   }

                   ak_int32 nj=kNAF[i];
                   if (nj>0) {
                       ak_epoint_add(&Q,R_+(nj-1-(nj/l)),ec);
                   }
                   else if(nj<0){
                       ak_invert_ep(&inv,R_+((-nj)-1-(-nj)/l),ec);
                       ak_epoint_add(&Q,&inv,ec);
                   }
               }
               ak_epoint_set_epoint( epoint1, &Q, ec );
               printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
           }
       }
       base=(base==3)?4:base;
   }
   for(int j=2;j<6;j++){
       printf("\nNAF с основаниями: 2,3 и одинаковым окном %d все значения предвычислены",j);

       ak_uint32 l[2]={2,3};
       ak_uint32 w[2]={j,j};
       ak_int32 len=2;
       for (int i = 0; i <= 10; i++) {
           ak_epoint_set_epoint(epoint1, &ec->point, ec);
           ak_mpznmax pows;
           ak_int32 kNAF[512],bases[512];
           ak_int32 A=1;
           for(int j=0;j<len;j++) A*=pow(l[j],w[j]);
           ak_int32 i=ak_n_to_Ext_wmb_NAF(k,kNAF,bases,l,w,len,size,ec);
           struct epoint Q, R_[13500],inv;
           ak_epoint_set_as_unit(&Q, ec );
           ak_epoint_set_epoint(R_, epoint1, ec );
           for (ak_int32 j=2;j<A/2;j++){
               ak_mpzn_set_ui(pows,size,j);
               ak_epoint_pow_binary(R_+j-1,R_,pows,size,ec);
           }

           clock_t cl = clock( );
           for (i=i-1; i>=0; i--) {

               switch (bases[i]) {
               case 2:
                   ak_epoint_double(&Q,ec); break;
               case 3:
                   ak_epoint_triple(&Q,ec); break;
               case 5:
                   ak_epoint_quintuple(&Q,ec); break;
               default:
                   ak_mpzn_set_ui(pows,size,bases[i]);
                   ak_epoint_pow_binary(&Q,&Q,pows,size,ec);
                   break;
               }

               ak_int32 nj=kNAF[i];
               if (nj>0) {
                   ak_epoint_add(&Q,R_+(nj-1),ec);
               }
               else if(nj<0){
                   ak_invert_ep(&inv,R_+(-nj-1),ec);
                   ak_epoint_add(&Q,&inv,ec);
               }
           }

           ak_epoint_set_epoint( epoint1, &Q, ec );
           printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
       }
   }
   for(int j=2;j<4;j++){
       printf("\nNAF с основаниями: 2,3,5 и одинаковым окном %d все значения предвычислены",j);

       ak_uint32 l[3]={2,3,5};
       ak_uint32 w[3]={j,j,j};
       ak_int32 len=3;
       for (int i = 0; i <= 10; i++) {
           ak_epoint_set_epoint(epoint1, &ec->point, ec);
           ak_mpznmax pows;
           ak_int32 kNAF[512],bases[512];
           ak_int32 A=1;
           for(int j=0;j<len;j++) A*=pow(l[j],w[j]);
           ak_int32 i=ak_n_to_Ext_wmb_NAF(k,kNAF,bases,l,w,len,size,ec);
           struct epoint Q, R_[13500],inv;
           ak_epoint_set_as_unit(&Q, ec );
           ak_epoint_set_epoint(R_, epoint1, ec );
           for (ak_int32 j=2;j<A/2;j++){
               ak_mpzn_set_ui(pows,size,j);
               ak_epoint_pow_binary(R_+j-1,R_,pows,size,ec);
           }
           clock_t cl = clock( );
           for (i=i-1; i>=0; i--) {

               switch (bases[i]) {
               case 2:
                   ak_epoint_double(&Q,ec); break;
               case 3:
                   ak_epoint_triple(&Q,ec); break;
               case 5:
                   ak_epoint_quintuple(&Q,ec); break;
               default:
                   ak_mpzn_set_ui(pows,size,bases[i]);
                   ak_epoint_pow_binary(&Q,&Q,pows,size,ec);
                   break;
               }

               ak_int32 nj=kNAF[i];
               if (nj>0) {
                   ak_epoint_add(&Q,R_+(nj-1),ec);
               }
               else if(nj<0){
                   ak_invert_ep(&inv,R_+(-nj-1),ec);
                   ak_epoint_add(&Q,&inv,ec);
               }
           }

           ak_epoint_set_epoint( epoint1, &Q, ec );
           printf(";%f",(double)(clock()-cl) / (double)CLOCKS_PER_SEC);
       }
   }



    ak_random_context_destroy(&r);
}
