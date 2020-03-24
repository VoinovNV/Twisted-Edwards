#include <stdio.h>
#include "ecurve.h"
#include <assert.h>
#include "tests.h"
#include "ak_curves.h"
#include "esign.h"
#include <stdlib.h>
#include <string.h>
#include <ak_parameters.h>
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
void test_pow() {
	struct ecurve gost = gost_3410_2012_512_paramSetC;
	ak_ecurve ec = &gost;
	struct epoint ep1, ep2, ep3, ep4;
	ak_epoint epoint1 = &ep1, epoint2 = &ep2, epoint3 = &ep3, epoint4 = &ep4;

	ak_epoint_set_epoint(epoint1, &ec->point, ec);
	ak_epoint_set_as_unit(epoint2, ec);
	ak_epoint_set_epoint(epoint3, epoint1, ec);
	ak_epoint_set_epoint(epoint4, epoint1, ec);

	ak_epoint_pow(epoint1, epoint1, ec->q, ec->size, ec);
	ak_epoint_reduce(epoint1, ec);

	assert(ak_epoint_is_ok(epoint1, ec) && "test pow : qP=0 : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint1, epoint2, ec) && "test pow : qP=0 : Something is wrong!");

	ak_mpzn512 one = ak_mpzn512_one, two, sum;
	ak_mpzn_add_montgomery(one, one, ec->q, ec->p, ec->size);
	ak_epoint_pow(epoint3, epoint4, one, ec->size, ec);
	ak_epoint_reduce(epoint3, ec);
	assert(ak_epoint_is_ok(epoint3, ec) && "test pow : (q+1)P=P : Point is not on a curve!\n");
	assert(!cmp_epoint(epoint3, epoint4, ec) && "test pow : (q+1)P=P : Something is wrong!");

	struct random r;
	int abc = ak_random_context_create_lcg(&r);
	for (int i = 0; i <= 50; i++) {
		ak_epoint_set_epoint(epoint3, &ec->point, ec);
		ak_epoint_set_epoint(epoint4, &ec->point, ec);
		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
		ak_epoint_pow(epoint3, epoint3, one, ec->size, ec);
		ak_epoint_double(epoint3, ec);
		ak_epoint_reduce(epoint3, ec);
		ak_epoint_double(epoint4, ec);
		ak_epoint_pow(epoint4, epoint4, one, ec->size, ec);
		ak_epoint_reduce(epoint4, ec);
		assert(ak_epoint_is_ok(epoint3, ec) && "test pow : 2((r1)P)=r1((2)P) : Point 1 is not on a curve!\n");
		assert(ak_epoint_is_ok(epoint4, ec) && "test pow : 2((r1)P)=r1((2)P) : Point 2 is not on a curve!\n");
		assert(!cmp_epoint(epoint3, epoint4, ec) && "test pow : 2((r1)P)=r1((2)P) : Something is wrong!");

		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
		ak_mpzn_set_random_modulo(two, ec->q, ec->size, &r);
		ak_epoint_set_epoint(epoint3, &ec->point, ec);
		ak_epoint_set_epoint(epoint4, &ec->point, ec);
		ak_mpzn_add_montgomery(sum, one, two, ec->p, ec->size);


		ak_epoint_pow(epoint4, epoint4, sum, ec->size, ec);
		ak_epoint_reduce(epoint4, ec);
		ak_epoint_pow(epoint1, epoint3, one, ec->size, ec);
		ak_epoint_pow(epoint2, epoint3, two, ec->size, ec);
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
	struct epoint ep1, ep2, ep3, ep4;
	ak_epoint epoint1 = &ep1, epoint2 = &ep2, epoint3 = &ep3, epoint4 = &ep4;
	struct wpoint wp1, wp2, wp3, wp4;
	ak_wpoint wpoint1 = &wp1, wpoint2 = &wp2, wpoint3 = &wp3, wpoint4 = &wp4;

	ak_epoint_set_epoint(epoint1, &ec->point, ec);
	ak_epoint_set_as_unit(epoint2, ec);
	ak_epoint_set_epoint(epoint3, epoint1, ec);
	ak_epoint_set_epoint(epoint4, epoint1, ec);

	ak_mpzn512 one=ak_mpzn512_one;
	
	struct random r;
	int abc = ak_random_context_create_lcg(&r);
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

int ak_signkey_test_edwards() {

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
	ak_uint8 sign512[128] =
	{ 0x36, 0xAE, 0x73, 0xE1, 0x44, 0x93, 0xE1, 0x17, 0x33, 0x5C, 0x9C, 0xCD, 0xCB, 0x3B, 0xC9, 0x60,
	  0x02, 0x85, 0x99, 0x06, 0xC9, 0x97, 0xC1, 0x9E, 0x1C, 0x0F, 0xB2, 0x86, 0x84, 0x55, 0x92, 0x54,
	  0xD3, 0xAC, 0xFC, 0xA8, 0xEE, 0x78, 0x3C, 0x64, 0xC2, 0xDC, 0xE0, 0x2E, 0xC8, 0xA3, 0x12, 0xE5,
	  0x9E, 0x68, 0x3C, 0x1E, 0x5E, 0x79, 0xDD, 0x23, 0x1A, 0x09, 0x81, 0xA0, 0x60, 0xFA, 0x86, 0x2F,
	  0x4A, 0x5B, 0x3E, 0xE7, 0xBD, 0x53, 0x98, 0x2A, 0xB9, 0x9C, 0x91, 0x56, 0x1F, 0xEB, 0x6E, 0x6A,
	  0x40, 0xCE, 0x70, 0x7F, 0xDF, 0x80, 0x60, 0x52, 0x62, 0xF3, 0xC4, 0xE8, 0x88, 0xE2, 0x3C, 0x82,
	  0xF5, 0x2F, 0xD5, 0x33, 0xE9, 0xFB, 0x0B, 0x1C, 0x08, 0xBC, 0xAD, 0x8A, 0x77, 0x56, 0x5F, 0x32,
	  0xB6, 0x26, 0x2D, 0x36, 0xA9, 0xE7, 0x85, 0x65, 0x8E, 0xFE, 0x6F, 0x69, 0x94, 0xB3, 0x81, 0x10 };

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
	/* 2. Второй пример из приложения А ГОСТ Р 34.10-2012. */
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

	/* поскольку проверка выполняется для результата функции хеширования, то мы
	   преобразуем последовательность 64х битных интов в последовательность байт */
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
int my_ak_signkey_test(void)
{
	
	/* d = BA6048AADAE241BA40936D47756D7C93091A0E8514669700EE7508E508B102072E8123B2200A0563322DAD2827E2714A2636B7BFD18AADFC62967821FA18DD4 */
	ak_uint8 key512[64] = {
	 0xd4, 0x8d, 0xa1, 0x1f, 0x82, 0x67, 0x29, 0xc6, 0xdf, 0xaa, 0x18, 0xfd, 0x7b, 0x6b, 0x63, 0xa2,
	 0x14, 0x27, 0x7e, 0x82, 0xd2, 0xda, 0x22, 0x33, 0x56, 0xa0, 0x00, 0x22, 0x3b, 0x12, 0xe8, 0x72,
	 0x20, 0x10, 0x8b, 0x50, 0x8e, 0x50, 0xe7, 0x0e, 0x70, 0x69, 0x46, 0x51, 0xe8, 0xa0, 0x91, 0x30,
	 0xc9, 0xd7, 0x56, 0x77, 0xd4, 0x36, 0x09, 0xa4, 0x1b, 0x24, 0xae, 0xad, 0x8a, 0x04, 0xa6, 0x0b };
/* е = 3754F3CFACC9E0615C4F4A7C4D8DAB531B09B6F9C170C533A71D147035B0C5917184EE536593F4414339976C647C5D5A407ADEDB1D560C4FC6777D2972075B8C */
	ak_uint64 e512[ak_mpzn512_size] =
	{ 0xC6777D2972075B8CLL, 0x407ADEDB1D560C4FLL, 0x4339976C647C5D5ALL, 0x7184EE536593F441,
	  0xA71D147035B0C591LL, 0x1B09B6F9C170C533LL, 0x5C4F4A7C4D8DAB53LL, 0x3754F3CFACC9E061 };
	/* k = 359E7F4B1410FEACC570456C6801496946312120B39D019D455986E364F365886748ED7A44B3E794434006011842286212273A6D14CF70EA3AF71BB1AE679F1 */
	ak_uint64 k512[ak_mpzn512_size] =
	{ 0xA3AF71BB1AE679F1LL, 0x212273A6D14CF70ELL, 0x4434006011842286LL, 0x86748ED7A44B3E79LL,
	  0xD455986E364F3658LL, 0x946312120B39D019LL, 0xCC570456C6801496LL, 0x0359E7F4B1410FEALL };

	ak_uint8 sign512[128] =
	{ 0x36, 0xAE, 0x73, 0xE1, 0x44, 0x93, 0xE1, 0x17, 0x33, 0x5C, 0x9C, 0xCD, 0xCB, 0x3B, 0xC9, 0x60,
	  0x02, 0x85, 0x99, 0x06, 0xC9, 0x97, 0xC1, 0x9E, 0x1C, 0x0F, 0xB2, 0x86, 0x84, 0x55, 0x92, 0x54,
	  0xD3, 0xAC, 0xFC, 0xA8, 0xEE, 0x78, 0x3C, 0x64, 0xC2, 0xDC, 0xE0, 0x2E, 0xC8, 0xA3, 0x12, 0xE5,
	  0x9E, 0x68, 0x3C, 0x1E, 0x5E, 0x79, 0xDD, 0x23, 0x1A, 0x09, 0x81, 0xA0, 0x60, 0xFA, 0x86, 0x2F,
	  0x4A, 0x5B, 0x3E, 0xE7, 0xBD, 0x53, 0x98, 0x2A, 0xB9, 0x9C, 0x91, 0x56, 0x1F, 0xEB, 0x6E, 0x6A,
	  0x40, 0xCE, 0x70, 0x7F, 0xDF, 0x80, 0x60, 0x52, 0x62, 0xF3, 0xC4, 0xE8, 0x88, 0xE2, 0x3C, 0x82,
	  0xF5, 0x2F, 0xD5, 0x33, 0xE9, 0xFB, 0x0B, 0x1C, 0x08, 0xBC, 0xAD, 0x8A, 0x77, 0x56, 0x5F, 0x32,
	  0xB6, 0x26, 0x2D, 0x36, 0xA9, 0xE7, 0x85, 0x65, 0x8E, 0xFE, 0x6F, 0x69, 0x94, 0xB3, 0x81, 0x10 };

#ifndef LIBAKRYPT_LITTLE_ENDIAN
	int i = 0;
#endif
	ak_uint64 e[64];
	char *str = NULL;
	struct signkey sk;
	ak_uint8 sign[128];
	struct verifykey pk;
	int error = ak_error_ok, audit = ak_log_get_level();

	if (audit >= ak_log_maximum)
		ak_error_message(ak_error_ok, __func__, "testing digital signatures started");

	/* 2. Второй пример из приложения А ГОСТ Р 34.10-2012. */
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
	ak_signkey_context_sign_const_values(&sk, k512, e512, sign);
	
	if ((error = ak_verifykey_context_create_from_signkey(&pk, &sk)) != ak_error_ok) {
		ak_error_message(error, __func__, "incorrect creation of digital signature public key");
		ak_signkey_context_destroy(&sk);
		return ak_false;
	}
	ak_signkey_context_destroy(&sk);

#ifdef LIBAKRYPT_LITTLE_ENDIAN
	memcpy(e, e512, 64);
#else
	for (i = 0; i < ak_mpzn512_size; i++) ((ak_uint64 *)e)[i] = bswap_64(e512[i]);
#endif

	if (ak_verifykey_context_verify_hash(&pk, e, 64, sign)) {
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
    int abc = ak_random_context_create_lcg(&r);
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
    int abc = ak_random_context_create_lcg(&r);
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
