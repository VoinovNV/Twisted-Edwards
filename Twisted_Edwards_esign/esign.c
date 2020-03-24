#include "esign.h"
#include <string.h>
void ak_signkey_context_sign_const_values_edwards(ak_signkey sctx,
	ak_uint64 *k, ak_uint64 *e, ak_pointer out, ak_ecurve ec)
{
#ifndef LIBAKRYPT_LITTLE_ENDIAN
	int i = 0;
#endif
	struct wpoint wr; struct epoint er;
	ak_wcurve wc = (ak_wcurve)sctx->key.data;
	ak_uint64 *r = (ak_uint64 *)out, *s = (ak_uint64 *)out + wc->size;

	/* поскольку функция не экспортируется, мы оставляем все проверки функциям верхнего уровня */
	/* вычисляем r */
	ak_epoint_pow(&er, &ec->point, k, ec->size, ec);
	ak_epoint_to_wpoint(&er, ec, &wr, wc);
	ak_wpoint_reduce(&wr, wc);
	ak_mpzn_rem(r, wr.x, wc->q, wc->size);

	/* приводим r к виду Монтгомери и помещаем во временную переменную wr.x <- r */
	ak_mpzn_mul_montgomery(wr.x, r, wc->r2q, wc->q, wc->nq, wc->size);

	/* вычисляем значение s <- r*d (mod q) (сначала домножаем на ключ, потом на его маску) */
	ak_mpzn_mul_montgomery(s, wr.x, sctx->key.key.data, wc->q, wc->nq, wc->size);
	ak_mpzn_mul_montgomery(s, s, sctx->key.mask.data, wc->q, wc->nq, wc->size);

	/* приводим k к виду Монтгомери и помещаем во временную переменную wr.y <- k */
	ak_mpzn_mul_montgomery(wr.y, k, wc->r2q, wc->q, wc->nq, wc->size);

	/* приводим e к виду Монтгомери и помещаем во временную переменную wr.z <- e */
	ak_mpzn_rem(wr.z, e, wc->q, wc->size);
	if (ak_mpzn_cmp_ui(wr.z, wc->size, 0)) ak_mpzn_set_ui(wr.z, wc->size, 1);
	ak_mpzn_mul_montgomery(wr.z, wr.z, wc->r2q, wc->q, wc->nq, wc->size);

	/* вычисляем k*e (mod q) и вычисляем s = r*d + k*e (mod q) (в форме Монтгомери) */
	ak_mpzn_mul_montgomery(wr.y, wr.y, wr.z, wc->q, wc->nq, wc->size); /* wr.y <- k*e */
	ak_mpzn_add_montgomery(s, s, wr.y, wc->q, wc->size);

	/* приводим s к обычной форме */
	ak_mpzn_mul_montgomery(s, s, wc->point.z, /* для экономии памяти пользуемся равенством z = 1 */
		wc->q, wc->nq, wc->size);
#ifndef LIBAKRYPT_LITTLE_ENDIAN
	for (i = 0; i < 2 * wc->size; i++) ((ak_uint64*)out)[i] = bswap_64(((ak_uint64*)out)[i]);
#endif

	/* завершаемся */
	memset(&wr, 0, sizeof(struct wpoint));
	sctx->key.set_mask(&sctx->key);
}

bool_t ak_verifykey_context_verify_hash_edwards(ak_verifykey pctx,
	const ak_pointer hash, const size_t hsize, ak_pointer sign,ak_ecurve ec)
{
#ifndef LIBAKRYPT_LITTLE_ENDIAN
	int i = 0;
#endif
	ak_mpzn512 v, z1, z2, u, r, s, h;
	struct wpoint cpoint;

	if (pctx == NULL) {
		ak_error_message(ak_error_null_pointer, __func__,
			"using a null pointer to secret key context");
		return ak_false;
	}
	if (hash == NULL) {
		ak_error_message(ak_error_null_pointer, __func__, "using a null pointer to hash value");
		return ak_false;
	}
	if (hsize != sizeof(ak_uint64)*(pctx->wc->size)) {
		ak_error_message(ak_error_wrong_length, __func__, "using hash value with wrong length");
		return ak_false;
	}
	if (sign == NULL) {
		ak_error_message(ak_error_null_pointer, __func__, "using a null pointer to sign value");
		return ak_false;
	}

	memcpy(r, (ak_uint64*)sign, sizeof(ak_uint64)*pctx->wc->size);
	memcpy(s, (ak_uint64 *)sign + pctx->wc->size, sizeof(ak_uint64)*pctx->wc->size);
	memcpy(h, hash, sizeof(ak_uint64)*pctx->wc->size);

#ifndef LIBAKRYPT_LITTLE_ENDIAN
	for (i = 0; i < pctx->wc->size; i++) {
		r[i] = bswap_64(r[i]);
		s[i] = bswap_64(s[i]);
		h[i] = bswap_64(h[i]);
	}
#endif

	ak_mpzn_set(v, h, pctx->wc->size);
	ak_mpzn_rem(v, v, pctx->wc->q, pctx->wc->size);
	if (ak_mpzn_cmp_ui(v, pctx->wc->size, 0)) ak_mpzn_set_ui(v, pctx->wc->size, 1);
	ak_mpzn_mul_montgomery(v, v, pctx->wc->r2q, pctx->wc->q, pctx->wc->nq, pctx->wc->size);

	/* вычисляем v (в представлении Монтгомери) */
	ak_mpzn_set_ui(u, pctx->wc->size, 2);
	ak_mpzn_sub(u, pctx->wc->q, u, pctx->wc->size);
	ak_mpzn_modpow_montgomery(v, v, u, pctx->wc->q, pctx->wc->nq, pctx->wc->size); // v <- v^{q-2} (mod q)

	/* вычисляем z1 */
	ak_mpzn_mul_montgomery(z1, s, pctx->wc->r2q, pctx->wc->q, pctx->wc->nq, pctx->wc->size);
	ak_mpzn_mul_montgomery(z1, z1, v, pctx->wc->q, pctx->wc->nq, pctx->wc->size);
	ak_mpzn_mul_montgomery(z1, z1, pctx->wc->point.z, pctx->wc->q, pctx->wc->nq, pctx->wc->size);

	/* вычисляем z2 */
	ak_mpzn_mul_montgomery(z2, r, pctx->wc->r2q, pctx->wc->q, pctx->wc->nq, pctx->wc->size);
	ak_mpzn_sub(z2, pctx->wc->q, z2, pctx->wc->size);
	ak_mpzn_mul_montgomery(z2, z2, v, pctx->wc->q, pctx->wc->nq, pctx->wc->size);
	ak_mpzn_mul_montgomery(z2, z2, pctx->wc->point.z, pctx->wc->q, pctx->wc->nq, pctx->wc->size);
	struct epoint ecpoint, etpoint;
	/* сложение точек и проверка */
	ak_epoint_pow(&ecpoint, &ec->point, z1, pctx->wc->size, ec);
	ak_wpoint_to_epoint(&pctx->qpoint, pctx->wc, &etpoint, ec);
	ak_epoint_pow(&etpoint, &etpoint,z2, pctx->wc->size, ec);
	ak_epoint_add(&etpoint, &ecpoint, ec);
	ak_epoint_to_wpoint(&etpoint, ec, &cpoint, pctx->wc);
	ak_wpoint_reduce(&cpoint, pctx->wc);
	ak_mpzn_rem(cpoint.x, cpoint.x, pctx->wc->q, pctx->wc->size);

	if (ak_mpzn_cmp(cpoint.x, r, pctx->wc->size)) return ak_false;
	return ak_true; 
}
